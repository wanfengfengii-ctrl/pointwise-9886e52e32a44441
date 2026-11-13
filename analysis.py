"""Exact AER-48 review analysis: registered adaptive schedule, integer path mass,
joint meter completions, exact two-sided p-values, joint 95% confidence set.

Usage from the repository root:
    python3 analysis.py run      # recompute result.json and REPORT.md
    python3 analysis.py verify   # independent recomputation + tiny exhaustive checks
    python3 analysis.py demo     # print the actual joint set, support and decision
"""
import hashlib
import json
import math
import os
import resource
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

from design import Design
from source import load, completions, baseline, summarize

ROOT = Path(__file__).resolve().parent
RESULT = ROOT / "result.json"
REPORT = ROOT / "REPORT.md"
ALPHA_NUM, ALPHA_DEN = 1, 20  # registered inclusion rule: p_max > 1/20


# ---------------------------------------------------------------- path-mass DP
def dp_distribution(design, y, on, tau, gamma, allowed_cache=None, mode="registered",
                    frozen=None):
    """Return {score: integer path mass} over all legal assignment paths.

    State: (last, run, used, hot_used, a[i-4], a[i-3], a[i-2], a[i-1]).
    The feasibility summary (last,run,used,hot_used) decides which next bits
    extend to a full legal path; the four previous bits decide the last three
    outcome values, hence the propensity; the score is accumulated separately
    and merged by value.  The score increment bit*b_i depends only on the
    observed y, the observed on vector and (tau,gamma), never on the path.
    """
    n = design.n
    pre = design.p["pre_y"]
    threshold = design.p["threshold"]
    if allowed_cache is None:
        allowed_cache = {}

    def allowed(state):
        value = allowed_cache.get(state)
        if value is None:
            value = design.allowed(state)
            allowed_cache[state] = value
        return value

    # key -> {score: mass}; initial last=-1, run=0, used=0, hot=0, four prior bits 0
    layer = {(-1, 0, 0, 0, 0, 0, 0, 0): {0: 1}}
    for i in range(n):
        nxt = {}
        for key, hist in layer.items():
            last, run, used, hot, c0, c1, c2, c3 = key
            state = (i, last, run, used, hot)
            window_sum = 0
            # j = i-3, i-2, i-1 with bit pairs (a_j, a_{j-1})
            for jj, (bit_j, bit_prev) in enumerate(((c1, c0), (c2, c1), (c3, c2))):
                j = i - 3 + jj
                if j < 0:
                    window_sum += pre[j]
                else:
                    window_sum += (y[j] + tau * (bit_j - on[j])
                                   + gamma * (bit_prev - (on[j - 1] if j >= 1 else 0)))
            choices = allowed(state)
            if mode == "uniform":
                weights = {bit: 1 for bit in choices if bit in (0, 1)}
            elif mode == "frozen":
                weights = {bit: frozen[i].get(bit, 0) for bit in choices}
                if not any(weights.values()):
                    continue
            elif len(choices) == 1:
                weights = {choices[0]: 3}
            elif len(choices) == 2:
                one = 2 if window_sum >= threshold else 1
                weights = {0: 3 - one, 1: one}
            else:
                continue
            for bit, weight in weights.items():
                ns = design.step(state, bit)
                if ns is None:
                    continue
                _, nl, nr, nu, nh = ns
                increment = bit * (y[i] - tau * on[i]
                                   - gamma * (on[i - 1] if i else 0))
                nk = (nl, nr, nu, nh, c1, c2, c3, bit)
                bucket = nxt.get(nk)
                if bucket is None:
                    bucket = {}
                    nxt[nk] = bucket
                for score, mass in hist.items():
                    ns_score = score + increment
                    bucket[ns_score] = bucket.get(ns_score, 0) + weight * mass
        layer = nxt

    total = {}
    for key, hist in layer.items():
        if key[2] == 0 and key[3] == design.hot_on:  # block counter reset, hot quota met
            for score, mass in hist.items():
                total[score] = total.get(score, 0) + mass
    return total


def count_legal_paths(design):
    """Number of complete legal assignments, from feasibility states only (no scores)."""
    from functools import lru_cache

    @lru_cache(None)
    def count(state):
        if state[0] == design.n:
            return 1 if (state[3] == 0 and state[4] == design.hot_on) else 0
        total = 0
        for bit in (0, 1):
            nxt = design.step(state, bit)
            if nxt is not None:
                total += count(nxt)
        return total

    return count(design.initial)


def distribution_digest(dist):
    payload = ";".join("%d:%d" % (s, dist[s]) for s in sorted(dist))
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def exact_two_sided(dist, s_obs, D):
    """Return (M, tail mass, reduced p numerator, reduced denominator)."""
    M = sum(s * mass for s, mass in dist.items())
    criterion = abs(s_obs * D - M)
    tail = sum(mass for s, mass in dist.items()
               if abs(s * D - M) >= criterion)  # equality belongs to the tail
    g = math.gcd(tail, D)
    return M, tail, tail // g, D // g


# ------------------------------------------------------------------- analysis
def _cell_task(args):
    ci, tau, gamma, y, on = args
    design = _WORKER["design"]
    mode = _WORKER.get("mode", "registered")
    frozen = _WORKER.get("frozen")
    if mode == "frozen" and frozen is None:
        frozen = observed_schedule(design, y, on)
    dist = dp_distribution(design, y, on, tau, gamma, _WORKER["allowed"],
                           mode=mode, frozen=frozen)
    total = sum(dist.values())
    D = _WORKER["D"]
    if total == 0:
        return (ci, tau, gamma, distribution_digest(dist), 0, 0, 0, 1, 1, 0)
    b = baseline(_WORKER["archive"], y, tau, gamma)
    s_obs = sum(on[i] * b[i] for i in range(len(on)))
    M = sum(s * m for s, m in dist.items())
    crit = abs(s_obs * total - M)
    tail = sum(m for s, m in dist.items() if abs(s * total - M) >= crit)
    g = math.gcd(tail, total)
    return (ci, tau, gamma, distribution_digest(dist), M, s_obs, tail,
            tail // g, total // g, total)


def _worker_init(archive, protocol, on, D, mode="registered", frozen=None):
    _WORKER["archive"] = archive
    _WORKER["design"] = Design(protocol)
    _WORKER["on"] = on
    _WORKER["D"] = D
    _WORKER["mode"] = mode
    _WORKER["frozen"] = frozen
    _WORKER["allowed"] = {}


_WORKER = {}


def compute_analysis(archive, parallel=True, evidence=True, mode="registered", frozen=None):
    """Exact cells, region, projection and decision; evidence is full and rerunnable."""
    protocol = archive["protocol"]
    design = Design(protocol)
    D = 3 ** protocol["n"]
    cases = list(completions(archive))
    taus = list(protocol["tau_grid"])
    gammas = list(protocol["gamma_grid"])
    if not cases:
        # conflicting materials: stop inference, emit no partial confidence set
        return {
            "cells": [], "region": [], "total_effect_support": [],
            "decision": "incompatible",
            "evidence": {"status": "NO_LEGAL_COMPLETION",
                         "message": "no integer completion satisfies all materials simultaneously; "
                                    "inference stopped, no partial confidence set emitted",
                         "conflicts": material_conflicts(ROOT),
                         "summaries": summarize(archive)},
        }

    tasks = [(ci, tau, gamma, c["y"], archive["on"])
             for ci, c in enumerate(cases) for tau in taus for gamma in gammas]
    results = []
    if parallel and len(tasks) >= 64 and os.cpu_count() and os.cpu_count() > 1:
        import multiprocessing as mp
        ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context()
        with ctx.Pool(processes=min(os.cpu_count(), 10), initializer=_worker_init,
                      initargs=(archive, protocol, archive["on"], D, mode, frozen)) as pool:
            results = pool.map(_cell_task, tasks, chunksize=8)
    else:
        _WORKER.update({"archive": archive, "design": design, "on": archive["on"],
                        "D": D, "allowed": {}, "mode": mode, "frozen": frozen})
        allowed_cache = {}
        for ci, tau, gamma, y, on in tasks:
            _WORKER["allowed"] = allowed_cache
            results.append(_cell_task((ci, tau, gamma, y, on)))

    p_by = {(r[1], r[2], r[0]): (r[7], r[8], r[4], r[5]) for r in results}
    digests = {(r[1], r[2], r[0]): (r[3], r[6], r[4], r[5]) for r in results}
    # registered mode must use exactly 3**n total mass
    if mode == "registered" and any(r[9] != D for r in results):
        raise AssertionError("PATH_MASS_NOT_NORMALIZED")
    cells = []
    per_cell_evidence = []
    for tau in taus:
        for gamma in gammas:
            entries = []
            for ci, c in enumerate(cases):
                num, den, M, s_obs = p_by[(tau, gamma, ci)]
                entries.append({"id": c["id"], "p": [num, den], "M": M,
                                "s_obs": s_obs,
                                "tail": digests[(tau, gamma, ci)][1],
                                "dist_sha256": digests[(tau, gamma, ci)][0]})
            worst = max(entries, key=lambda e: (Fraction(*e["p"]), -int(e["id"][1:])))
            worst_id = worst["id"]
            p_max = worst["p"]
            accepted = ALPHA_DEN * p_max[0] > ALPHA_NUM * p_max[1]
            cells.append({
                "tau": tau, "gamma": gamma,
                "p_by_completion": [e["p"] for e in entries],
                "worst_id": worst_id, "p_max": p_max, "accepted": accepted,
            })
            if evidence:
                per_cell_evidence.append({"tau": tau, "gamma": gamma,
                                          "completions": entries, "worst_id": worst_id,
                                          "p_max": p_max})
    region = [[c["tau"], c["gamma"]] for c in cells if c["accepted"]]
    support = sorted({c["tau"] + c["gamma"] for c in cells if c["accepted"]})
    if not region:
        decision = "incompatible"
    elif all(v < 0 for v in support):
        decision = "saving"
    else:
        decision = "inconclusive"
    return {
        "cells": cells, "region": region, "total_effect_support": support,
        "decision": decision,
        "evidence": {
            "path_mass_total": str(D), "n": protocol["n"],
            "alpha": [ALPHA_NUM, ALPHA_DEN],
            "include_rule": "accepted iff p_max > 1/20 (strict, exact integer comparison)",
            "protocol": protocol, "on": archive["on"],
            "n_completions": len(cases), "completion_ids": [c["id"] for c in cases],
            "completions": [{"id": c["id"], "y": c["y"]} for c in cases],
            "recovery_records": archive["totals"],
            "latest_revisions": {
                str(i): {"rev": archive["latest"][i]["rev"],
                         "line": archive["latest"][i]["line"],
                         "kwh": archive["latest"][i]["value"][0],
                         "range": list(archive["latest"][i]["value"][1]) or None}
                for i in sorted(archive["latest"])
            },
            "dp": {
                "state": "(i,last,run,used,hot_used,a[i-4],a[i-3],a[i-2],a[i-1])",
                "merge": "same state => same feasible continuations, same weights, same score increment; masses merged by score",
                "score_increment": "a_i * (y_i - tau*on_i - gamma*on_(i-1))",
                "terminal": "used==0 and hot_used==hot_on after period n",
                "mass_check": "sum H(s) == 3**n for every cell and completion",
            },
            "per_cell": per_cell_evidence,
            "runtime": runtime_stats(),
            "summaries": summarize(archive),
        },
    }


def observed_schedule(design, y, on):
    """Weights the scheduler used on the observed path (observed outcomes, fixed)."""
    history = list(design.p["pre_y"])
    state = design.initial
    schedule = []
    for i, bit in enumerate(on):
        weights = design.weights(state, history)
        schedule.append({0: weights.get(0, 0), 1: weights.get(1, 0)})
        history.append(y[i])
        history = history[-3:]
        state = design.step(state, bit)
    return schedule


def _variant_cells(archive, mode):
    """Run every (cell, completion) under a non-registered weight rule.

    Masses are normalised by the rule's own total path mass (uniform/frozen
    masses need not equal 3**n).
    """
    variant = compute_analysis(archive, evidence=False, mode=mode)
    p_max_all = [Fraction(*c["p_max"]) for c in variant["cells"]]
    return {"n_cells": len(variant["region"]), "region": variant["region"],
            "support": variant["total_effect_support"], "decision": variant["decision"],
            "max_p_at_most_alpha": sum(1 for p in p_max_all if p <= Fraction(ALPHA_NUM, ALPHA_DEN))}


def diagnostics(archive):
    """Non-registered comparisons used only to explain what each material feature does.

    They are clearly marked as invalid decompositions and never feed the decision.
    """
    notes = {}
    design = Design(archive["protocol"])
    notes["legal_paths"] = str(count_legal_paths(design))
    # stale lower-revision exact readings that the highest revision revokes:
    # take them from the raw export audit, never from the cleaned archive.
    records = scan_revisions(ROOT)
    audit = {rec["period"]: rec for rec in records if rec["status"] == "ACTIVE_BASE"}
    revoked = [rec for rec in records
               if rec["status"] == "ACTIVE_OVERRIDES"
               and rec["kwh"] is not None
               and audit.get(rec["period"], {}).get("kwh") is not None]
    if revoked:
        stale = dict(archive)
        stale["latest"] = {i: dict(row) for i, row in archive["latest"].items()}
        periods = {}
        for rec in revoked:
            base = audit[rec["period"]]
            stale["latest"][rec["period"]] = {"rev": base["rev"],
                                              "value": (base["kwh"], ()),
                                              "line": base["line"]}
            periods[str(rec["period"])] = base["kwh"]
        r = compute_analysis(stale, evidence=False)
        notes["stale_readings"] = {
            "periods": periods,
            "n_cells": len(r["region"]), "region": r["region"],
            "support": r["total_effect_support"], "decision": r["decision"],
            "note": "uses revoked readings; invalid under the revision rule; shown for impact only",
        }
    for mode in ("uniform", "frozen"):
        notes[mode + "_weights"] = _variant_cells(archive, mode)
    return notes


def material_conflicts(root):
    """Locate conflicting/unsatisfied material positions without guessing."""
    root = Path(root)
    notes = []
    try:
        protocol = json.loads((root / "data/protocol.json").read_text())
        lines = (root / "data/meter.jsonl").read_text().splitlines()
        latest = {}
        for line_no, line in enumerate(lines, 1):
            row = json.loads(line)
            i, rev = row.get("period"), row.get("rev")
            value = (row.get("kwh"), tuple(row.get("range", [])))
            old = latest.get(i)
            if old and old["rev"] == rev and old["value"] != value:
                notes.append({"path": "data/meter.jsonl", "line": line_no,
                              "period": i, "rev": rev,
                              "issue": "same-rev conflict with line %d" % old["line"]})
            elif not old or rev > old["rev"]:
                latest[i] = {"rev": rev, "value": value, "line": line_no}
        missing = sorted(set(range(1, protocol["n"] + 1)) - set(latest))
        if missing:
            notes.append({"path": "data/meter.jsonl", "issue": "missing periods",
                          "periods": missing})
    except Exception as exc:  # materials unreadable: report, do not infer
        notes.append({"path": "data/", "issue": "MATERIAL_UNREADABLE", "detail": str(exc)})
    if not notes:
        notes.append({"path": "data/meter.jsonl + data/recovery.json",
                      "issue": "ALL_REGISTER_CONSTRAINTS_JOINTLY_UNSATISFIABLE",
                      "detail": "no integer completion satisfies the ranges and every recovery total"})
    return notes


def runtime_stats():
    self_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child_rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    peak_mb = max(self_rss, child_rss) / 1024.0
    return {"peak_rss_mb": round(peak_mb, 1)}


# ----------------------------------------------------------------- fault model
def scan_revisions(root):
    """Read meter export order, mirroring source.load's revision rules exactly."""
    root = Path(root)
    protocol = json.loads((root / "data/protocol.json").read_text())
    records = []
    latest = {}
    for line_no, line in enumerate((root / "data/meter.jsonl").read_text().splitlines(), 1):
        row = json.loads(line)
        i, rev = row["period"], row["rev"]
        value = (row.get("kwh"), tuple(row.get("range", [])))
        old = latest.get(i)
        if old and old["rev"] == rev and old["value"] == value:
            status = "DUPLICATE_SAME_REV"
        elif old and old["rev"] == rev:
            status = "REVISION_CONFLICT"
        elif old and rev < old["rev"]:
            status = "IGNORED_LATE_LOWER_REV"
        else:
            status = "ACTIVE_OVERRIDES" if old else "ACTIVE_BASE"
            latest[i] = {"rev": rev, "value": value, "line": line_no}
        records.append({"line": line_no, "period": i, "rev": rev,
                        "kwh": row.get("kwh"), "range": row.get("range"),
                        "status": status})
    return records


# ------------------------------------------------------------------ reporting
def frac(pair):
    return "%d/%d" % (pair[0], pair[1])


def build_report(results, runtime=None, archive=None):
    summary = results["evidence"].get("summaries") or {}
    cases = results["evidence"].get("completion_ids", [])
    cells = results["cells"]
    region = results["region"]
    support = results["total_effect_support"]
    decision = results["decision"]
    runtime = runtime or results["evidence"].get("runtime", {})
    lines = []
    a = lines.append

    a("# AER-48 曝气试验封存复核报告")
    a("")
    a("本报告只分析 `data/` 封存材料；所有数字由 `analysis.py` 精确整数递推重算。"
      "所谓“节能确认”仅是登记模型和登记置信规则下的结论，不是对现实机组因果结构的证明。")
    a("")
    a("## 0. 结论")
    a("")
    if decision == "incompatible":
        if not cells:
            a("材料约束无合法联合补全：**停止推断**，不输出部分置信集合（incompatible）。")
        else:
            a("联合 95% 置信集合为空：登记模型与本批材料不相容（incompatible），不作节能判定。")
    elif decision == "saving":
        a("联合 95% 置信集合非空，其 `tau+gamma` 精确投影全部为负："
          "在登记模型下可确认长期每期耗电下降（saving）。")
    else:
        a("联合 95% 置信集合非空，但 `tau+gamma` 投影并非全部为负："
          "不能确认节能（inconclusive）。")
    a("")
    a("| 项目 | 值 |")
    a("|---|---|")
    if summary:
        a("| 期数 n | %d |" % summary["periods"])
        a("| on 期数 | %d |" % summary["on_blocks"])
        a("| 不确定电表期 | %s |" % "、".join(str(i) for i in summary["uncertain_periods"]))
        a("| 恢复合计记录 | %d |" % summary["recovery_records"])
        a("| 联合补全数 | %d |" % summary["joint_completions"])
        a("| 联合置信格数 | %d / %d |" % (len(region), len(cells)))
        a("| tau+gamma 支持集 | %s |" % (support if support else "（空）"))
        a("| 决策 | **%s** |" % decision)
    a("")

    status = results["evidence"].get("status")
    if status in ("NO_LEGAL_COMPLETION", "MATERIAL_LOAD_FAILED"):
        a("## 材料冲突（停止推断）")
        a("")
        if status == "NO_LEGAL_COMPLETION":
            a("无任何整数联合补全同时满足区间与全部恢复合计。依据 `result.json` 的 "
              "`evidence.conflicts`，停止推断，不输出部分置信集合：")
        else:
            a("密封材料无法按登记规则一致加载。依据 `result.json` 的 "
              "`evidence.conflicts`，停止推断，不输出部分置信集合：")
        a("")
        for note in results["evidence"].get("conflicts", []):
            a("- `%s`：%s（%s）" % (note.get("path"), note.get("issue"),
                                   note.get("detail", "见结构化证据")))
        a("")
        return "\n".join(lines) + "\n"

    a("## 1. 材料、修订与恢复合计")
    a("")
    p = (archive or {}).get("protocol", {})
    a("- 设计格点：`data/protocol.json` 键 `tau_grid`、`gamma_grid`（各 %d 个值），"
      "`n=%d`、`block=%d`、`on_per_block=%d`、`max_run=%d`、`hot_on=%d`、`threshold=%d`、`pre_y=%s`。"
      % (len(p.get("tau_grid", [])), p.get("n", summary["periods"]), p.get("block", 0),
         p.get("on_per_block", 0), p.get("max_run", 0), p.get("hot_on", 0),
         p.get("threshold", 0), p.get("pre_y")))
    a("- 实际分配：`data/schedule.json` 键 `on`，共 %d 位、%d 个 on 期。"
      % (summary["periods"], summary["on_blocks"]))
    a("- 电表：`data/meter.jsonl` 按导出到达顺序逐行；每期仅最高 `rev` 有效。审计结果：")
    a("")
    a("| 行 | 期 | rev | kwh/range | 状态 |")
    a("|---|---|---|---|---|")
    records = scan_revisions(ROOT)
    for rec in records:
        if rec["status"] != "ACTIVE_BASE":
            val = rec["kwh"] if rec["kwh"] is not None else str(rec["range"])
            a("| %d | %d | %d | %s | %s |" % (rec["line"], rec["period"], rec["rev"], val, rec["status"]))
    a("")
    a("关键覆盖关系：")
    for rec in records:
        if rec["status"] in ("ACTIVE_OVERRIDES", "IGNORED_LATE_LOWER_REV"):
            same = [r for r in records if r["period"] == rec["period"]]
            base = [min(same, key=lambda r: r["rev"])] if same else []
            if rec["status"] == "IGNORED_LATE_LOWER_REV":
                a("- 第 %d 期：第 %d 行 rev %d（`kwh=%s`）晚到且版本更低，按规则忽略；"
                  "有效值取自最高修订（见 `evidence.latest_revisions`），**不得用旧读数填补**。"
                  % (rec["period"], rec["line"], rec["rev"], rec["kwh"]))
            elif base:
                a("- 第 %d 期：第 %d 行 rev %d（%s）为最高修订，覆盖第 %d 行 rev %d（%s）。"
                  % (rec["period"], rec["line"], rec["rev"],
                     rec["kwh"] if rec["kwh"] is not None else rec["range"],
                     base[0]["line"], base[0]["rev"],
                     base[0]["kwh"] if base[0]["kwh"] is not None else base[0]["range"]))
    ranges = [(rec["line"], rec["period"], rec["range"]) for rec in scan_revisions(ROOT)
              if rec["status"] == "ACTIVE_OVERRIDES" and rec["range"]]
    if ranges:
        a("- 最高修订为区间的期共 %d 个（行 %s），对应期 %s，范围 %s；这些期只能取区间内的整数。"
          % (len(ranges),
             "、".join(str(r[0]) for r in ranges),
             "、".join(str(r[1]) for r in ranges),
             "、".join(str(r[2]) for r in ranges)))
    totals = (archive or {}).get("totals") or []
    a("- 恢复合计：`data/recovery.json` 键 `totals`，%s，所有记录同时成立。" %
      "、".join("%s(%s)=%d" % (row["id"], ",".join(str(x) for x in row["periods"]), row["kwh"])
                for row in totals))
    a("")
    a("### 联合补全")
    a("")
    uncertain = summary["uncertain_periods"]
    product_size = 3 ** len(uncertain)
    a("%d 个不确定期的区间乘积为 3^%d=%d；同时满足全部 %d 条恢复合计的整数向量为 %d 个，"
      "按 `source.completions` 的字典序编号 %s–%s。补全只限制观测 `y`，"
      "**不筛选反事实分配**。" % (len(uncertain), len(uncertain), product_size,
                                  summary["recovery_records"],
                                  summary["joint_completions"],
                                  cases[0], cases[-1]))
    a("")

    a("## 2. 登记模型与精确检验")
    a("")
    a("登记恒定加性模型：`Y_i(a) = b_i + tau*a_i + gamma*a_(i-1)`，`a_0=0`，`b_i` 不随分配改变。"
      "对每个联合补全 `y` 和每对 `(tau,gamma)`，由实际 `on` 反推 `b`。")
    a("")
    n_periods = p.get("n") or summary.get("periods") or 0
    D = 3 ** n_periods
    a("对合法反事实分配 `a`，`S(a)=sum_i a_i*b_i`，路径概率为 "
      "`design.path_mass(a,b,tau,gamma)/3**n`（n=%d）。记 `D=3**%d`，`H(s)` 为得分 `s` 的整数路径质量和，"
      "`M=sum_s s*H(s)`，`s_obs=S(实际on)`。精确双侧 p 值：" % (n_periods, n_periods))
    a("")
    a("`p = sum(H(s) for |s*D-M| >= |s_obs*D-M|) / D`，等号计入尾部；全整数运算，无浮点阈值、无蒙特卡洛。")
    a("")
    a("每格的 `p_max` 是 %d 个联合补全 p 值的最大值；`p_max > 1/20`（严格）才纳入联合 95%% 置信集合；"
      "同分见证取最小补全 ID（按补全字典序编号）。" % summary["joint_completions"])
    a("")

    a("## 3. 路径质量递推与无错误合并")
    a("")
    a("递推状态 `(i,last,run,used,hot_used,a[i-4..i-1])`：")
    a("- `(i,last,run,used,hot_used)` 决定后续可行性（块配额、连续同值、hot 配额），即 `Design.possible` 的完整信息；")
    a("- 最近 4 个分配位与观测基线决定第 `i-3,i-2,i-1` 期的 `Y` 值，从而决定下一步倾向；基线 `b_i` 由观测 `y`、"
      "`on`、`tau`、`gamma` 固定，与反事实路径无关；")
    a("- 得分增量 `a_i*(y_i - tau*on_i - gamma*on_(i-1))` 与路径无关，递推对每个得分分别累计质量。")
    a("")
    a("因此同一状态的两条前缀，其全部后续合法延伸集合、每一步倾向和后续得分增量相同；"
      "按得分合并质量是等价合并，没有错误合并。终态仅保留 `used==0` 且 `hot_used==hot_on` 的完整合法路径。"
      "每次 `(格, 补全)` 都断言 `sum_s H(s)==3**%d`。复杂度为状态数×得分桶数，未枚举任何完整 %d 期分配。"
      % (n_periods, n_periods))
    a("")
    a("## 4. 逐格精确双侧 p 值、最大值与见证")
    a("")
    a("下表每行一个 `(tau,gamma)` 格，各列按补全 c0000–c%04d 顺序；`p_max` 为最大精确 p 值，"
      "`worst` 为其最小 ID 见证。" % (len(cases) - 1))
    a("")
    header = "| tau | gamma | " + " | ".join(cases) + " | p_max | worst | accepted |"
    a(header)
    a("|" + "---|" * (len(cases) + 5))
    for cell in cells:
        row = "| %d | %d | " % (cell["tau"], cell["gamma"])
        row += " | ".join(frac(p) for p in cell["p_by_completion"])
        row += " | %s | %s | %s |" % (frac(cell["p_max"]), cell["worst_id"],
                                       "yes" if cell["accepted"] else "no")
        a(row)
    a("")

    a("## 5. 联合 95% 置信集合")
    a("")
    if region:
        a("纳入格（按 (tau,gamma) 升序）：%s。" % "、".join("(%d,%d)" % (t, g) for t, g in region))
        comps = connectivity_components(region)
        a("")
        a("网格连通分量（共享边为连边）：%d 个。" % len(comps))
        for idx, comp in enumerate(comps, 1):
            a("- 分量 %d：%s" % (idx, "、".join("(%d,%d)" % (t, g) for t, g in comp)))
        holes = bounding_box_holes(region)
        if len(comps) > 1 or holes:
            a("")
            if len(comps) > 1:
                a("注意：该集合**不连通**；报告保留每一个格，不填洞、不取包络。")
            if holes:
                a("集合包围盒内被剔除的格（洞，不填）：%s。"
                  % "、".join("(%d,%d)" % h for h in holes))
            else:
                a("集合的包围盒内没有洞。")
    else:
        a("没有格满足 `p_max > 1/20`：联合置信集合为空。")
    a("")
    a("该集合是登记规则下的联合集合，不是把各参数独立区间拼起来的包络；"
      "判断只看 `tau+gamma` 的精确集合投影：")
    a("")
    a("`total_effect_support = {%s}`" % ", ".join(str(v) for v in support))
    a("")
    a("只有该投影非空且每个值都 < 0 才判为 saving；否则 inconclusive；集合为空则 incompatible。"
      "本次决策：**%s**。" % decision)
    a("")

    a("## 6. 补全极值完整性与最坏情形")
    a("")
    a("每个格对全部 %d 个补全逐一递推并取最大 p；没有任何补全被丢弃，也没有把补全概率加权或抽样。"
      "这是对“未知但联合可行的观测补全”的保守取最大（least favorable）处理。" % summary["joint_completions"])
    a("")
    accepted_cells = [c for c in cells if c["accepted"]]
    worst_cells = sorted(accepted_cells, key=lambda c: (Fraction(*c["p_max"]), c["tau"], c["gamma"]))[:5]
    a("联合集合中最接近纳入边界的 5 个格（p_max 最小但仍 > 1/20）：")
    a("")
    a("| tau | gamma | p_max | worst |")
    a("|---|---|---|---|")
    for c in worst_cells:
        a("| %d | %d | %s | %s |" % (c["tau"], c["gamma"], frac(c["p_max"]), c["worst_id"]))
    a("")
    a("每个 `p_max` 都能用 `result.json` 中同格的 `p_by_completion` 与 `worst_id` 直接复核；"
      "`evidence.per_cell` 另存每补全的 `s_obs`、`M`、尾质量和分布 SHA-256。")
    a("")

    a("## 7. 置信覆盖论证")
    a("")
    a("登记规则对每格取 `p_max`（%d 个补全中的最大精确双侧 p 值），仅当 `p_max > 1/20` 保留该格。"
      "这是保守的并集式覆盖：只要真实补全在这 %d 个联合可行向量中，其对应 p 值不超过最大值；"
      "若格被剔除，则每个补全都以 ≤1/20 的精确尾概率被拒绝。`p_max` 取最大意味着不依赖补全先验，"
      "也不做多重补全加权；联合集合的 95%% 名义水平由登记规则定义，不额外声称频率覆盖率。"
      "集合允许不连通，投影按集合精确取象，不以区间填洞。"
      % (summary["joint_completions"], summary["joint_completions"]))
    a("")

    a("## 7b. 各材料特征对结论的影响（仅解释，不参与决策）")
    a("")
    diag = results["evidence"].get("diagnostics") or {}
    if diag.get("status") == "SKIPPED":
        a("诊断计算跳过：%s。" % diag.get("detail"))
    else:
        a("- 完整合法分配数（仅按可行性计数）：%s。本分析从不枚举这些分配，只累计其路径质量。" %
          diag.get("legal_paths"))
        stale = diag.get("stale_readings") or {}
        a("- 若误用被撤销的旧读数（%s），联合集合为 %d 格，支持集 %s，决策 %s。"
          "该做法破坏修订规则，只用于说明“不得用旧读数填补”。"
          % ("；".join("第 %s 期用 %s" % (k, v) for k, v in stale.get("periods", {}).items()),
             stale.get("n_cells", 0), stale.get("support"), stale.get("decision")))
        a("- 滞后项 `gamma` 通过 `a_(i-1)` 进入反事实 `Y`，并改变下一步倾向；"
          "本批联合集合跨多个 `gamma`，说明滞后与同期效应在登记模型下不能分开确认。")
        for mode, label in (("uniform", "均匀置换（每步对可行位等权）"),
                            ("frozen", "冻结观测路径权重（不以模拟耗电重算）")):
            var = diag.get(mode + "_weights") or {}
            a("- 若改用%s：联合集合 %d 格，支持集 %s，决策 %s。"
              "这类替代不遵循登记的自适应倾向语义，只用于说明结论对分配概率模型的敏感性。"
              % (label, var.get("n_cells", 0), var.get("support"), var.get("decision")))
    a("")

    a("## 8. 运行实测与适用边界")
    a("")
    a("- 实测墙钟时间：%s 秒（`analysis.py run` 本次运行，含登记推断与第 7b 节诊断变体）。"
      % runtime.get("wall_seconds", "见 result.json"))
    a("- 实测峰值常驻内存：%s MB（`ru_maxrss`，取本进程与子进程最大值）。"
      % runtime.get("peak_rss_mb", "见 result.json"))
    a("- 适用边界：结论仅限本封存设计、登记恒定加性模型与登记置信规则；不声称模型符合现实，"
      "不外推到其他机组或其他调度。补全取最大 p 是保守处理，不表示每个补全同等可信。"
      "若 `data/meter.jsonl`、`data/recovery.json` 或 `data/schedule.json` 改动，必须重跑。")
    a("")

    a("## 9. 原始证据追溯")
    a("")
    a("- 材料：`data/protocol.json`（n、格点、hot、threshold、pre_y）、`data/schedule.json`（`on`）、"
      "`data/meter.jsonl` 第 1–57 行、`data/recovery.json`（`totals`）。")
    a("- 清洗：`source.load` 的最高修订规则、同版冲突拒绝、低版迟到忽略、区间整数与恢复合计联立；"
      "`source.completions` 的 %d 个联合补全。" % summary["joint_completions"])
    a("- 调度：`design.Design` 的 `step/possible/allowed/weights/path_mass`；"
      "`checks.py normal` 仍保留原有业务回归断言，`analysis.py verify` 会运行它。")
    a("- 推断：`result.json` 的 `cells`、`region`、`total_effect_support`、`decision` 与 `evidence.per_cell`"
      "（含每格每补全的精确 p、`M`、`s_obs`、尾质量及路径质量分布摘要）；`analysis.py verify` 重算并核对，"
      "另用 ≤10 期小设计的全路径枚举独立核对质量分布、尾部和补全极值。")
    a("")
    return "\n".join(lines) + "\n"


def adjacency(a, b):
    return abs(a[0] - b[0]) + abs(a[1] - b[1]) == 1


def bounding_box_holes(region):
    """Cells inside the (tau,gamma) bounding box that are not in the joint set."""
    taus = sorted({t for t, _ in region})
    gammas = sorted({g for _, g in region})
    box = {(t, g) for t in taus for g in gammas}
    return sorted(box - {tuple(x) for x in region})


def connectivity_components(region):
    remaining = {tuple(x) for x in region}
    comps = []
    while remaining:
        seed = min(remaining)
        stack = [seed]
        comp = set()
        while stack:
            cur = stack.pop()
            if cur not in remaining:
                continue
            remaining.remove(cur)
            comp.add(cur)
            for other in list(remaining):
                if adjacency(cur, other):
                    stack.append(other)
        comps.append(sorted(comp))
    return sorted(comps)


# ---------------------------------------------------------------------- entry
def produce():
    started = time.time()
    archive = load()
    results = compute_analysis(archive)
    try:
        results["evidence"]["diagnostics"] = diagnostics(archive)
    except Exception as exc:  # diagnostics never block the registered analysis
        results["evidence"]["diagnostics"] = {"status": "SKIPPED", "detail": str(exc)}
    results["evidence"]["runtime"] = {"wall_seconds": round(time.time() - started, 2),
                                      "peak_rss_mb": runtime_stats()["peak_rss_mb"]}
    return results


def write_outputs(results, archive=None):
    runtime = results["evidence"]["runtime"]
    report = build_report(results, runtime, archive)
    REPORT.write_text(report)
    RESULT.write_text(json.dumps(results, ensure_ascii=False, indent=1, sort_keys=True) + "\n")


def run():
    try:
        archive = load()
    except Exception as exc:
        # conflicting/unreadable materials: stop inference, do not emit a partial set
        try:
            notes = material_conflicts(ROOT)
        except Exception:
            notes = []
        notes = [n for n in notes
                 if n.get("issue") != "ALL_REGISTER_CONSTRAINTS_JOINTLY_UNSATISFIABLE"]
        conflicts = [{"path": "data/", "issue": "MATERIAL_LOAD_FAILED", "detail": str(exc)}] + notes
        results = {"cells": [], "region": [], "total_effect_support": [],
                   "decision": "incompatible",
                   "evidence": {"status": "MATERIAL_LOAD_FAILED", "conflicts": conflicts,
                                "runtime": {"wall_seconds": 0,
                                            "peak_rss_mb": runtime_stats()["peak_rss_mb"]}}}
        write_outputs(results, None)
        print("materials conflict: inference stopped; %s and %s record the conflict"
              % (RESULT.name, REPORT.name))
        raise SystemExit("MATERIAL_CONFLICT: %s" % exc)
    results = produce()
    write_outputs(results, archive)
    print("wrote %s and %s" % (RESULT.name, REPORT.name))
    print("decision=%s region=%d support=%s" % (results["decision"], len(results["region"]),
                                                results["total_effect_support"]))
    return results


def compare_json(a, b, path="$"):
    if type(a) is not type(b) and not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
        raise AssertionError("TYPE_MISMATCH %s: %r != %r" % (path, a, b))
    if isinstance(a, dict):
        if set(a) != set(b):
            raise AssertionError("KEY_MISMATCH %s: %s != %s" % (path, sorted(a), sorted(b)))
        for k in a:
            compare_json(a[k], b[k], path + "." + k)
    elif isinstance(a, list):
        if len(a) != len(b):
            raise AssertionError("LEN_MISMATCH %s" % path)
        for i, (x, y) in enumerate(zip(a, b)):
            compare_json(x, y, path + "[%d]" % i)
    elif isinstance(a, float) or isinstance(b, float):
        if abs(float(a) - float(b)) > 1e-9:
            raise AssertionError("VALUE_MISMATCH %s: %r != %r" % (path, a, b))
    elif a != b:
        raise AssertionError("VALUE_MISMATCH %s: %r != %r" % (path, a, b))


def strip_runtime(results):
    clone = json.loads(json.dumps(results))
    clone.get("evidence", {}).pop("runtime", None)
    return clone


def tiny_archives():
    """Small designs with missing ranges and recovery constraints (n <= 10)."""
    defs = [
        {"n": 8, "block": 4, "on_per_block": 2,
         "hot": [0, 1, 0, 0, 1, 0, 0, 0], "hot_on": 2, "max_run": 3,
         "pre_y": [20, 20, 20], "threshold": 60,
         "tau_grid": [-2, 0, 2], "gamma_grid": [-2, 0, 2],
         "exact": {1: 18, 2: 20, 4: 21, 5: 19, 7: 20, 8: 18},
         "ranges": {3: [19, 21], 6: [17, 19]},
         "totals": [{"id": "T1", "periods": [3, 6], "kwh": 39}]},
        {"n": 10, "block": 5, "on_per_block": 2,
         "hot": [0, 1, 0, 0, 0, 1, 0, 0, 0, 1], "hot_on": 3, "max_run": 3,
         "pre_y": [20, 20, 20], "threshold": 60,
         "tau_grid": [-1, 0, 1], "gamma_grid": [-1, 0, 1],
         "exact": {1: 18, 2: 20, 3: 21, 5: 19, 6: 18, 7: 20, 8: 21, 9: 19},
         "ranges": {4: [17, 19], 10: [20, 22]},
         "totals": [{"id": "T2", "periods": [4, 10], "kwh": 40}]},
    ]
    archives = []
    for spec in defs:
        protocol = {k: spec[k] for k in ("n", "block", "on_per_block", "hot", "hot_on",
                                         "max_run", "pre_y", "threshold", "tau_grid", "gamma_grid")}
        design = Design(protocol)
        if not design.possible(design.initial):
            raise AssertionError("TINY_DESIGN_INFEASIBLE")
        paths = list(design.tiny_paths())
        # actual on: a middle legal path (tiny_paths already enforce the hot quota)
        on = list(paths[len(paths) // 2])
        latest = {}
        for i in range(1, protocol["n"] + 1):
            if i in spec["ranges"]:
                lo, hi = spec["ranges"][i]
                latest[i] = {"rev": 1, "value": (None, (lo, hi)), "line": i}
            else:
                latest[i] = {"rev": 1, "value": (spec["exact"][i], ()), "line": i}
        archives.append({"protocol": protocol, "on": on, "latest": latest,
                         "totals": spec["totals"]})
    return archives


def enumerate_distribution(design, y, on, tau, gamma, paths):
    """Independent reference: enumerate every legal path and use design.path_mass."""
    base = baseline({"on": on}, y, tau, gamma)
    H = {}
    for a in paths:
        mass = design.path_mass(a, base, tau, gamma)
        s = sum(a[i] * base[i] for i in range(design.n))
        H[s] = H.get(s, 0) + mass
    return H


def tiny_verify():
    """Full enumeration check of path-mass distribution, tails and completion extrema."""
    checked = []
    for archive in tiny_archives():
        protocol = archive["protocol"]
        design = Design(protocol)
        D = 3 ** protocol["n"]
        paths = list(design.tiny_paths())
        assert sum(1 for _ in paths) == len(paths)
        cases = list(completions(archive))
        if not cases:
            raise AssertionError("TINY_NO_COMPLETION")
        if len(cases) < 2:
            raise AssertionError("TINY_COMPLETION_NOT_EXERCISED")
        # DP cells (forced serial inside this process)
        results = compute_analysis(archive, parallel=False)
        assert results["decision"] != "incompatible", "TINY_UNEXPECTED_DECISION"
        n_cells = len(protocol["tau_grid"]) * len(protocol["gamma_grid"])
        assert len(results["cells"]) == n_cells
        for cell, ev in zip(results["cells"], results["evidence"]["per_cell"]):
            tau, gamma = cell["tau"], cell["gamma"]
            direct = {}
            for ci, c in enumerate(cases):
                y = c["y"]
                dist = enumerate_distribution(design, y, archive["on"], tau, gamma, paths)
                assert sum(dist.values()) == D, "TINY_MASS"
                b = baseline(archive, y, tau, gamma)
                s_obs = sum(archive["on"][i] * b[i] for i in range(protocol["n"]))
                M, tail, num, den = exact_two_sided(dist, s_obs, D)
                # compare exact distribution, M, s_obs, tail and p with DP outputs
                dp_entry = ev["completions"][ci]
                assert dp_entry["dist_sha256"] == distribution_digest(dist), "TINY_DIST_DIFF"
                assert dp_entry["M"] == M and dp_entry["s_obs"] == s_obs, "TINY_MOMENT_DIFF"
                assert dp_entry["tail"] == tail, "TINY_TAIL_DIFF"
                assert cell["p_by_completion"][ci] == [num, den], "TINY_P_DIFF"
                direct[c["id"]] = (num / den, [num, den])
            # completion extrema: max p and smallest-ID witness recomputed directly
            best = max(direct.items(), key=lambda kv: (Fraction(*kv[1][1]), -int(kv[0][1:])))
            pm = best[1][1]
            assert cell["worst_id"] == best[0], "TINY_WORST_ID"
            assert cell["p_max"] == pm, "TINY_P_MAX"
            assert cell["accepted"] == (ALPHA_DEN * pm[0] > ALPHA_NUM * pm[1]), "TINY_ACCEPT"
        region = [[c["tau"], c["gamma"]] for c in results["cells"] if c["accepted"]]
        support = sorted({c["tau"] + c["gamma"] for c in results["cells"] if c["accepted"]})
        assert results["region"] == region and results["total_effect_support"] == support
        checked.append({"n": protocol["n"], "paths": len(paths),
                        "completions": len(cases), "cells": n_cells,
                        "tau_sums": sorted({c["tau"] + c["gamma"] for c in results["cells"]})})
    return checked


def verify():
    # 1) preserve and run the original business regression
    proc = subprocess.run([sys.executable, "checks.py", "normal"], cwd=str(ROOT),
                          capture_output=True, text=True)
    if proc.returncode != 0 or "NORMAL_REGRESSION_OK" not in proc.stdout:
        raise AssertionError("CHECKS_NORMAL_FAILED\n%s\n%s" % (proc.stdout, proc.stderr))
    print("checks.py normal: OK")

    # 2) independent tiny exhaustive verification (distribution, tail, completion extrema)
    tiny = tiny_verify()
    print("tiny exhaustive: OK", json.dumps(tiny, ensure_ascii=False, sort_keys=True))

    # 3) full recomputation from raw materials, exact structural comparison
    results = produce()
    stored = json.loads(RESULT.read_text())
    compare_json(strip_runtime(stored), strip_runtime(results))
    # evidence digests / tails / moments are part of the stored evidence and thus
    # compared above; additionally recompute each (cell, completion) independently.
    for cell in results["cells"]:
        entries = next(e for e in results["evidence"]["per_cell"]
                       if e["tau"] == cell["tau"] and e["gamma"] == cell["gamma"])
        for cid, entry, pair in zip(results["evidence"]["completion_ids"],
                                    entries["completions"], cell["p_by_completion"]):
            assert entry["id"] == cid and entry["p"] == pair, "VERIFY_EVIDENCE_ID"
    # direct DP recomputation of every (cell, completion) digest, tail, p and extremum
    archive = load()
    design = Design(archive["protocol"])
    cases = list(completions(archive))
    D = 3 ** archive["protocol"]["n"]
    ev_cells = {(e["tau"], e["gamma"]): e for e in results["evidence"]["per_cell"]}
    stored_max = {}
    for tau in archive["protocol"]["tau_grid"]:
        for gamma in archive["protocol"]["gamma_grid"]:
            best_p = None
            best_id = None
            cell_index = archive["protocol"]["tau_grid"].index(tau) * len(
                archive["protocol"]["gamma_grid"]) + archive["protocol"]["gamma_grid"].index(gamma)
            cell = results["cells"][cell_index]
            ev = ev_cells[(tau, gamma)]
            for ci, c in enumerate(cases):
                entry = ev["completions"][ci]
                dist = dp_distribution(design, c["y"], archive["on"], tau, gamma)
                assert distribution_digest(dist) == entry["dist_sha256"], "VERIFY_DIGEST"
                assert sum(dist.values()) == D, "VERIFY_MASS"
                b = baseline(archive, c["y"], tau, gamma)
                s_obs = sum(archive["on"][i] * b[i] for i in range(len(archive["on"])))
                M, tail, num, den = exact_two_sided(dist, s_obs, D)
                assert (M, s_obs, tail) == (entry["M"], entry["s_obs"], entry["tail"]), "VERIFY_TAIL"
                assert cell["p_by_completion"][ci] == [num, den], "VERIFY_P"
                assert entry["p"] == [num, den] and entry["id"] == c["id"], "VERIFY_EVIDENCE_ID"
                if best_p is None or Fraction(num, den) > best_p or (
                        Fraction(num, den) == best_p and int(c["id"][1:]) < int(best_id[1:])):
                    best_p, best_id = Fraction(num, den), c["id"]
            assert cell["p_max"] == [best_p.numerator, best_p.denominator], "VERIFY_P_MAX"
            assert cell["worst_id"] == best_id, "VERIFY_WORST"
            assert cell["accepted"] == (ALPHA_DEN * best_p.numerator > ALPHA_NUM * best_p.denominator), \
                "VERIFY_ACCEPT"
            stored_max[(tau, gamma)] = (best_id, cell["p_max"])
    # region/projection reconstructed from the recomputed per-cell maxima
    region = [[t, g] for (t, g), (_, pm) in sorted(stored_max.items())
              if ALPHA_DEN * pm[0] > ALPHA_NUM * pm[1]]
    assert results["region"] == region, "VERIFY_REGION"
    support = sorted({t + g for t, g in region})
    assert results["total_effect_support"] == support, "VERIFY_SUPPORT"
    print("full DP recomputation and evidence cross-check: OK")

    # 4) report must be exactly reproducible from the stored results/runtime
    runtime = stored["evidence"]["runtime"]
    expected = build_report(results, runtime, load())
    if REPORT.read_text() != expected:
        raise AssertionError("REPORT_MISMATCH")
    print("REPORT.md matches recomputed structured evidence: OK")

    # 5) decision consistency
    assert results["decision"] in ("saving", "inconclusive", "incompatible")
    if not results["region"]:
        assert results["decision"] == "incompatible"
    elif all(v < 0 for v in results["total_effect_support"]):
        assert results["decision"] == "saving"
    else:
        assert results["decision"] == "inconclusive"
    print("decision consistency: OK")
    print("verify: PASS")


def demo():
    results = produce()
    print("AER-48 joint set (%d cells):" % len(results["region"]))
    for tau, gamma in results["region"]:
        cell = next(c for c in results["cells"] if c["tau"] == tau and c["gamma"] == gamma)
        print("  (%d,%d)  p_max=%s  worst=%s" % (tau, gamma, frac(cell["p_max"]), cell["worst_id"]))
    print("total_effect_support (tau+gamma):", results["total_effect_support"])
    print("decision:", results["decision"])
    worst = sorted(results["cells"], key=lambda c: (-Fraction(*c["p_max"]), c["tau"], c["gamma"]))
    if worst:
        c = worst[0]
        print("global max p_max: (%d,%d) p_max=%s worst=%s" % (
            c["tau"], c["gamma"], frac(c["p_max"]), c["worst_id"]))
    accepted = [c for c in results["cells"] if c["accepted"]]
    if accepted:
        c = min(accepted, key=lambda c: (Fraction(*c["p_max"]), c["tau"], c["gamma"]))
        print("closest accepted cell to the 1/20 boundary: (%d,%d) p_max=%s worst=%s" % (
            c["tau"], c["gamma"], frac(c["p_max"]), c["worst_id"]))
    if results["decision"] == "saving":
        print("结论：登记模型下支持长期每期耗电下降（节能）")
    elif results["decision"] == "inconclusive":
        print("结论：不能确认节能（联合集合非空但总效应支持并非全负）")
    else:
        print("结论：模型与材料不相容，不作节能判定")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "run":
        run()
    elif mode == "verify":
        verify()
    elif mode == "demo":
        demo()
    else:
        raise SystemExit("use run, verify or demo")
