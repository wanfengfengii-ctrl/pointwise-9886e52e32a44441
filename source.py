"""Read the sealed archive and enumerate meter completions, not assignments."""
import json
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def load(root=ROOT):
    root = Path(root)
    protocol = json.loads((root / "data/protocol.json").read_text())
    schedule = json.loads((root / "data/schedule.json").read_text())
    latest = {}
    for line_no, line in enumerate((root / "data/meter.jsonl").read_text().splitlines(), 1):
        row = json.loads(line)
        i, rev = row["period"], row["rev"]
        if type(i) is not int or not 1 <= i <= protocol["n"] or type(rev) is not int or rev < 1:
            raise ValueError("METER_KEY")
        value = (row.get("kwh"), tuple(row.get("range", [])))
        if (value[0] is None) == (not value[1]):
            raise ValueError("METER_SHAPE")
        if value[0] is not None and (type(value[0]) is not int or value[0] < 0):
            raise ValueError("METER_VALUE")
        if value[1] and (len(value[1]) != 2 or any(type(x) is not int for x in value[1])
                         or not 0 <= value[1][0] <= value[1][1]):
            raise ValueError("METER_RANGE")
        old = latest.get(i)
        if old and old["rev"] == rev and old["value"] != value:
            raise ValueError("METER_REVISION_CONFLICT")
        if not old or rev > old["rev"]:
            latest[i] = {"rev": rev, "value": value, "line": line_no}
    if set(latest) != set(range(1, protocol["n"] + 1)):
        raise ValueError("METER_COVERAGE")
    if len(schedule["on"]) != protocol["n"] or any(type(x) is not int or x not in (0, 1) for x in schedule["on"]):
        raise ValueError("SCHEDULE_SHAPE")
    totals = json.loads((root / "data/recovery.json").read_text())["totals"]
    for row in totals:
        if not row["periods"] or any(i not in latest for i in row["periods"]):
            raise ValueError("RECOVERY_PERIOD")
    return {"protocol": protocol, "on": schedule["on"], "latest": latest, "totals": totals}


def completions(archive):
    """Each returned vector jointly satisfies all register recovery records."""
    n = archive["protocol"]["n"]
    domains = []
    for i in range(1, n + 1):
        exact, bounds = archive["latest"][i]["value"]
        domains.append([exact] if exact is not None else range(bounds[0], bounds[1] + 1))
    number = 0
    for values in product(*domains):
        if all(sum(values[i - 1] for i in row["periods"]) == row["kwh"]
               for row in archive["totals"]):
            yield {"id": "c%04d" % number, "y": list(values)}
            number += 1


def baseline(archive, y, tau, gamma):
    on = archive["on"]
    return [value - tau * on[i] - gamma * (on[i-1] if i else 0)
            for i, value in enumerate(y)]


def summarize(archive):
    uncertain = [i for i, row in archive["latest"].items() if row["value"][0] is None]
    return {"periods": archive["protocol"]["n"], "on_blocks": sum(archive["on"]),
            "uncertain_periods": sorted(uncertain),
            "recovery_records": len(archive["totals"]),
            "joint_completions": sum(1 for _ in completions(archive))}
