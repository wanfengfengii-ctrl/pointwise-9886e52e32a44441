import json
import sys
from design import Design
from source import load, completions, baseline, summarize


def normal():
    p = {"n": 2, "block": 2, "on_per_block": 1, "hot": [0, 0],
         "hot_on": 0, "max_run": 3, "pre_y": [0, 0, 0], "threshold": 1}
    d = Design(p)
    paths = list(d.tiny_paths())
    assert paths == [(0, 1), (1, 0)]
    assert [d.path_mass(a, [0, 0], 0, 0) for a in paths] == [6, 3]
    assert sum(d.path_mass(a, [0, 0], 0, 0) for a in paths) == 3 ** 2
    hot_history = Design(dict(p, pre_y=[9, 9, 9]))
    assert [hot_history.path_mass(a, [0, 0], 0, 0) for a in paths] == [3, 6]
    p4 = dict(p, n=4, block=4, on_per_block=2, hot=[0, 0, 0, 0])
    d4 = Design(p4)
    paths4 = list(d4.tiny_paths())
    low = [d4.path_mass(a, [0]*4, 0, 0) for a in paths4]
    high = [d4.path_mass(a, [9]*4, 0, 0) for a in paths4]
    assert sum(low) == sum(high) == 3**4 and low != high
    a = load()
    d = Design(a["protocol"])
    assert d.possible(d.initial)
    cases = list(completions(a))
    assert len(cases) == 27
    assert len({tuple(c["y"]) for c in cases}) == 27
    for c in cases:
        assert d.path_mass(a["on"], baseline(a, c["y"], 0, 0), 0, 0) > 0
    assert a["latest"][9]["rev"] == 3
    assert a["latest"][9]["value"][0] == 20
    assert a["latest"][18]["rev"] == 2
    print("NORMAL_REGRESSION_OK")
    print(json.dumps(summarize(a), ensure_ascii=False, sort_keys=True))


def demo():
    a = load()
    print(json.dumps(summarize(a), ensure_ascii=False, indent=2))
    print("highest revision at period 9:", a["latest"][9])
    for c in list(completions(a))[:2]:
        print(json.dumps({"id": c["id"], "missing_y": {str(i): c["y"][i-1]
              for i, row in a["latest"].items() if row["value"][0] is None}}, sort_keys=True))
    print("No treatment-effect inference is implemented in this archive reader.")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "normal"
    if mode == "normal":
        normal()
    elif mode == "demo":
        demo()
    else:
        raise SystemExit("use normal or demo")
