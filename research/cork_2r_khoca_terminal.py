#!/usr/bin/env python3
from __future__ import annotations

import json
import time
from pathlib import Path

from khoca import InteractiveCalculator

CANDIDATES = {
    "positive_trefoil": [[4, 2, 5, 1], [6, 4, 1, 3], [2, 6, 3, 5]],
    "figure_eight": [[4, 2, 5, 1], [2, 7, 3, 8], [8, 3, 1, 4], [6, 5, 7, 6]],
    "K14n3411_exact_cocore": [
        [20, 26, 21, 25], [19, 10, 20, 11], [17, 5, 18, 4],
        [13, 7, 14, 6], [1, 14, 2, 15], [15, 0, 16, 1],
        [9, 22, 10, 23], [8, 27, 9, 0], [12, 3, 13, 4],
        [26, 22, 27, 21], [2, 7, 3, 8], [5, 17, 6, 16],
        [24, 11, 25, 12], [23, 19, 24, 18],
    ],
    "m_n2_b_p1_T2_3_keep1_exact_cocore": [
        [29, 48, 30, 49], [37, 33, 38, 32], [38, 21, 39, 22],
        [41, 13, 42, 12], [16, 31, 17, 32], [28, 4, 29, 3],
        [0, 25, 1, 26], [15, 23, 16, 22], [24, 1, 25, 2],
        [46, 8, 47, 7], [20, 14, 21, 13], [19, 40, 20, 41],
        [33, 45, 34, 44], [34, 9, 35, 10], [11, 43, 12, 42],
        [26, 6, 27, 5], [2, 0, 3, 49], [8, 46, 9, 45],
        [47, 37, 48, 36], [14, 40, 15, 39], [6, 35, 7, 36],
        [30, 17, 31, 18], [23, 19, 24, 18], [43, 11, 44, 10],
        [4, 28, 5, 27],
    ],
}


def extract_s(result):
    active_q = []
    for t, q, torsion, coefficient in result[1]:
        coefficient = int(coefficient)
        if coefficient:
            active_q.extend([int(q)] * abs(coefficient))
    active_q.sort()
    if len(active_q) != 2 or active_q[1] - active_q[0] != 2:
        return {"ok": False, "active_q": active_q, "s": None}
    return {
        "ok": True,
        "active_q": active_q,
        "s": (active_q[0] + active_q[1]) // 2,
    }


def run(field, pd):
    started = time.time()
    calculator = InteractiveCalculator(field, (0, -1), 0)
    result, messages = calculator(pd, print_messages=True, progress=False)
    return {
        "field": field,
        "elapsed_seconds": time.time() - started,
        "extracted": extract_s(result),
        "result": result,
        "messages": messages,
    }


def main():
    output = {
        "schema": "cork-amplified-2r-khoca-terminal-v2",
        "frobenius_algebra": "F[X]/(X^2-X)",
        "candidates": {},
    }
    for name, pd in CANDIDATES.items():
        record = {"crossings": len(pd), "pd": pd}
        for label, field in (("Q", 1), ("F3", 3)):
            print("START", name, label, flush=True)
            record[label] = run(field, pd)
            print("DONE", name, label, record[label]["extracted"], flush=True)
        output["candidates"][name] = record
        Path("cork-2r-khoca-result.json").write_text(
            json.dumps(output, indent=2, sort_keys=True)
        )

    controls = output["candidates"]
    assert controls["positive_trefoil"]["Q"]["extracted"]["s"] in (2, -2)
    assert controls["figure_eight"]["Q"]["extracted"]["s"] == 0
    output["terminal"] = {
        "calibrated": True,
        "targets": {},
        "nonzero_ordinary_rasmussen_targets": [],
    }
    for name in CANDIDATES:
        if name in ("positive_trefoil", "figure_eight"):
            continue
        target = controls[name]
        target_result = {
            "s_Q": target["Q"]["extracted"]["s"],
            "s_F3": target["F3"]["extracted"]["s"],
            "nonzero_ordinary_rasmussen": target["Q"]["extracted"]["s"] not in (None, 0),
        }
        output["terminal"]["targets"][name] = target_result
        if target_result["nonzero_ordinary_rasmussen"]:
            output["terminal"]["nonzero_ordinary_rasmussen_targets"].append(name)
    Path("cork-2r-khoca-result.json").write_text(
        json.dumps(output, indent=2, sort_keys=True)
    )
    print(json.dumps(output["terminal"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
