#!/usr/bin/env python3
"""Proof-carrying search for a nonslice/ribbon equal-zero-surgery pair.

The public files are self-describing CSV tables.  A row is promoted only when
one side has a terminal smooth nonsliceness witness and the other side has a
replayable ribbon payload (a direct band to the unknot or a stored full ribbon
certificate).  Mere `slice=1`, `ribbon=1`, or `base_slice=1` flags are retained
as weak leads but never called terminal.
"""
from __future__ import annotations

import ast
import bz2
import csv
import json
import os
import re
from collections import Counter
from pathlib import Path
from typing import Iterator

ROOT = Path(os.environ.get("PUBLIC_DATA_ROOT", "public-data/Data"))
OUT = Path(os.environ.get("COLLISION_OUTPUT", "collision-output"))
OUT.mkdir(parents=True, exist_ok=True)

ZERO_FILES = ["zero_friends.csv.bz2", "more_zero_friends.csv.bz2"]
BLANKS = {"", "0", "0.0", "none", "null", "false", "[]", "{}", "nan", "?"}


def open_text(path: Path):
    if path.suffix == ".bz2":
        return bz2.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def rows(path: Path) -> Iterator[dict[str, str]]:
    with open_text(path) as fh:
        yield from csv.DictReader(fh)


def norm(x: object) -> str:
    return "" if x is None else str(x).strip()


def integer(x: object) -> int | None:
    s = norm(x)
    if not s:
        return None
    try:
        return int(float(s))
    except Exception:
        return None


def numbers(x: object) -> list[float]:
    s = norm(x)
    if not s:
        return []
    return [float(y) for y in re.findall(r"[-+]?\d+(?:\.\d+)?", s)]


def invariant_nonzero(x: object) -> bool:
    return any(abs(y) > 1e-12 for y in numbers(x))


def parse_mapping(x: object) -> dict:
    s = norm(x)
    if not s or s == "{}":
        return {}
    try:
        value = ast.literal_eval(s)
        return value if isinstance(value, dict) else {}
    except Exception:
        return {}


def ribbon_payload(rec: dict[str, str], source: str) -> tuple[bool, list[str]]:
    """Return only proof-bearing ribbon payloads, never a bare status flag."""
    witnesses: list[str] = []
    mapping = parse_mapping(rec.get("ribbon_to"))
    if "unknot" in mapping:
        witnesses.append(f"{source}:ribbon_to_unknot={mapping['unknot']}")
    cert = norm(rec.get("ribbon_cert"))
    if cert:
        witnesses.append(f"{source}:ribbon_cert")
    return bool(witnesses), witnesses


def nonslice_payload(rec: dict[str, str], source: str) -> tuple[bool, list[str]]:
    witnesses: list[str] = []
    if integer(rec.get("slice")) == -1:
        witnesses.append(f"{source}:slice=-1")
    for field in ("tau", "nu", "epsilon", "s_0", "s_2", "s_3", "CLS_red", "bls_odd", "Sq1_sum"):
        value = rec.get(field)
        if invariant_nonzero(value):
            witnesses.append(f"{source}:{field}={norm(value)}")
    # Donaldson and HKL fields are also rigorous nonsliceness obstructions,
    # but they are recorded separately so the smooth-vs-topological provenance
    # remains visible.
    if integer(rec.get("bifac_obs")) == 1:
        witnesses.append(f"{source}:Donaldson_bifactorability")
    for field in ("HKL_basic", "HKL_fancy", "HKL_direct"):
        if norm(rec.get(field)) not in ("", "0", "None"):
            witnesses.append(f"{source}:{field}={norm(rec.get(field))}")
    return bool(witnesses), witnesses


def lead_status(rec: dict[str, str]) -> dict:
    return {
        "slice": integer(rec.get("slice")),
        "ribbon": integer(rec.get("ribbon")),
        "base_slice": integer(rec.get("base_slice")),
        "ribbon_by_two_bands": integer(rec.get("ribbon_by_two_bands")),
        "ribbon_by_three_bands": integer(rec.get("ribbon_by_three_bands")),
    }


def load_table(path: Path) -> dict[str, dict[str, str]]:
    ans: dict[str, dict[str, str]] = {}
    if path.exists():
        for rec in rows(path):
            name = norm(rec.get("name"))
            if name:
                ans[name] = rec
    return ans


def compact_row(rec: dict[str, str]) -> dict:
    keep = [
        "name", "base_knot", "tri", "PD_code", "num_cross", "core_len",
        "slice", "ribbon", "base_slice", "ribbon_to",
        "ribbon_by_two_bands", "ribbon_by_three_bands", "s_2", "s_3",
        "Sq1_odd", "bls_odd", "Sq1_sum",
    ]
    return {k: rec.get(k, "") for k in keep if k in rec}


def main() -> None:
    unknown = load_table(ROOT / "plausibly_unknown.csv")
    slice16 = load_table(ROOT / "plausibly_slice_16.csv")
    auxiliary = {**slice16, **unknown}

    cls_names = []
    unknown_counts = Counter()
    for name, rec in unknown.items():
        if invariant_nonzero(rec.get("CLS_red")):
            cls_names.append(name)
        for f in ("slice", "ribbon", "s_0", "s_2", "s_3", "CLS_red", "tau", "nu", "epsilon"):
            if invariant_nonzero(rec.get(f)):
                unknown_counts[f] += 1
    (OUT / "CLS_red_nonzero_names.txt").write_text("\n".join(sorted(cls_names)) + "\n")

    terminal: list[dict] = []
    weak: list[dict] = []
    table_stats: list[dict] = []
    names_with_friends: set[str] = set()
    toxic_friend_rows = 0
    replayable_friend_ribbons = 0

    for filename in ZERO_FILES:
        path = ROOT / filename
        if not path.exists():
            continue
        count = 0
        for rec in rows(path):
            count += 1
            friend = norm(rec.get("name"))
            base = norm(rec.get("base_knot"))
            if not friend or not base:
                continue
            names_with_friends.update((friend, base))

            friend_aux = auxiliary.get(friend, {})
            base_aux = auxiliary.get(base, {})

            ft0, fw0 = nonslice_payload(rec, f"{filename}:friend")
            ft1, fw1 = nonslice_payload(friend_aux, "auxiliary:friend")
            friend_toxic = ft0 or ft1
            friend_toxic_w = fw0 + fw1
            if friend_toxic:
                toxic_friend_rows += 1

            bt, bw = nonslice_payload(base_aux, "auxiliary:base")
            # `base_slice=-1` is a terminal status in the friend table even if
            # the complete base record is not present in the 16-crossing file.
            if integer(rec.get("base_slice")) == -1:
                bt = True
                bw.append(f"{filename}:base_slice=-1")

            fr0, frw0 = ribbon_payload(rec, f"{filename}:friend")
            fr1, frw1 = ribbon_payload(friend_aux, "auxiliary:friend")
            friend_ribbon = fr0 or fr1
            friend_ribbon_w = frw0 + frw1
            if friend_ribbon:
                replayable_friend_ribbons += 1

            br, brw = ribbon_payload(base_aux, "auxiliary:base")

            common = {
                "friend_name": friend,
                "base_name": base,
                "friend_table": filename,
                "friend_exterior_tri": rec.get("tri", ""),
                "friend_core_len": rec.get("core_len", ""),
                "friend_PD_code": rec.get("PD_code", ""),
                "raw_friend_row": compact_row(rec),
            }

            if friend_toxic and br:
                terminal.append({
                    **common,
                    "toxic_name": friend,
                    "ribbon_name": base,
                    "toxic_side": "friend",
                    "nonslice_witnesses": friend_toxic_w,
                    "ribbon_witnesses": brw,
                    "zero_surgery_basis": "public zero-friend construction",
                })
            if bt and friend_ribbon:
                terminal.append({
                    **common,
                    "toxic_name": base,
                    "ribbon_name": friend,
                    "toxic_side": "base",
                    "nonslice_witnesses": bw,
                    "ribbon_witnesses": friend_ribbon_w,
                    "zero_surgery_basis": "public zero-friend construction",
                })

            # Weak leads preserve rows where database statuses suggest a
            # collision but a replayable ribbon payload is missing.
            friend_flags = lead_status(rec)
            if friend_toxic and integer(rec.get("base_slice")) == 1 and not br:
                weak.append({
                    **common,
                    "lead_type": "toxic_friend__base_slice_flag_without_public_ribbon_payload",
                    "nonslice_witnesses": friend_toxic_w,
                    "status": friend_flags,
                })
            if bt and (friend_flags["slice"] == 1 or friend_flags["ribbon"] == 1) and not friend_ribbon:
                weak.append({
                    **common,
                    "lead_type": "toxic_base__friend_slice_flag_without_replayable_payload",
                    "nonslice_witnesses": bw,
                    "status": friend_flags,
                })
            if friend_toxic and integer(rec.get("base_slice")) == 0:
                weak.append({
                    **common,
                    "lead_type": "toxic_friend__unresolved_base",
                    "nonslice_witnesses": friend_toxic_w,
                    "status": friend_flags,
                })

        table_stats.append({"file": filename, "rows": count})

    def dedup(records: list[dict], keys: tuple[str, ...]) -> list[dict]:
        out: dict[tuple, dict] = {}
        for r in records:
            out[tuple(r.get(k, "") for k in keys)] = r
        return sorted(out.values(), key=lambda r: tuple(r.get(k, "") for k in keys))

    terminal = dedup(terminal, ("toxic_name", "ribbon_name", "friend_exterior_tri"))
    weak = dedup(weak, ("lead_type", "friend_name", "base_name", "friend_exterior_tri"))

    summary = {
        "public_unknown_rows": len(unknown),
        "public_slice16_rows": len(slice16),
        "CLS_red_nonzero_count": len(cls_names),
        "CLS_red_matches_published_890_count": len(cls_names) == 890,
        "unknown_nonzero_field_counts": dict(unknown_counts),
        "friend_tables": table_stats,
        "distinct_names_with_friends": len(names_with_friends),
        "toxic_friend_rows": toxic_friend_rows,
        "replayable_friend_ribbon_rows": replayable_friend_ribbons,
        "terminal_collision_count": len(terminal),
        "weak_lead_count": len(weak),
        "determination": (
            "TERMINAL_RIBBON_FRIEND_COLLISION_FOUND" if terminal
            else "NO_TERMINAL_COLLISION_IN_PUBLIC_TABLES"
        ),
        "acceptance_rule": (
            "same row must join a terminal nonsliceness witness to a replayable "
            "ribbon payload; bare slice/ribbon flags are not terminal"
        ),
    }

    (OUT / "terminal_collisions.json").write_text(json.dumps(terminal, indent=2))
    (OUT / "weak_leads.json").write_text(json.dumps(weak, indent=2))
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))

    with (OUT / "terminal_collisions.csv").open("w", newline="", encoding="utf-8") as fh:
        fields = ["toxic_name", "ribbon_name", "toxic_side", "friend_table", "friend_exterior_tri", "nonslice_witnesses", "ribbon_witnesses"]
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for rec in terminal:
            row = {k: rec.get(k, "") for k in fields}
            row["nonslice_witnesses"] = "; ".join(rec["nonslice_witnesses"])
            row["ribbon_witnesses"] = "; ".join(rec["ribbon_witnesses"])
            w.writerow(row)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
