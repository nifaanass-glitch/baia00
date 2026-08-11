#!/usr/bin/env python3
"""Search entire 0-surgery friend clusters for a toxic/ribbon collision.

Each row in zero_friends.csv describes a knot (by triangulation and PD code)
whose oriented 0-surgery agrees with that of `base_knot`.  Therefore two
*different* rows with the same base already have the same oriented 0-surgery.
This scanner looks for one row with an explicit smooth nonsliceness invariant
and another row with a replayable ribbon-to-unknot payload.

Fail-closed rules:
  * bare slice/ribbon flags are leads, never terminal certificates;
  * terminal nonsliceness requires a displayed invariant value;
  * terminal ribbonness requires `ribbon_to` containing `unknot` or a complete
    stored ribbon certificate;
  * a table-level collision still requires independent replay of the oriented
    0-surgery and ribbon certificates before publication.
"""
from __future__ import annotations

import ast
import bz2
import csv
import hashlib
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterator

ROOT = Path(os.environ.get("PUBLIC_DATA_ROOT", "public-data/Data"))
OUT = Path(os.environ.get("COLLISION_OUTPUT", "collision-output"))
OUT.mkdir(parents=True, exist_ok=True)
ZERO_FILES = ["zero_friends.csv.bz2", "more_zero_friends.csv.bz2"]


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
    try:
        return int(float(norm(x)))
    except Exception:
        return None


def numbers(x: object) -> list[float]:
    return [float(y) for y in re.findall(r"[-+]?\d+(?:\.\d+)?", norm(x))]


def nonzero(x: object) -> bool:
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


def load_named(path: Path) -> dict[str, dict[str, str]]:
    ans: dict[str, dict[str, str]] = {}
    if path.exists():
        for rec in rows(path):
            name = norm(rec.get("name"))
            if name:
                ans[name] = rec
    return ans


def invariant_witnesses(rec: dict[str, str], source: str) -> list[str]:
    witnesses: list[str] = []
    # Explicit smooth concordance/Khovanov witnesses only.  A bare `slice=-1`
    # is retained separately as a provenance lead.
    for field in (
        "tau", "nu", "epsilon", "s_0", "s_2", "s_3", "CLS_red",
        "bls_odd", "Sq1_sum",
    ):
        value = rec.get(field)
        if nonzero(value):
            witnesses.append(f"{source}:{field}={norm(value)}")
    if integer(rec.get("bifac_obs")) == 1:
        witnesses.append(f"{source}:Donaldson_bifactorability")
    for field in ("HKL_basic", "HKL_fancy", "HKL_direct"):
        value = norm(rec.get(field))
        if value not in ("", "0", "None", "[]", "{}"):
            witnesses.append(f"{source}:{field}={value}")
    return witnesses


def ribbon_witnesses(rec: dict[str, str], source: str) -> list[str]:
    witnesses: list[str] = []
    mapping = parse_mapping(rec.get("ribbon_to"))
    if "unknot" in mapping:
        witnesses.append(f"{source}:ribbon_to_unknot={mapping['unknot']}")
    cert = norm(rec.get("ribbon_cert"))
    if cert:
        witnesses.append(f"{source}:complete_ribbon_cert")
    return witnesses


def friend_id(rec: dict[str, str], filename: str) -> str:
    tri = norm(rec.get("tri"))
    pd = norm(rec.get("PD_code"))
    digest = hashlib.sha256(pd.encode()).hexdigest()[:16]
    return f"{filename}:{tri}:{digest}"


def compact_friend(rec: dict[str, str], filename: str) -> dict:
    return {
        "id": friend_id(rec, filename),
        "source_table": filename,
        "base_knot": norm(rec.get("base_knot")),
        "tri": norm(rec.get("tri")),
        "PD_code": norm(rec.get("PD_code")),
        "num_cross": norm(rec.get("num_cross")),
        "core_len": norm(rec.get("core_len")),
        "slice_flag": integer(rec.get("slice")),
        "ribbon_flag": integer(rec.get("ribbon")),
        "base_slice_flag": integer(rec.get("base_slice")),
        "ribbon_to": norm(rec.get("ribbon_to")),
        "ribbon_by_two_bands": integer(rec.get("ribbon_by_two_bands")),
        "ribbon_by_three_bands": integer(rec.get("ribbon_by_three_bands")),
        "s_2": norm(rec.get("s_2")),
        "s_3": norm(rec.get("s_3")),
        "bls_odd": norm(rec.get("bls_odd")),
        "Sq1_sum": norm(rec.get("Sq1_sum")),
    }


def dedup(records: list[dict], keys: tuple[str, ...]) -> list[dict]:
    out: dict[tuple, dict] = {}
    for rec in records:
        out[tuple(rec.get(k, "") for k in keys)] = rec
    return sorted(out.values(), key=lambda r: tuple(r.get(k, "") for k in keys))


def main() -> None:
    unknown = load_named(ROOT / "plausibly_unknown.csv")
    slice16 = load_named(ROOT / "plausibly_slice_16.csv")
    auxiliary = {**slice16, **unknown}

    clusters: dict[str, list[dict]] = defaultdict(list)
    table_stats: list[dict] = []
    status_counts = Counter()

    for filename in ZERO_FILES:
        path = ROOT / filename
        if not path.exists():
            continue
        count = 0
        for raw in rows(path):
            count += 1
            base = norm(raw.get("base_knot"))
            if not base:
                continue
            friend = compact_friend(raw, filename)
            friend["nonslice_witnesses"] = invariant_witnesses(raw, f"{filename}:friend")
            friend["ribbon_witnesses"] = ribbon_witnesses(raw, f"{filename}:friend")
            clusters[base].append(friend)
            if friend["nonslice_witnesses"]:
                status_counts["toxic_friend_rows"] += 1
            if friend["ribbon_witnesses"]:
                status_counts["replayable_ribbon_friend_rows"] += 1
            if friend["slice_flag"] == -1 and not friend["nonslice_witnesses"]:
                status_counts["nonslice_status_without_explicit_witness"] += 1
        table_stats.append({"file": filename, "rows": count})

    terminal: list[dict] = []
    provenance_leads: list[dict] = []
    unresolved_leads: list[dict] = []

    for base, friends in clusters.items():
        toxic_friends = [f for f in friends if f["nonslice_witnesses"]]
        ribbon_friends = [f for f in friends if f["ribbon_witnesses"]]

        # Strongest join: two distinct friends of one base.
        for toxic in toxic_friends:
            for ribbon in ribbon_friends:
                if toxic["id"] == ribbon["id"]:
                    continue
                terminal.append({
                    "collision_type": "friend_friend_same_base",
                    "base_knot": base,
                    "toxic_object": toxic,
                    "ribbon_object": ribbon,
                    "nonslice_witnesses": toxic["nonslice_witnesses"],
                    "ribbon_witnesses": ribbon["ribbon_witnesses"],
                    "zero_surgery_basis": (
                        "both rows are oriented 0-friends of the same base_knot"
                    ),
                })

        base_rec = auxiliary.get(base, {})
        base_nonslice = invariant_witnesses(base_rec, "auxiliary:base")
        base_ribbon = ribbon_witnesses(base_rec, "auxiliary:base")

        # Toxic named base + ribbon friend.
        if base_nonslice:
            for ribbon in ribbon_friends:
                terminal.append({
                    "collision_type": "base_friend",
                    "base_knot": base,
                    "toxic_object": {"name": base, "auxiliary_record": base_rec},
                    "ribbon_object": ribbon,
                    "nonslice_witnesses": base_nonslice,
                    "ribbon_witnesses": ribbon["ribbon_witnesses"],
                    "zero_surgery_basis": "ribbon row is an oriented 0-friend of base_knot",
                })

        # Ribbon named base + toxic friend.
        if base_ribbon:
            for toxic in toxic_friends:
                terminal.append({
                    "collision_type": "friend_base",
                    "base_knot": base,
                    "toxic_object": toxic,
                    "ribbon_object": {"name": base, "auxiliary_record": base_rec},
                    "nonslice_witnesses": toxic["nonslice_witnesses"],
                    "ribbon_witnesses": base_ribbon,
                    "zero_surgery_basis": "toxic row is an oriented 0-friend of base_knot",
                })

        # Preserve rows that could close after recovering missing provenance.
        if toxic_friends and not ribbon_friends:
            direct_flagged = [
                f for f in friends
                if f["slice_flag"] == 1 or f["ribbon_flag"] == 1
                or f["ribbon_by_two_bands"] == 1 or f["ribbon_by_three_bands"] == 1
            ]
            if direct_flagged:
                provenance_leads.append({
                    "base_knot": base,
                    "lead_type": "toxic_friend_plus_ribbon_flag_without_replayable_payload",
                    "toxic_friends": toxic_friends,
                    "flagged_friends": direct_flagged,
                })
            elif any(f["base_slice_flag"] == 0 for f in toxic_friends):
                unresolved_leads.append({
                    "base_knot": base,
                    "lead_type": "toxic_friend_cluster_with_unresolved_base_and_no_ribbon_friend",
                    "toxic_friends": toxic_friends,
                    "cluster_size": len(friends),
                })

        if ribbon_friends and not toxic_friends:
            status_toxic = [f for f in friends if f["slice_flag"] == -1]
            if status_toxic:
                provenance_leads.append({
                    "base_knot": base,
                    "lead_type": "ribbon_friend_plus_nonslice_status_without_explicit_witness",
                    "ribbon_friends": ribbon_friends,
                    "status_toxic_friends": status_toxic,
                })

    terminal = dedup(
        terminal,
        ("collision_type", "base_knot", "nonslice_witnesses", "ribbon_witnesses"),
    )
    provenance_leads = dedup(provenance_leads, ("base_knot", "lead_type"))
    unresolved_leads = dedup(unresolved_leads, ("base_knot", "lead_type"))

    cls_names = [name for name, rec in unknown.items() if nonzero(rec.get("CLS_red"))]
    (OUT / "CLS_red_nonzero_names.txt").write_text("\n".join(sorted(cls_names)) + "\n")
    (OUT / "terminal_cluster_collisions.json").write_text(json.dumps(terminal, indent=2))
    (OUT / "provenance_leads.json").write_text(json.dumps(provenance_leads, indent=2))
    (OUT / "unresolved_toxic_clusters.json").write_text(json.dumps(unresolved_leads, indent=2))

    summary = {
        "public_unknown_rows": len(unknown),
        "public_slice16_rows": len(slice16),
        "friend_tables": table_stats,
        "base_clusters": len(clusters),
        "friend_rows_total": sum(len(v) for v in clusters.values()),
        "status_counts": dict(status_counts),
        "terminal_cluster_collision_count": len(terminal),
        "provenance_lead_count": len(provenance_leads),
        "unresolved_toxic_cluster_count": len(unresolved_leads),
        "CLS_red_nonzero_count_in_old_public_unknown_table": len(cls_names),
        "determination": (
            "TERMINAL_TABLE_LEVEL_COLLISION_FOUND" if terminal
            else "NO_TERMINAL_COLLISION_IN_PUBLIC_FRIEND_CLUSTERS"
        ),
        "next_validation_if_hit": [
            "replay the ribbon movie",
            "certify the orientation-preserving zero-surgery homeomorphism",
            "replay the nonsliceness invariant",
            "assemble the trace-plus-disk-exterior homotopy sphere",
        ],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
