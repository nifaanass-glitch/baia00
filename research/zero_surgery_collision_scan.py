#!/usr/bin/env python3
"""Fail-closed search for a nonslice/ribbon same-zero-surgery collision.

The scanner consumes the public Dunfield--Gong data repository.  It never
promotes a hit merely from a name match: it records separately

  * the nonsliceness field and its raw value,
  * the ribbon/slice certificate field and its raw value,
  * the zero-friend row and orientation/signature metadata.

All parsing is streaming so the multi-million-row tables remain manageable.
"""
from __future__ import annotations

import bz2
import csv
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator

ROOT = Path(os.environ.get("PUBLIC_DATA_ROOT", "public-data/Data"))
OUT = Path(os.environ.get("COLLISION_OUTPUT", "collision-output"))
OUT.mkdir(parents=True, exist_ok=True)

PL_COLUMNS = [
    "name", "pd", "kind", "DT", "PD", "alexander", "s0", "s2", "s3",
    "thin", "LMO_degree1", "Knot_Floer", "known_slice", "CLS_red",
]
ZF_COLUMNS = [
    "base_name", "base_pd", "base_alexander", "base_s0", "base_s2",
    "base_s3", "base_LMO_degree1", "base_HFK", "red", "green", "blue",
    "framings", "friend_name", "friend_pd", "friend_alexander", "friend_s0",
    "friend_s2", "friend_s3", "friend_LMO_degree1", "friend_HFK", "iso_sig",
    "linking_matrix",
]

BLANKS = {"", "0", "0.0", "none", "null", "false", "[]", "{}", "nan", "?"}


def open_text(path: Path):
    if path.suffix == ".bz2":
        return bz2.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def first_existing(stems: Iterable[str]) -> Path | None:
    for stem in stems:
        p = ROOT / stem
        if p.exists():
            return p
    return None


def normalized(value: object) -> str:
    return "" if value is None else str(value).strip()


def meaningful(value: object) -> bool:
    return normalized(value).lower() not in BLANKS


def nonzero_numeric(value: object) -> bool:
    s = normalized(value)
    if s.lower() in BLANKS:
        return False
    nums = re.findall(r"(?<![A-Za-z])[-+]?\d+(?:\.\d+)?", s)
    return any(abs(float(x)) > 1e-12 for x in nums)


def hfk_nonzero(value: object) -> tuple[bool, str]:
    s = normalized(value)
    if s.lower() in BLANKS:
        return False, ""
    witnesses: list[str] = []
    for key in ("tau", "nu", "epsilon", "eps"):
        for match in re.finditer(rf"['\"]?{key}['\"]?\s*[:=]\s*([-+]?\d+)", s, re.I):
            if int(match.group(1)) != 0:
                witnesses.append(f"{key}={match.group(1)}")
    return bool(witnesses), ";".join(witnesses)


def detect_header(first: list[str], columns: list[str]) -> bool:
    a = [x.strip() for x in first]
    return a == columns or (a and a[0] == columns[0] and len(a) == len(columns))


def rows(path: Path, columns: list[str]) -> Iterator[dict[str, str]]:
    with open_text(path) as fh:
        reader = csv.reader(fh)
        try:
            first = next(reader)
        except StopIteration:
            return
        if detect_header(first, columns):
            header = [x.strip() for x in first]
        else:
            header = columns
            yield dict(zip(header, first))
        for row in reader:
            if not row:
                continue
            if len(row) < len(header):
                row = row + [""] * (len(header) - len(row))
            yield dict(zip(header, row[: len(header)]))


def toxic_from_side(rec: dict[str, str], prefix: str) -> tuple[bool, list[str]]:
    witnesses: list[str] = []
    for field in ("s0", "s2", "s3"):
        value = rec.get(f"{prefix}_{field}", "")
        if nonzero_numeric(value):
            witnesses.append(f"{field}={value}")
    hfk = rec.get(f"{prefix}_HFK", "")
    nz, why = hfk_nonzero(hfk)
    if nz:
        witnesses.append(f"HFK({why})")
    return bool(witnesses), witnesses


def ribbon_certificate(rec: dict[str, str]) -> tuple[bool, str]:
    raw = normalized(rec.get("known_slice"))
    kind = normalized(rec.get("kind"))
    # `known_slice` is the proof-bearing field in the public table.  We record
    # the raw payload and do not infer ribbonness from a bare database name.
    certified = meaningful(raw)
    return certified, raw or kind


def field_stats(path: Path, columns: list[str], max_examples: int = 5) -> dict:
    counts = {c: 0 for c in columns}
    nonzero = {c: 0 for c in columns}
    values: dict[str, Counter[str]] = {c: Counter() for c in columns}
    n = 0
    examples: list[dict[str, str]] = []
    for rec in rows(path, columns):
        n += 1
        if len(examples) < max_examples:
            examples.append(rec)
        for c in columns:
            v = normalized(rec.get(c))
            if meaningful(v):
                counts[c] += 1
                if len(values[c]) < 50 or v in values[c]:
                    values[c][v] += 1
            if nonzero_numeric(v):
                nonzero[c] += 1
    return {
        "path": str(path),
        "rows": n,
        "nonblank_counts": counts,
        "numeric_nonzero_counts": nonzero,
        "top_values": {c: values[c].most_common(20) for c in columns},
        "examples": examples,
    }


def main() -> None:
    file_inventory = [
        {"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size}
        for p in sorted(ROOT.rglob("*")) if p.is_file()
    ]
    (OUT / "file_inventory.json").write_text(json.dumps(file_inventory, indent=2))

    unknown_path = first_existing([
        "plausibly_unknown.csv", "plausibly_unknown.csv.bz2",
    ])
    slice_path = first_existing([
        "plausibly_slice.csv", "plausibly_slice.csv.bz2",
    ])
    zero_paths = [p for p in [
        first_existing(["zero_friends.csv", "zero_friends.csv.bz2"]),
        first_existing(["more_zero_friends.csv", "more_zero_friends.csv.bz2"]),
    ] if p is not None]

    summary: dict = {
        "root": str(ROOT),
        "unknown_path": str(unknown_path) if unknown_path else None,
        "slice_path": str(slice_path) if slice_path else None,
        "zero_friend_paths": [str(p) for p in zero_paths],
    }

    # 1. Inspect the mystery table and test whether CLS_red is exactly the
    #    new 890-knot obstruction field.
    leo_names: set[str] = set()
    unknown_by_name: dict[str, dict[str, str]] = {}
    if unknown_path:
        stats = field_stats(unknown_path, PL_COLUMNS)
        (OUT / "plausibly_unknown_stats.json").write_text(json.dumps(stats, indent=2))
        cls_nonzero = 0
        cls_nonblank = 0
        for rec in rows(unknown_path, PL_COLUMNS):
            name = normalized(rec.get("name"))
            if name:
                unknown_by_name[name] = rec
            raw = rec.get("CLS_red", "")
            if meaningful(raw):
                cls_nonblank += 1
            if nonzero_numeric(raw):
                cls_nonzero += 1
                if name:
                    leo_names.add(name)
        summary["CLS_red_nonblank"] = cls_nonblank
        summary["CLS_red_numeric_nonzero"] = cls_nonzero
        summary["CLS_red_matches_published_890_count"] = cls_nonzero == 890
        (OUT / "CLS_red_nonzero_names.txt").write_text("\n".join(sorted(leo_names)) + "\n")

    # 2. Read all friend pairs first.  Keep raw rows: 277k + 256k is small.
    friend_rows: list[dict[str, str]] = []
    all_friend_names: set[str] = set()
    zero_stats: list[dict] = []
    for path in zero_paths:
        count = 0
        for rec in rows(path, ZF_COLUMNS):
            count += 1
            friend_rows.append(rec)
            for key in ("base_name", "friend_name"):
                name = normalized(rec.get(key))
                if name:
                    all_friend_names.add(name)
        zero_stats.append({"path": str(path), "rows": count})
    summary["zero_friend_tables"] = zero_stats
    summary["zero_friend_pairs_total"] = len(friend_rows)
    summary["distinct_names_in_friend_tables"] = len(all_friend_names)

    # 3. Stream the 2.3M slice/ribbon table and retain only knots occurring in
    #    the friend tables.
    slice_hits: dict[str, dict[str, str]] = {}
    slice_rows_total = 0
    known_slice_payloads = Counter()
    if slice_path:
        for rec in rows(slice_path, PL_COLUMNS):
            slice_rows_total += 1
            name = normalized(rec.get("name"))
            cert, payload = ribbon_certificate(rec)
            if cert:
                known_slice_payloads[payload[:120]] += 1
            if name in all_friend_names and cert:
                slice_hits[name] = rec
    summary["slice_rows_total"] = slice_rows_total
    summary["friend_names_with_slice_certificate"] = len(slice_hits)
    summary["known_slice_payload_examples"] = known_slice_payloads.most_common(20)

    # 4. Join: existing standard invariants and, separately, the CLS_red field.
    hits: list[dict] = []
    for rec in friend_rows:
        for toxic_side, ribbon_side in (("base", "friend"), ("friend", "base")):
            toxic_name = normalized(rec.get(f"{toxic_side}_name"))
            ribbon_name = normalized(rec.get(f"{ribbon_side}_name"))
            if not toxic_name or ribbon_name not in slice_hits:
                continue

            standard_toxic, standard_witness = toxic_from_side(rec, toxic_side)
            cls_rec = unknown_by_name.get(toxic_name, {})
            cls_raw = normalized(cls_rec.get("CLS_red"))
            cls_toxic = nonzero_numeric(cls_raw)
            if not standard_toxic and not cls_toxic:
                continue

            ribbon_rec = slice_hits[ribbon_name]
            hit = {
                "toxic_name": toxic_name,
                "ribbon_name": ribbon_name,
                "toxic_side": toxic_side,
                "standard_nonslice": standard_toxic,
                "standard_witnesses": standard_witness,
                "CLS_red_nonslice": cls_toxic,
                "CLS_red_raw": cls_raw,
                "ribbon_known_slice_raw": ribbon_rec.get("known_slice", ""),
                "ribbon_kind": ribbon_rec.get("kind", ""),
                "zero_surgery_iso_sig": rec.get("iso_sig", ""),
                "framings": rec.get("framings", ""),
                "linking_matrix": rec.get("linking_matrix", ""),
                "raw_zero_friend_row": rec,
            }
            hits.append(hit)

    # Deduplicate by oriented name pair and signature.
    dedup: dict[tuple[str, str, str], dict] = {}
    for hit in hits:
        key = (hit["toxic_name"], hit["ribbon_name"], hit["zero_surgery_iso_sig"])
        dedup[key] = hit
    hits = sorted(dedup.values(), key=lambda x: (x["toxic_name"], x["ribbon_name"]))

    terminal_standard = [h for h in hits if h["standard_nonslice"]]
    tentative_cls = [h for h in hits if h["CLS_red_nonslice"]]
    summary.update({
        "joined_hits_total": len(hits),
        "standard_invariant_plus_slice_certificate_hits": len(terminal_standard),
        "CLS_red_plus_slice_certificate_hits": len(tentative_cls),
        "determination": (
            "STANDARD_INVARIANT_RIBBON_COLLISION_FOUND" if terminal_standard else
            "CLS_RED_RIBBON_COLLISION_FOUND_REQUIRING_FIELD_IDENTIFICATION" if tentative_cls else
            "NO_COLLISION_IN_COMPLETED_PUBLIC_TABLE_JOIN"
        ),
    })

    (OUT / "collision_hits.json").write_text(json.dumps(hits, indent=2))
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))

    with (OUT / "collision_hits.csv").open("w", newline="", encoding="utf-8") as fh:
        fields = [
            "toxic_name", "ribbon_name", "toxic_side", "standard_nonslice",
            "standard_witnesses", "CLS_red_nonslice", "CLS_red_raw",
            "ribbon_known_slice_raw", "ribbon_kind", "zero_surgery_iso_sig",
            "framings", "linking_matrix",
        ]
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for hit in hits:
            row = {k: hit.get(k, "") for k in fields}
            row["standard_witnesses"] = "; ".join(hit["standard_witnesses"])
            writer.writerow(row)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
