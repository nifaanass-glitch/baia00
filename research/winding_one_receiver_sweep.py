#!/usr/bin/env python3
"""Proof-carrying screen for the winding-one balanced receiver family."""

from __future__ import annotations

import json
import math
import os
import time
import traceback
from pathlib import Path

import knot_floer_homology as hfk
import snappy
from snappy.exterior_to_link import exterior_to_link
from snappy.exterior_to_link.main import filled_is_3sphere
from spherogram import ClosedBraid, Link


def inverse_word(word: list[int]) -> list[int]:
    return [-x for x in reversed(word)]


def simplify_link(link: Link, rounds: int = 6) -> Link:
    for _ in range(rounds):
        for mode in ("basic", "pickup", "global"):
            try:
                link.simplify(mode)
            except BaseException:
                pass
    return link


def component_sum(link: Link, component: int, knot: Link) -> Link:
    first = link.copy()
    second = knot.copy()
    endpoint = first.link_components[component][0]
    f1, i1 = endpoint.crossing, endpoint.strand_index
    f2, i2 = f1.adjacent[i1]
    g1, j1 = second.crossings[0], 0
    g2, j2 = g1.adjacent[j1]
    f1[i1] = g2[j2]
    f2[i2] = g1[j1]
    first.crossings.extend(second.crossings)
    return Link(first.crossings)


def companion(name: str) -> Link | None:
    if name == "unknot":
        return None
    if name.startswith("T2_"):
        n = int(name.split("_")[1])
        return ClosedBraid(*([1] * n))
    if name.startswith("mirror_T2_"):
        n = int(name.split("_")[-1])
        return ClosedBraid(*([-1] * n))
    if name == "T3_4":
        return ClosedBraid(*([1, 2] * 4))
    if name == "trefoil_sum2":
        trefoil = Link(ClosedBraid(*([1] * 3)).PD_code())
        return trefoil.connected_sum(trefoil)
    raise ValueError(f"unknown companion {name!r}")


def find_s3_slopes(exterior: snappy.Manifold) -> list[dict]:
    slopes: list[tuple[int, int]] = []
    try:
        for slope in exterior.short_slopes(length=8.0)[0]:
            value = tuple(map(int, slope))
            if value not in slopes:
                slopes.append(value)
    except BaseException:
        pass
    for radius in (4, 8, 12):
        for p in range(-radius, radius + 1):
            for q in range(-radius, radius + 1):
                if math.gcd(p, q) == 1 and (p, q) not in slopes:
                    slopes.append((p, q))

    answer: list[dict] = []
    for slope in slopes:
        filling = exterior.copy()
        filling.dehn_fill(slope)
        try:
            if str(filling.homology()) != "0":
                continue
        except BaseException:
            continue
        methods: list[str] = []
        try:
            if filled_is_3sphere(filling):
                methods.append("filled_is_3sphere")
        except BaseException:
            pass
        try:
            group = filling.fundamental_group()
            group.simplify()
            if group.num_generators() == 0:
                methods.append("trivial_simplified_group")
        except BaseException:
            pass
        if methods:
            answer.append({"slope": slope, "methods": methods})
    return answer


def extract_diagrams(exterior: snappy.Manifold, s3_records: list[dict]) -> tuple[list[dict], dict | None]:
    diagrams: list[dict] = []
    best: dict | None = None
    for record in s3_records[:3]:
        slope = tuple(record["slope"])
        filling = exterior.copy()
        filling.dehn_fill(slope)
        for seed in range(6):
            started = time.time()
            try:
                knot = exterior_to_link(
                    filling,
                    verbose=False,
                    check_input=True,
                    check_answer=True,
                    careful_perturbation=True,
                    simplify_link=True,
                    pachner_search_tries=40,
                    seed=seed,
                )
                simplify_link(knot, 10)
                pd = [list(x) for x in knot.PD_code()]
                morse = hfk.pd_to_morse(pd)
                item = {
                    "slope": slope,
                    "seed": seed,
                    "crossings": len(knot.crossings),
                    "pd": pd,
                    "dt": str(knot.DT_code()),
                    "morse_girth": morse.get("girth"),
                    "morse_events": len(morse.get("events", [])),
                    "elapsed_seconds": time.time() - started,
                }
                diagrams.append(item)
                print("DIAGRAM", slope, seed, item["crossings"], item["morse_girth"], flush=True)
                if best is None or (item["morse_girth"], item["crossings"]) < (
                    best["morse_girth"],
                    best["crossings"],
                ):
                    best = item
            except BaseException as exc:
                diagrams.append(
                    {
                        "slope": slope,
                        "seed": seed,
                        "error": repr(exc),
                        "traceback": traceback.format_exc(),
                        "elapsed_seconds": time.time() - started,
                    }
                )
                print("DIAGRAM_ERROR", slope, seed, repr(exc), flush=True)
    return diagrams, best


def main() -> None:
    name = os.environ["COMPANION"]
    output_dir = Path(os.environ.get("OUTPUT_DIR", f"out/winding-one-v2-{name}"))
    output_dir.mkdir(parents=True, exist_ok=True)
    out: dict = {
        "companion": name,
        "construction": "beta_11=A*B*[A,B]^2",
        "receiver_chain_vector_expected": [1, 1],
    }
    started = time.time()

    try:
        A = [2, 1, 1, -2]
        B = [2, 2]
        commutator = A + B + inverse_word(A) + inverse_word(B)
        braid_word = A + B + commutator + commutator
        base = Link(ClosedBraid(*braid_word).PD_code())
        knot = companion(name)
        link = base if knot is None else component_sum(base, 0, knot)
        simplify_link(link, 4)

        linking_matrix = link.linking_matrix()
        out["braid_word"] = braid_word
        out["handle_link_crossings"] = len(link.crossings)
        out["linking_matrix"] = linking_matrix
        print("HANDLE", name, len(link.crossings), linking_matrix, flush=True)
        if len(link.link_components) != 3:
            raise RuntimeError(f"expected three components, got {len(link.link_components)}")

        graph_degrees = [sum(1 for value in row if value != 0) for row in linking_matrix]
        handles = [i for i, degree in enumerate(graph_degrees) if degree == 2]
        if len(handles) != 1:
            raise RuntimeError(f"expected one component linking both dotted circles, degrees={graph_degrees}")
        handle = handles[0]
        dotted = [i for i in range(3) if i != handle]
        attaching_vector = [linking_matrix[handle][d] for d in dotted]
        out["two_handle_component"] = handle
        out["dotted_components"] = dotted
        out["attaching_vector"] = attaching_vector
        if sorted(abs(x) for x in attaching_vector) != [1, 1]:
            raise RuntimeError(f"expected attaching vector (±1,±1), got {attaching_vector}")

        dotted_link = link.sublink(dotted)
        dotted_initial_components = len(dotted_link.link_components) + getattr(
            dotted_link, "unlinked_unknot_components", 0
        )
        simplify_link(dotted_link, 8)
        out["dotted_pair_initial_components"] = dotted_initial_components
        out["dotted_pair_crossings_after_simplify"] = len(dotted_link.crossings)
        out["dotted_pair_drawn_components_after"] = len(dotted_link.link_components)
        out["dotted_pair_detached_counter_after"] = getattr(
            dotted_link, "unlinked_unknot_components", 0
        )
        out["dotted_unlink_reidemeister_certificate"] = bool(
            dotted_initial_components == 2 and len(dotted_link.crossings) == 0
        )
        print("DOTTED", dotted, out["dotted_unlink_reidemeister_certificate"], flush=True)
        if not out["dotted_unlink_reidemeister_certificate"]:
            raise RuntimeError("dotted pair did not Reidemeister-reduce to a crossing-free 2-component diagram")

        marked_duals: list[dict] = []
        for keep in dotted:
            candidate: dict = {
                "keep_cusp": keep,
                "fill_cusps": [i for i in range(3) if i != keep],
            }
            print("KEEP", keep, flush=True)
            try:
                manifold = link.exterior()
                for cusp in candidate["fill_cusps"]:
                    manifold.dehn_fill((0, 1), cusp)
                exterior = snappy.Manifold(manifold.filled_triangulation())
                candidate["candidate_tets"] = exterior.num_tetrahedra()
                candidate["candidate_homology"] = str(exterior.homology())
                candidate["candidate_cusps"] = exterior.num_cusps()
                try:
                    candidate["candidate_signature"] = str(exterior.isometry_signature(of_link=True))
                except BaseException as exc:
                    candidate["candidate_signature_error"] = repr(exc)
                print(
                    "EXTERIOR",
                    keep,
                    candidate["candidate_cusps"],
                    candidate["candidate_tets"],
                    candidate["candidate_homology"],
                    flush=True,
                )

                s3_records = find_s3_slopes(exterior)
                candidate["s3_slopes"] = s3_records
                for record in s3_records:
                    print("S3", keep, record["slope"], record["methods"], flush=True)
                if not s3_records:
                    raise RuntimeError("no certified S3 meridian found")

                diagrams, best = extract_diagrams(exterior, s3_records)
                candidate["diagrams"] = diagrams
                candidate["best_diagram"] = best
                if best is None:
                    raise RuntimeError("no marked diagram extracted")
                if best["morse_girth"] > 20:
                    raise RuntimeError(
                        f"Morse girth {best['morse_girth']} exceeds the terminal HFK threshold"
                    )

                answer = hfk.pd_to_hfk(best["pd"], prime=2, complex=False)
                candidate["hfk_summary"] = {
                    key: answer.get(key)
                    for key in (
                        "tau",
                        "nu",
                        "epsilon",
                        "seifert_genus",
                        "fibered",
                        "L_space_knot",
                        "total_rank",
                        "modulus",
                    )
                }
                candidate["terminal_nonslice_hfk"] = bool(
                    candidate["hfk_summary"].get("tau")
                    or candidate["hfk_summary"].get("nu")
                    or candidate["hfk_summary"].get("epsilon")
                )
                print("HFK", keep, json.dumps(candidate["hfk_summary"], sort_keys=True), flush=True)
            except BaseException as exc:
                candidate["fatal_error"] = repr(exc)
                candidate["traceback"] = traceback.format_exc()
                print("CANDIDATE_FATAL", keep, repr(exc), flush=True)
            marked_duals.append(candidate)

        out["marked_duals"] = marked_duals
        out["terminal_hits"] = [
            candidate for candidate in marked_duals if candidate.get("terminal_nonslice_hfk")
        ]
    except BaseException as exc:
        out["fatal_error"] = repr(exc)
        out["traceback"] = traceback.format_exc()
        print("FATAL", repr(exc), flush=True)

    out["elapsed_seconds"] = time.time() - started
    with (output_dir / "result.json").open("w") as stream:
        json.dump(out, stream, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
