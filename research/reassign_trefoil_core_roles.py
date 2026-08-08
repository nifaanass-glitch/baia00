#!/usr/bin/env python3
"""Search all dotted/2-handle reassignments after trefoil insertion.

The earlier audit kept the original handle labels fixed and rejected a
38-crossing knot with tau=nu=epsilon=1.  A Kirby diagram only depends on the
chosen dotted components and framed 2-handle components, so after local
knotting we must test every role assignment.  This program does so with
fail-closed homology, longitude, S3-filling and HFK gates.
"""
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
from spherogram import ClosedBraid, Link

import balanced_rbg_word_sweep_v2 as util
from balanced_rbg_word_sweep_v3 import inherited_longitude_certificate

A = [2, 1, 1, -2]
B = [2, 2]
COMM = A + B + util.inv(A) + util.inv(B)
BETA = A + COMM + COMM


def unlink_certificate(parent, components):
    sub = parent.sublink(components)
    before = len(sub.crossings)
    util.simplify(sub, 12)
    detached = int(getattr(sub, "unlinked_unknot_components", 0))
    certified = len(components) == 2 and len(sub.crossings) == 0 and (
        detached == 2 or before == 0
    )
    return {
        "components": components,
        "crossings_before": before,
        "crossings_after": len(sub.crossings),
        "drawn_components_after": len(sub.link_components),
        "detached_unlinks_after": detached,
        "certified": certified,
        "rule": (
            "two selected parent components; after Reidemeister simplification the "
            "induced sublink is crossing-free, with either two detached counters or "
            "an originally crossing-free two-circle diagram"
        ),
    }


def search_marked_dual(link, red, dotted, keep):
    record = {
        "red_component": red,
        "dotted_components": dotted,
        "keep_component": keep,
        "fill_components": [i for i in range(3) if i != keep],
        "inherited_longitude": (0, 1),
    }
    manifold = link.exterior()
    for cusp in record["fill_components"]:
        manifold.dehn_fill((0, 1), cusp)
    exterior = snappy.Manifold(manifold.filled_triangulation())
    record["exterior_homology"] = str(exterior.homology())
    record["exterior_tetrahedra"] = exterior.num_tetrahedra()
    record["longitude_certificate"] = inherited_longitude_certificate(exterior, (0, 1))
    try:
        record["decorated_signature"] = str(exterior.isometry_signature(of_link=True))
    except BaseException as exc:
        record["signature_error"] = repr(exc)

    s3 = util.find_s3(exterior, (0, 1))
    record["s3_slopes"] = s3
    if not s3:
        raise RuntimeError("no distance-one certified S3 meridian")

    best = None
    diagrams = []
    for srec in s3[:2]:
        filling = exterior.copy()
        filling.dehn_fill(tuple(srec["slope"]))
        for seed in range(6):
            t0 = time.time()
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
                util.simplify(knot, 12)
                pd = [list(x) for x in knot.PD_code()]
                morse = hfk.pd_to_morse(pd)
                item = {
                    "slope": tuple(srec["slope"]),
                    "seed": seed,
                    "crossings": len(knot.crossings),
                    "pd": pd,
                    "dt": str(knot.DT_code()),
                    "morse_girth": morse.get("girth"),
                    "morse_events": len(morse.get("events", [])),
                    "elapsed_seconds": time.time() - t0,
                }
                try:
                    isos = exterior.is_isometric_to(knot.exterior(), return_isometries=True)
                    item["exterior_isometries"] = [str(iso) for iso in isos]
                except BaseException as exc:
                    item["isometry_error"] = repr(exc)
                diagrams.append(item)
                key = (
                    item["morse_girth"] if item["morse_girth"] is not None else 10**9,
                    item["crossings"],
                )
                old = None if best is None else (
                    best["morse_girth"] if best["morse_girth"] is not None else 10**9,
                    best["crossings"],
                )
                if old is None or key < old:
                    best = item
            except BaseException as exc:
                diagrams.append(
                    {"seed": seed, "error": repr(exc), "traceback": traceback.format_exc()}
                )
    record["diagrams"] = diagrams
    record["best_diagram"] = best
    if best is None:
        raise RuntimeError("no marked knot diagram extracted")
    if best["morse_girth"] is None or best["morse_girth"] > 20:
        raise RuntimeError(f"Morse girth {best['morse_girth']} exceeds terminal threshold")

    answer = hfk.pd_to_hfk(best["pd"], prime=2, complex=False)
    record["hfk_summary"] = {
        key: answer.get(key)
        for key in (
            "tau", "nu", "epsilon", "seifert_genus", "fibered",
            "L_space_knot", "total_rank", "modulus",
        )
    }
    record["terminal_nonslice_hfk"] = bool(
        record["hfk_summary"].get("tau")
        or record["hfk_summary"].get("nu")
        or record["hfk_summary"].get("epsilon")
    )
    return record


def main():
    infected_component = int(os.environ.get("INFECT_COMPONENT", "1"))
    output_dir = Path(
        os.environ.get("OUTPUT_DIR", f"out/reassign-trefoil-ci{infected_component}")
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    out = {
        "schema": "trefoil-infected-receiver-role-reassignment-v1",
        "infected_original_component": infected_component,
        "roles": [],
    }

    try:
        original = Link(ClosedBraid(*BETA).PD_code())
        link = util.simplify(
            util.component_sum(original, infected_component, ClosedBraid(1, 1, 1)), 6
        )
        if len(link.link_components) != 3:
            raise RuntimeError(f"infected link has {len(link.link_components)} components")
        lm = link.linking_matrix()
        out.update(
            handle_PD=[list(x) for x in link.PD_code()],
            handle_crossings=len(link.crossings),
            linking_matrix=lm,
        )

        for red in range(3):
            dotted = [i for i in range(3) if i != red]
            role = {
                "red_component": red,
                "dotted_components": dotted,
                "attaching_vector": [int(lm[red][i]) for i in dotted],
                "marked_duals": [],
            }
            vector = role["attaching_vector"]
            role["primitive_attaching_vector"] = math.gcd(*[abs(x) for x in vector]) == 1
            role["dotted_unlink"] = unlink_certificate(link, dotted)
            role["receiver_structurally_valid"] = bool(
                role["primitive_attaching_vector"] and role["dotted_unlink"]["certified"]
            )
            if role["receiver_structurally_valid"]:
                for keep in dotted:
                    try:
                        record = search_marked_dual(link, red, dotted, keep)
                    except BaseException as exc:
                        record = {
                            "red_component": red,
                            "dotted_components": dotted,
                            "keep_component": keep,
                            "fatal_error": repr(exc),
                            "traceback": traceback.format_exc(),
                        }
                    role["marked_duals"].append(record)
            out["roles"].append(role)

        out["terminal_hits"] = [
            {"role": role, "marked_dual": cand}
            for role in out["roles"]
            if role.get("receiver_structurally_valid")
            for cand in role.get("marked_duals", [])
            if cand.get("terminal_nonslice_hfk")
            and cand.get("longitude_certificate")
            and cand.get("s3_slopes")
            and "fatal_error" not in cand
        ]
    except BaseException as exc:
        out["fatal_error"] = repr(exc)
        out["traceback"] = traceback.format_exc()

    out["elapsed_seconds"] = time.time() - started
    with (output_dir / "result.json").open("w") as stream:
        json.dump(out, stream, indent=2, sort_keys=True)
    print(json.dumps(
        {
            "infected_component": infected_component,
            "structurally_valid_roles": sum(
                bool(r.get("receiver_structurally_valid")) for r in out.get("roles", [])
            ),
            "terminal_hits": len(out.get("terminal_hits", [])),
            "fatal_error": out.get("fatal_error"),
        },
        sort_keys=True,
    ))


if __name__ == "__main__":
    main()
