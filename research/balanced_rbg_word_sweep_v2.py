#!/usr/bin/env python3
"""Fail-closed screen for balanced RBG receiver words A B^eps [A,B]^m.

Version 2 repairs the crossing-free dotted-sublink bug in the first sweep and
requires a certified S3 filling at distance one from the actual homological
longitude of the retained cusp.
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
from snappy.exterior_to_link.main import filled_is_3sphere
from spherogram import ClosedBraid, Link


def inv(word):
    return [-x for x in reversed(word)]


def power(word, exponent):
    return word * exponent if exponent >= 0 else inv(word) * (-exponent)


def simplify(link, rounds=7):
    for _ in range(rounds):
        for mode in ("basic", "pickup", "global"):
            try:
                link.simplify(mode)
            except BaseException:
                pass
    return link


def component_sum(link, component, knot):
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


def det(a, b):
    return int(a[0]) * int(b[1]) - int(a[1]) * int(b[0])


def find_s3(exterior, longitude):
    slopes = [(1, 0), (-1, 0)]
    try:
        for slope in exterior.short_slopes(length=7.0)[0]:
            value = tuple(map(int, slope))
            if value not in slopes:
                slopes.append(value)
    except BaseException:
        pass
    for radius in (4, 8):
        for p in range(-radius, radius + 1):
            for q in range(-radius, radius + 1):
                if math.gcd(p, q) == 1 and (p, q) not in slopes:
                    slopes.append((p, q))
    answer = []
    for slope in slopes:
        # A knot meridian must meet the homological longitude exactly once.
        if abs(det(slope, longitude)) != 1:
            continue
        filled = exterior.copy()
        filled.dehn_fill(slope)
        try:
            if str(filled.homology()) != "0":
                continue
        except BaseException:
            continue
        methods = []
        try:
            if filled_is_3sphere(filled):
                methods.append("filled_is_3sphere")
        except BaseException:
            pass
        try:
            group = filled.fundamental_group()
            group.simplify()
            if group.num_generators() == 0:
                methods.append("trivial_group")
        except BaseException:
            pass
        if methods:
            answer.append(
                {
                    "slope": slope,
                    "methods": methods,
                    "intersection_with_longitude": abs(det(slope, longitude)),
                }
            )
    return answer


def main():
    m = int(os.environ["COMMUTATOR_EXPONENT"])
    eps = int(os.environ["B_SIGN"])
    if eps not in (-1, 1) or m == 0:
        raise ValueError("B_SIGN must be +/-1 and COMMUTATOR_EXPONENT nonzero")
    label = f"m{m:+d}_b{eps:+d}".replace("+", "p").replace("-", "n")
    output_dir = Path(os.environ.get("OUTPUT_DIR", f"out/rbg-v2-{label}"))
    output_dir.mkdir(parents=True, exist_ok=True)
    out = {
        "schema": "balanced-rbg-word-sweep-v2",
        "label": label,
        "m": m,
        "B_sign": eps,
        "companion": "positive_trefoil",
    }
    started = time.time()

    try:
        A = [2, 1, 1, -2]
        B = [2, 2]
        comm = A + B + inv(A) + inv(B)
        word = A + power(B, eps) + power(comm, m)
        base = Link(ClosedBraid(*word).PD_code())
        lm0 = base.linking_matrix()
        degrees0 = [sum(value != 0 for value in row) for row in lm0]
        handles0 = [i for i, degree in enumerate(degrees0) if degree == 2]
        if len(handles0) != 1:
            raise RuntimeError(f"base word has no unique red component: {lm0}")
        red0 = handles0[0]

        trefoil = ClosedBraid(1, 1, 1)
        link = simplify(component_sum(base, red0, trefoil), 5)
        lm = link.linking_matrix()
        degrees = [sum(value != 0 for value in row) for row in lm]
        handles = [i for i, degree in enumerate(degrees) if degree == 2]
        if len(handles) != 1:
            raise RuntimeError(f"infected link has no unique red component: {lm}")
        red = handles[0]
        dotted = [i for i in range(3) if i != red]
        attaching_vector = [int(lm[red][d]) for d in dotted]
        if math.gcd(*[abs(x) for x in attaching_vector]) != 1:
            raise RuntimeError(f"nonprimitive attaching vector: {attaching_vector}")
        out.update(
            {
                "braid_word": word,
                "handle_link_crossings": len(link.crossings),
                "linking_matrix": lm,
                "red_component": red,
                "dotted_components": dotted,
                "attaching_vector": attaching_vector,
            }
        )
        print("LINK", label, len(link.crossings), lm, flush=True)

        # Spherogram represents a crossing-free sublink by an empty crossing list
        # and may discard its component counter.  The selected components are known
        # from the parent link.  Thus two selected parent components plus zero
        # surviving sublink crossings is an exact two-component unlink certificate.
        dotted_link = link.sublink(dotted)
        crossings_before = len(dotted_link.crossings)
        simplify(dotted_link, 9)
        out["dotted_unlink"] = {
            "selected_parent_components": len(dotted),
            "crossings_before": crossings_before,
            "crossings_after": len(dotted_link.crossings),
            "drawn_components_after": len(dotted_link.link_components),
            "detached_counter_after": getattr(dotted_link, "unlinked_unknot_components", 0),
            "certified": len(dotted) == 2 and len(dotted_link.crossings) == 0,
            "certificate_rule": "two distinct parent components and crossing-free induced sublink",
        }
        if not out["dotted_unlink"]["certified"]:
            raise RuntimeError("blue-green sublink is not certified as the two-component unlink")

        marked = []
        for keep in dotted:
            record = {
                "keep_cusp": keep,
                "fill_cusps": [i for i in range(3) if i != keep],
                "inherited_zero_surgery_slope": (0, 1),
            }
            try:
                manifold = link.exterior()
                for cusp in record["fill_cusps"]:
                    manifold.dehn_fill((0, 1), cusp)
                exterior = snappy.Manifold(manifold.filled_triangulation())
                record["tets"] = exterior.num_tetrahedra()
                record["homology"] = str(exterior.homology())
                if record["homology"] != "Z":
                    raise RuntimeError(f"retained exterior homology is {record['homology']}, not Z")
                longitude = exterior.homological_longitude()
                if longitude is None:
                    raise RuntimeError("retained cusp has no homological longitude")
                longitude = tuple(map(int, longitude))
                record["homological_longitude"] = longitude
                record["inherited_slope_is_homological_longitude"] = longitude in ((0, 1), (0, -1))
                if not record["inherited_slope_is_homological_longitude"]:
                    raise RuntimeError(
                        f"inherited zero-surgery slope (0,1) is not the homological longitude {longitude}"
                    )
                try:
                    record["signature"] = str(exterior.isometry_signature(of_link=True))
                except BaseException as exc:
                    record["signature_error"] = repr(exc)

                s3 = find_s3(exterior, longitude)
                record["s3_slopes"] = s3
                print("S3", label, keep, s3, flush=True)
                if not s3:
                    raise RuntimeError("no certified distance-one S3 meridian")

                best = None
                diagrams = []
                for srec in s3[:2]:
                    filling = exterior.copy()
                    filling.dehn_fill(tuple(srec["slope"]))
                    for seed in range(5):
                        t0 = time.time()
                        try:
                            knot = exterior_to_link(
                                filling,
                                verbose=False,
                                check_input=True,
                                check_answer=True,
                                careful_perturbation=True,
                                simplify_link=True,
                                pachner_search_tries=35,
                                seed=seed,
                            )
                            simplify(knot, 10)
                            pd = [list(x) for x in knot.PD_code()]
                            morse = hfk.pd_to_morse(pd)
                            item = {
                                "slope": tuple(srec["slope"]),
                                "seed": seed,
                                "crossings": len(knot.crossings),
                                "morse_girth": morse.get("girth"),
                                "morse_events": len(morse.get("events", [])),
                                "pd": pd,
                                "dt": str(knot.DT_code()),
                                "elapsed_seconds": time.time() - t0,
                            }
                            try:
                                isos = exterior.is_isometric_to(knot.exterior(), return_isometries=True)
                                item["exterior_isometries"] = [str(iso) for iso in isos]
                            except BaseException as exc:
                                item["isometry_error"] = repr(exc)
                            diagrams.append(item)
                            print(
                                "DIAGRAM", label, keep, seed,
                                item["crossings"], item["morse_girth"], flush=True,
                            )
                            key = (
                                item["morse_girth"] if item["morse_girth"] is not None else 10**9,
                                item["crossings"],
                            )
                            if best is None:
                                best = item
                            else:
                                old = (
                                    best["morse_girth"] if best["morse_girth"] is not None else 10**9,
                                    best["crossings"],
                                )
                                if key < old:
                                    best = item
                        except BaseException as exc:
                            diagrams.append(
                                {
                                    "seed": seed,
                                    "error": repr(exc),
                                    "traceback": traceback.format_exc(),
                                }
                            )
                record["diagrams"] = diagrams
                record["best_diagram"] = best
                if best is None:
                    raise RuntimeError("no marked knot diagram")
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
                print("HFK", label, keep, json.dumps(record["hfk_summary"], sort_keys=True), flush=True)
            except BaseException as exc:
                record["fatal_error"] = repr(exc)
                record["traceback"] = traceback.format_exc()
                print("MARKED_ERROR", label, keep, repr(exc), flush=True)
            marked.append(record)
        out["marked_duals"] = marked
        out["terminal_hits"] = [
            r for r in marked
            if r.get("terminal_nonslice_hfk")
            and r.get("inherited_slope_is_homological_longitude")
            and r.get("s3_slopes")
            and "fatal_error" not in r
        ]
    except BaseException as exc:
        out["fatal_error"] = repr(exc)
        out["traceback"] = traceback.format_exc()
        print("FATAL", label, repr(exc), flush=True)

    out["elapsed_seconds"] = time.time() - started
    with (output_dir / "result.json").open("w") as stream:
        json.dump(out, stream, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
