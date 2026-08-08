#!/usr/bin/env python3
"""Corrected balanced RBG receiver sweep without Sage dependencies.

The retained cusp is in the original link basis.  Since its exterior has
H_1 = Z, a primitive slope is the homological longitude exactly when filling
that slope still has H_1 = Z.  We certify the inherited zero-framing (0,1)
by this criterion and then search only S3 fillings at distance one from it.
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

import balanced_rbg_word_sweep_v2 as base


def inherited_longitude_certificate(exterior, slope=(0, 1)):
    if math.gcd(abs(int(slope[0])), abs(int(slope[1]))) != 1:
        raise RuntimeError(f"inherited slope is not primitive: {slope}")
    if str(exterior.homology()) != "Z":
        raise RuntimeError(f"retained exterior homology is {exterior.homology()}, not Z")
    filled = exterior.copy()
    filled.dehn_fill(slope)
    filled_h1 = str(filled.homology())
    if filled_h1 != "Z":
        raise RuntimeError(
            f"inherited zero-surgery slope {slope} gives H1={filled_h1}, not Z"
        )
    # For H1(M)=Z, quotienting by the image of a primitive boundary slope
    # leaves Z iff that image is zero.  Hence the slope lies in the kernel of
    # H1(boundary)->H1(M), and being primitive it is the homological longitude.
    return {
        "slope": tuple(map(int, slope)),
        "exterior_homology": "Z",
        "longitude_filling_homology": "Z",
        "proof": (
            "H1(exterior)=Z and filling the primitive inherited slope leaves H1=Z; "
            "therefore its image in H1(exterior) is zero and it is the primitive "
            "homological longitude"
        ),
    }


def main():
    m = int(os.environ["COMMUTATOR_EXPONENT"])
    eps = int(os.environ["B_SIGN"])
    if eps not in (-1, 1) or m == 0:
        raise ValueError("B_SIGN must be +/-1 and COMMUTATOR_EXPONENT nonzero")
    label = f"m{m:+d}_b{eps:+d}".replace("+", "p").replace("-", "n")
    output_dir = Path(os.environ.get("OUTPUT_DIR", f"out/rbg-v3-{label}"))
    output_dir.mkdir(parents=True, exist_ok=True)
    out = {
        "schema": "balanced-rbg-word-sweep-v3",
        "label": label,
        "m": m,
        "B_sign": eps,
        "companion": "positive_trefoil",
    }
    started = time.time()

    try:
        A = [2, 1, 1, -2]
        B = [2, 2]
        comm = A + B + base.inv(A) + base.inv(B)
        word = A + base.power(B, eps) + base.power(comm, m)
        raw = Link(ClosedBraid(*word).PD_code())
        lm0 = raw.linking_matrix()
        degrees0 = [sum(value != 0 for value in row) for row in lm0]
        handles0 = [i for i, degree in enumerate(degrees0) if degree == 2]
        if len(handles0) != 1:
            raise RuntimeError(f"base word has no unique red component: {lm0}")

        link = base.simplify(
            base.component_sum(raw, handles0[0], ClosedBraid(1, 1, 1)), 5
        )
        lm = link.linking_matrix()
        degrees = [sum(value != 0 for value in row) for row in lm]
        handles = [i for i, degree in enumerate(degrees) if degree == 2]
        if len(handles) != 1:
            raise RuntimeError(f"infected link has no unique red component: {lm}")
        red = handles[0]
        dotted = [i for i in range(3) if i != red]
        vector = [int(lm[red][d]) for d in dotted]
        if math.gcd(*[abs(x) for x in vector]) != 1:
            raise RuntimeError(f"nonprimitive attaching vector: {vector}")
        out.update(
            braid_word=word,
            handle_link_crossings=len(link.crossings),
            linking_matrix=lm,
            red_component=red,
            dotted_components=dotted,
            attaching_vector=vector,
            handle_PD=[list(x) for x in link.PD_code()],
        )

        dotted_link = link.sublink(dotted)
        crossings_before = len(dotted_link.crossings)
        base.simplify(dotted_link, 9)
        out["dotted_unlink"] = {
            "selected_parent_components": len(dotted),
            "crossings_before": crossings_before,
            "crossings_after": len(dotted_link.crossings),
            "certified": len(dotted) == 2 and len(dotted_link.crossings) == 0,
            "certificate_rule": (
                "the induced sublink consists of two distinct parent components and "
                "simplifies to a crossing-free diagram"
            ),
        }
        if not out["dotted_unlink"]["certified"]:
            raise RuntimeError("dotted pair is not certified as the two-component unlink")

        marked = []
        inherited_longitude = (0, 1)
        for keep in dotted:
            record = {
                "keep_cusp": keep,
                "fill_cusps": [i for i in range(3) if i != keep],
                "inherited_zero_surgery_slope": inherited_longitude,
            }
            try:
                manifold = link.exterior()
                for cusp in record["fill_cusps"]:
                    manifold.dehn_fill((0, 1), cusp)
                exterior = snappy.Manifold(manifold.filled_triangulation())
                record["tets"] = exterior.num_tetrahedra()
                record["homology"] = str(exterior.homology())
                record["longitude_certificate"] = inherited_longitude_certificate(
                    exterior, inherited_longitude
                )
                try:
                    record["signature"] = str(exterior.isometry_signature(of_link=True))
                except BaseException as exc:
                    record["signature_error"] = repr(exc)

                s3 = base.find_s3(exterior, inherited_longitude)
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
                            base.simplify(knot, 10)
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
                                isos = exterior.is_isometric_to(
                                    knot.exterior(), return_isometries=True
                                )
                                item["exterior_isometries"] = [str(iso) for iso in isos]
                            except BaseException as exc:
                                item["isometry_error"] = repr(exc)
                            diagrams.append(item)
                            key = (
                                item["morse_girth"]
                                if item["morse_girth"] is not None
                                else 10**9,
                                item["crossings"],
                            )
                            old = None if best is None else (
                                best["morse_girth"]
                                if best["morse_girth"] is not None
                                else 10**9,
                                best["crossings"],
                            )
                            if old is None or key < old:
                                best = item
                            print(
                                "DIAGRAM", label, keep, seed,
                                item["crossings"], item["morse_girth"], flush=True,
                            )
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
                    raise RuntimeError(f"Morse girth {best['morse_girth']} exceeds threshold")
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
            and r.get("longitude_certificate")
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
