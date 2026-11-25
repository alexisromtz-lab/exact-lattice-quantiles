"""Command line: the same numbers the library computes, printable and testable.

Exit codes are part of the contract:

``0``  the request was answered
``1``  the answer is no: the margin requested is above the instance ceiling --
       impossible at any price, and the ceiling is printed so the caller learns
       the number -- or an ``acta`` anchor does not reproduce, or a claim the
       package makes about its own lemma is refuted by its own census
``2``  the margin requested exists but costs more points than the budget
``3``  the request itself is malformed
``4``  the acta could not be measured (no git, no such table, not exact): the
       boundary is ``NO_MEDIDO``, and ``NO_MEDIDO`` never degrades to ``0``
``5``  the instance has an answer and its table is above :data:`TABLE_BUDGET`:
       the implementation refused by name before allocating, so the shell can
       tell "too big to compute here" from "malformed" and from "impossible"
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction

from .ceiling import (InfeasibleAboveInstanceCeiling, ceiling_by_product,
                      ceiling_by_scan, ceiling_census, instance_ceiling)
from .acta import (ManifestError, SourceError, inspect_source, load_manifest,
                   parse_tables, render, shipped_manifest, verify_file,
                   verify_manifest)
from .coverage import (coverage_census, design_law, drift_law, drift_radius,
                       drift_tolerance, witness_margin)
from .arithmetic import primes_of
from .lattice import TABLE_BUDGET, MaskTableTooLarge, mask_table, universe
from .pricing import (antichain_saving_ceiling, block_catalogue,
                      block_dominance_census, block_of, block_saving,
                      carrier_block_of, carrier_block_saving,
                      exactness_frontier, pair_antichain_of,
                      pair_antichain_saving, price_exact,
                      price_form_census, price_form_sweep,
                      reciprocal_prime_sum, sign_regime, sign_theorem_census,
                      strong_lemma_frontier, worst_antichain_weight)
from .solve import OK, InfeasibleWithinBudget, ladder, solve_lattice_schedule

RC_OK = 0
RC_ABOVE_CEILING = 1
RC_OVER_BUDGET = 2
RC_BAD_REQUEST = 3
RC_NO_MEDIDO = 4
RC_BEYOND_TABLE_BUDGET = 5


def _margin(value: Fraction | None) -> str:
    """A margin that is not measured is not a margin of zero."""
    return "NO_MEDIDO (uncovered divisor, see below)" if value is None else str(value)


def _fraction(text: str) -> Fraction:
    try:
        return Fraction(text)
    except (ValueError, ZeroDivisionError) as exc:
        raise argparse.ArgumentTypeError(f"not an exact rational: {text!r}") from exc


def _print_ceiling(period: int, as_json: bool) -> int:
    c = instance_ceiling(period)
    if as_json:
        print(json.dumps({
            "period": c.period,
            "margin": str(c.margin),
            "critical_divisor": c.critical_divisor,
            "dropped": list(c.dropped),
            "prime_product": str(c.prime_product),
            "radical_bound_hit": c.radical_bound_hit,
            "universe": list(c.universe),
            "by_product": str(ceiling_by_product(period)),
            "by_scan": str(ceiling_by_scan(period)),
        }, indent=2, sort_keys=True))
        return RC_OK
    print(f"period            {c.period}")
    print(f"universe U(T)     {list(c.universe)}")
    print(f"ceiling nu(T)     {c.margin}  ({float(c.margin) * 100:.4f} %)")
    print(f"critical divisor  {c.critical_divisor}")
    print(f"P(T)              {c.prime_product}")
    print(f"dropped primes    {list(c.dropped) or 'none'}")
    print(f"radical bound hit {c.radical_bound_hit}")
    print(f"by product        {ceiling_by_product(period)}")
    print(f"by scan           {ceiling_by_scan(period)}")
    return RC_OK


def _print_solve(args: argparse.Namespace) -> int:
    try:
        schedule = solve_lattice_schedule(
            args.period, args.max_cost, args.probes,
            Fraction(0) if args.min_margin is None else args.min_margin)
    except InfeasibleAboveInstanceCeiling as exc:
        if args.json:
            print(json.dumps({
                "verdict": "INFEASIBLE_ABOVE_INSTANCE_CEILING",
                "period": exc.period,
                "needed": str(exc.needed),
                "ceiling": str(exc.ceiling),
            }, indent=2, sort_keys=True))
        else:
            print(f"INFEASIBLE_ABOVE_INSTANCE_CEILING: {exc}")
        return RC_ABOVE_CEILING
    except InfeasibleWithinBudget as exc:
        if args.json:
            print(json.dumps({
                "verdict": "OVER_BUDGET",
                "period": exc.period,
                "needed": str(exc.needed),
                "cheapest": exc.cheapest,
                "max_cost": exc.max_cost,
            }, indent=2, sort_keys=True))
        else:
            print(f"OVER_BUDGET: {exc}")
        return RC_OVER_BUDGET
    except MaskTableTooLarge as exc:
        if args.json:
            print(json.dumps({
                "verdict": "BEYOND_TABLE_BUDGET",
                "period": exc.period,
                "site": exc.site,
                "masks": f"2^{exc.count}",
                "width": exc.width,
                "cells": exc.cells,
                "budget": TABLE_BUDGET,
                "reason": str(exc),
            }, indent=2, sort_keys=True))
        else:
            print(f"BEYOND_TABLE_BUDGET: {exc}")
        return RC_BEYOND_TABLE_BUDGET
    if args.json:
        print(json.dumps(schedule.as_dict(), indent=2, sort_keys=True))
        return RC_OK if schedule.verdict == OK else RC_NO_MEDIDO
    print(f"period           {schedule.period}")
    print(f"design           {list(schedule.design)}")
    print(f"cost             {schedule.cost} points (budget {schedule.max_cost})")
    print(f"design margin    {schedule.design_margin}")
    print(f"sections ({len(schedule.probes)})  {[str(q) for q in schedule.probes]}")
    print(f"witness margin   {_margin(schedule.witness_margin)}")
    print(f"robust margin    {_margin(schedule.robust_margin)}")
    if schedule.uncovered:
        print(f"uncovered        {list(schedule.uncovered)} of {len(schedule.ceiling.universe)}"
              f" divisors -- no section lands in them, so no margin is measured")
    print(f"probes needed    {schedule.probes_needed}")
    print(f"ceiling          {schedule.ceiling.margin} "
          f"(critical divisor {schedule.ceiling.critical_divisor})")
    print(f"verdict          {schedule.verdict}")
    # A verdict that is not OK is a schedule not delivered: exiting 0 there
    # would let a shell caller read "we have a schedule" from a failed one.
    return RC_OK if schedule.verdict == OK else RC_NO_MEDIDO


def _print_drift(args: argparse.Namespace) -> int:
    """The trade-off the witness margin hides: slack is bought with tolerance."""
    try:
        schedule = solve_lattice_schedule(
            args.period, args.max_cost, args.probes,
            Fraction(0) if args.min_margin is None else args.min_margin)
    except InfeasibleAboveInstanceCeiling as exc:
        print(f"INFEASIBLE_ABOVE_INSTANCE_CEILING: {exc}", file=sys.stderr)
        return RC_ABOVE_CEILING
    except InfeasibleWithinBudget as exc:
        print(f"OVER_BUDGET: {exc}", file=sys.stderr)
        return RC_OVER_BUDGET
    except MaskTableTooLarge as exc:
        print(f"BEYOND_TABLE_BUDGET: {exc}", file=sys.stderr)
        return RC_BEYOND_TABLE_BUDGET
    U = universe(schedule.period)
    Q = schedule.probes
    if schedule.witness_margin is None:
        print(f"UNDERCOVERED: no section lands in {list(schedule.uncovered)} of "
              f"{U}, so there is no margin to trade against drift",
              file=sys.stderr)
        return RC_NO_MEDIDO
    curve = []
    for share in (Fraction(1), Fraction(3, 4), Fraction(1, 2), Fraction(1, 4), Fraction(0)):
        kept = schedule.witness_margin * share
        tol = min(drift_tolerance(schedule.design, Q, U, kept).values()) if Q else Fraction(0)
        curve.append((kept, tol))
    design = design_law(schedule.design, Q, U) if Q else None
    divisor_laws = drift_law(schedule.design, Q, U) if Q else {}
    if args.json:
        print(json.dumps({
            "period": schedule.period,
            "design": list(schedule.design),
            "cost": schedule.cost,
            "probes": [str(q) for q in Q],
            "witness_margin": None if schedule.witness_margin is None
            else str(schedule.witness_margin),
            "coverage_drift": str(design.radius if design else Fraction(0)),
            "linear_regime_end": str(design.radius if design else Fraction(0)),
            "first_divisor_break": str(design.first_divisor_break) if design else None,
            "divisor_regimes": {
                str(d): {"radius": str(law.radius),
                         "critical_margin": str(law.critical_margin),
                         "cliff": str(law.cliff)}
                for d, law in sorted(divisor_laws.items())
            },
            "kept_margin_to_drift": [[str(v), str(t)] for v, t in curve],
        }, indent=2, sort_keys=True))
        return RC_OK
    print(f"period           {schedule.period}")
    print(f"design           {list(schedule.design)}")
    print(f"cost             {schedule.cost} points (budget {schedule.max_cost})")
    print(f"sections ({len(Q)})  {[str(q) for q in Q]}")
    print(f"witness margin   {schedule.witness_margin}")
    if design is not None:
        print(f"coverage drift   {design.radius}  (captured, margin not kept)")
        print(f"tolerance(v)     {design.radius} - v for every v up to {design.radius}"
              f"  (one-for-one, no cliff)")
        if design.first_divisor_break < design.radius:
            print(f"narrowest        divisor regime ends at {design.first_divisor_break}"
                  f"  (pessimistic: the design pays nothing there)")
    print("margin kept      drift tolerated")
    for kept, tol in curve:
        print(f"  {str(kept):<15} {tol}")
    return RC_OK


def _print_ladder(period: int, as_json: bool) -> int:
    curve = ladder(period)
    if as_json:
        print(json.dumps({str(cost): str(margin) for cost, margin in sorted(curve.items())},
                         indent=2, sort_keys=True))
        return RC_OK
    for cost, margin in sorted(curve.items()):
        print(f"  {cost:5d} points -> margin {margin}")
    return RC_OK


def _print_census(low: int, high: int, as_json: bool) -> int:
    data = ceiling_census(low, high)
    if as_json:
        payload = dict(data)
        lowest = data.get("lowest")
        if isinstance(lowest, tuple):
            payload["lowest"] = {"margin": str(lowest[0]), "period": lowest[1]}
        print(json.dumps(payload, indent=2, sort_keys=True))
        return RC_OK
    print(f"periods {low}..{high}: {data['periods']}")
    print(f"  above 1/3  {data['above']}")
    print(f"  equal 1/3  {data['equal']}")
    print(f"  below 1/3  {data['below']}")
    lowest = data.get("lowest")
    if isinstance(lowest, tuple):
        print(f"  lowest     {lowest[0]} at period {lowest[1]}")
    return RC_OK


def _print_price(period: int, as_json: bool) -> int:
    rows = [price_form_census(period), price_form_census(period, reduced=False)]
    if as_json:
        print(json.dumps([row.as_dict() for row in rows], indent=2, sort_keys=True))
        return RC_OK
    for row in rows:
        print(f"{period} over {'U(T)' if row.reduced else 'every divisor'}"
              f" ({row.universe_size} elements, {row.consumers} requests)")
        print(f"  one grid per maximal    closed form, {row.frontiers} frontiers")
        print(f"  exact minimum           cheapest cover {row.cheapest_cover}"
              f" points, one grid alone {row.one_grid}")
        print(f"  requests overcharged    {row.discrepancies} ({row.share})"
              f" on {row.saving_frontiers} frontiers")
    return RC_OK


def _print_price_sweep(low: int, high: int, as_json: bool) -> int:
    data = price_form_sweep(low, high)
    if as_json:
        print(json.dumps(data, indent=2, sort_keys=True))
        return RC_OK
    print(f"even periods {low}..{high} over {data['universe']}: {data['periods']} censused")
    print(f"  closed form is the price   {data['adequate']}")
    print(f"  closed form overcharges    {data['inadequate']}")
    print(f"  cheaper than one grid      {data['cheaper_than_one_grid']}")
    cost, period = data["costliest"]
    print(f"  costliest cover            {cost} points at period {period}")
    skipped = data["skipped"]
    print(f"  skipped (unaffordable)     {len(skipped)}")
    return RC_OK


def _print_frontier(max_period: int, reduced: bool, as_json: bool) -> int:
    data = exactness_frontier(max_period, reduced=reduced)
    if as_json:
        print(json.dumps({**data, "witness": list(data["witness"])},
                         indent=2, sort_keys=True))
        return RC_OK
    if data["period"] is None:
        print(f"even periods 2..{max_period} over {data['universe']}: the closed"
              f" form is the price on every one of the {data['scanned']} censused")
        return RC_OK
    print(f"even periods 2..{max_period} over {data['universe']}: {data['breaks']}"
          f" of {data['scanned']} mispriced, first at T = {data['period']}")
    print(f"  witness                    {data['witness']}")
    print(f"  closed form                {data['bound']} points")
    print(f"  one grid {data['grid']:<18}{data['cover']} points")
    print(f"  saving                     {data['saving']} points")
    print(f"  decided by                 {data['law']}")
    return RC_OK


def _print_blocks(limit: int, as_json: bool) -> int:
    rows = block_catalogue(limit)
    if as_json:
        print(json.dumps([{"n": n, "saving": s, "block": list(b)}
                          for n, s, b in rows], indent=2, sort_keys=True))
        return RC_OK
    print(f"divisibility-minimal blocks up to {limit}: {len(rows)}")
    for n, s, b in rows:
        print(f"  {n:>6}  saves {s:>3}  {b}")
    return RC_OK


def _print_dominance(limit: int, as_json: bool) -> int:
    """The lemma behind the law, checked instead of assumed."""
    census = block_dominance_census(limit)
    if as_json:
        print(json.dumps(census, indent=2, sort_keys=True))
        return RC_OK if not census["violations"] else RC_ABOVE_CEILING
    print(f"antichains with a least common multiple up to {limit}: "
          f"{census['in_play']} periods where a saving is in play")
    if census["violations"]:
        for row in census["violations"]:
            print(f"  VIOLATION  {row['period']}: ceiling {row['ceiling']} "
                  f"on {tuple(row['witness'])} beats the block "
                  f"{row['block']} on {tuple(row['block_witness'])}")
        return RC_ABOVE_CEILING
    print("  every saving is carried by the drop-one block")
    print(f"  {census['sign_verdict']} (the reading the law needs)")
    return RC_OK


def _sci(value: int, places: int = 3) -> str:
    """Scientific notation without touching float: these savings have hundreds of digits."""
    text = str(abs(value))
    mantissa = text[0] + "." + text[1:1 + places].ljust(places, "0")
    return f"{'-' if value < 0 else ''}{mantissa}e+{len(text) - 1}"


def _print_strong_lemma(max_primes: int, as_json: bool) -> int:
    """The ranking reading of the lemma is false: publish where, in exact integers."""
    report = strong_lemma_frontier(max_primes)
    crossing = report["crossing_primes"]
    if as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return RC_OK if crossing else RC_ABOVE_CEILING
    print("ranking lemma: no antichain saves more than the block of its own "
          "least common multiple")
    if not crossing:
        print(f"  no crossing below {max_primes} primes: the claim is NOT_MEASURED,"
              " not confirmed")
        return RC_ABOVE_CEILING
    below = report["below"]
    print(f"  crossing at {crossing} primes, the largest being "
          f"{report['largest_prime']}")
    if below:
        print(f"  at {below['primes']:>3} primes (largest {below['largest_prime']}, "
              f"{below['digits']}-digit period)  the block wins by "
              f"{_sci(-below['excess'])}")
    row = report["crossing"]
    per_mille = row["excess"] * 100_000 // row["block_saving"]
    print(f"  at {crossing:>3} primes (largest {row['largest_prime']}, "
          f"{row['digits']}-digit period)  the pairs win by "
          f"{_sci(row['excess'])} "
          f"(+{per_mille // 1000}.{per_mille % 1000:03d} %)")
    print(f"  verdict  {report['verdict']}")
    return RC_OK


def _print_sign_theorem(max_period: int, as_json: bool) -> int:
    """The surviving half of the lemma: proved where the reciprocals allow it."""
    census = sign_theorem_census(max_period)
    broken = (census["lemma_violations"] or census["sign_violations"]
              or census["characterisation_mismatches"])
    if as_json:
        print(json.dumps(census, indent=2, sort_keys=True))
        return RC_ABOVE_CEILING if broken else RC_OK
    def line(label: str, value: str, note: str = "") -> str:
        return f"  {label:<32}{value:>6}  {note}"

    print("sign lemma: if no block saves, no antichain saves")
    print(line(f"periods 2..{max_period}", str(census["periods"]), "censused"))
    print(line("sum 1/p at most one", str(census["proven_zero"]),
               "proved: no antichain saves any point"))
    print(line("  the block is the heaviest", str(census["tight"]),
               "of those, the lemma is sharp"))
    print(line("the block saves already", str(census["block_saves"]),
               "nothing left to prove"))
    print(line("band, sum 1/p above one", str(census["band"]),
               "measured only, the proof does not reach"))
    for row in census["band_rows"]:
        equation = "  Giuga equation" if row["band_equation"] else ""
        print(f"    {row['period']:>6}  sum 1/p = {row['reciprocal_sum']:<10}"
              f" period*(sum-1) = {row['excess']:<2}  block {row['block']:>3}"
              f"  ceiling {row['ceiling']}{equation}")
    for row in census["lemma_violations"]:
        print(f"  VIOLATION  {row['period']}: weight {row['weight']} on "
              f"{tuple(row['witness'])} beats the bound {row['bound']}")
    for row in census["sign_violations"]:
        print(f"  VIOLATION  {row['period']}: an antichain saves "
              f"{row['ceiling']} where the block saves nothing"
              f" ({row['regime']})")
    for row in census["characterisation_mismatches"]:
        print(f"  MISMATCH  {row['period']}: band without the equation "
              f"sum 1/p = 1 + 1/period twice an odd number")
    print(f"  lemma    {census['lemma']}")
    print(f"  verdict  {census['sign_verdict']}")
    return RC_ABOVE_CEILING if broken else RC_OK


def _print_verify_acta(manifest: str, as_json: bool, source: str = "") -> int:
    """Recompute the anchors of a manifest: the boundary is a gate, not a memory."""
    try:
        data = load_manifest(manifest)
        checks = verify_manifest(data)
    except ManifestError as exc:
        print(f"BAD_REQUEST: {exc}", file=sys.stderr)
        return RC_BAD_REQUEST
    report = None
    if source:
        try:
            report = inspect_source(source, data)
        except SourceError as exc:
            print(f"NO_MEDIDO: {exc}", file=sys.stderr)
            return RC_NO_MEDIDO
        checks = checks + list(report.checks)
    divergent = [check for check in checks if not check.ok]
    anchors = len({check.anchor for check in checks})
    if as_json:
        payload = {
            "manifest": str(manifest),
            "anchors": anchors,
            "checks": len(checks),
            "divergent": len(divergent),
            "results": [{"anchor": c.anchor, "kind": c.kind, "field": c.field,
                         "manifest": render(c.expected),
                         "subject": c.subject,
                         c.subject: render(c.measured), "ok": c.ok}
                        for c in checks],
        }
        if report is not None:
            payload["source"] = {
                "path": report.path, "blob": report.blob, "lines": report.lines,
                "declared_blob": report.declared_blob,
                "declared_lines": report.declared_lines,
                "moved": report.moved, "tables": list(report.tables),
                "unanchored": list(report.unanchored),
            }
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        for check in checks:
            print(check)
        if report is not None:
            verdict = "MOVIDO" if report.moved else "IDENTICO"
            print(f"source {report.path}")
            print(f"  blob   {report.blob}  manifest {report.declared_blob}  {verdict}")
            print(f"  lines  {report.lines}  manifest {report.declared_lines}")
            print(f"  tables {', '.join(report.tables)} read by AST, never imported")
            if report.moved and not report.divergent:
                print("  the numbers hold; it is the citation of the source that is stale")
            if report.unanchored:
                shown = ", ".join(report.unanchored[:6])
                more = "" if len(report.unanchored) <= 6 else f" (+{len(report.unanchored) - 6} more)"
                print(f"  sin ancla: {len(report.unanchored)} cells no anchor declares: {shown}{more}")
        print(f"{len(checks)} checks over {anchors} anchors:"
              f" {len(checks) - len(divergent)} reproduce, {len(divergent)} divergent")
    return RC_OK if not divergent else RC_ABOVE_CEILING


def _self_test() -> int:
    """A handful of claims the package makes, checked in-process."""
    ok = True

    def check(label: str, got: object, want: object) -> None:
        nonlocal ok
        if got != want:
            ok = False
            print(f"FALLA  {label}: got {got!r}, want {want!r}")
        else:
            print(f"ok     {label}")

    c = instance_ceiling(180)
    check("nu(180) = 1/3", c.margin, Fraction(1, 3))
    check("critical divisor of 180", c.critical_divisor, 6)
    check("180 drops only the prime 5", c.dropped, (5,))
    check("1/3 = P(180) * 5/4", c.prime_product * Fraction(5, 4), c.margin)
    check("paths agree at 180", ceiling_by_product(180), ceiling_by_scan(180))
    check("nu(360) = 4/15", instance_ceiling(360).margin, Fraction(4, 15))
    check("nu(252) = 2/7", instance_ceiling(252).margin, Fraction(2, 7))
    try:
        solve_lattice_schedule(180, max_cost=10**6, k_probes=4,
                               min_margin=Fraction(2, 5))
    except InfeasibleAboveInstanceCeiling as exc:
        check("2/5 is above the ceiling of 180", exc.ceiling, Fraction(1, 3))
    else:
        check("2/5 is above the ceiling of 180", "no exception", "exception")
    s = solve_lattice_schedule(180, max_cost=64, k_probes=4, min_margin="1/3")
    check("cheapest design at 1/3", (s.cost, s.design),
          (64, (4, 6, 10, 12, 15, 18, 20, 45)))
    check("1/3 realised with 4 probes", s.witness_margin, Fraction(1, 3))
    U = universe(180)
    G, Q = s.design, s.probes
    check("no drift is exactly the coverage radius",
          min(drift_tolerance(G, Q, U, Fraction(0)).values()),
          min(drift_radius(G, Q, U).values()))
    check("the witness margin survives no drift at all",
          min(drift_tolerance(G, Q, U, s.witness_margin).values()), Fraction(0))
    d_edge, edge = min(drift_law(G, Q, U).items(),
                       key=lambda item: item[1].critical_margin)
    check("the price of a margin is one-for-one up to the critical margin",
          drift_tolerance(G, Q, U, edge.critical_margin / 2)[d_edge],
          edge.radius - edge.critical_margin / 2)
    check("and strictly worse at it",
          drift_tolerance(G, Q, U, edge.critical_margin)[d_edge],
          edge.radius - edge.critical_margin - edge.cliff)
    shape = design_law(G, Q, U)
    past = (shape.first_divisor_break + shape.radius) / 2
    check("the design itself pays one-for-one past that break",
          min(drift_tolerance(G, Q, U, past).values()), shape.radius - past)
    check("and all the way to zero tolerance",
          min(drift_tolerance(G, Q, U, shape.radius).values()), Fraction(0))
    half = drift_law(G, Q, U, doubled=False)
    d_half = min(half, key=lambda divisor: half[divisor].critical_margin)
    check("a critical margin feeds back in the undoubled view too",
          drift_tolerance(G, Q, U, half[d_half].critical_margin, doubled=False)[d_half],
          half[d_half].radius - half[d_half].critical_margin - half[d_half].cliff)
    check("block of 270 is {54, 90, 135}", block_of(270)[0], (54, 90, 135))
    check("block of 270 saves 4", block_saving(270), 4)
    check("block of 30 saves nothing", block_saving(30), 0)
    check("every divisor universe misprices first at 60",
          exactness_frontier(2000, reduced=False)["period"], 60)
    check("every carrier universe misprices first at 1080",
          exactness_frontier(2000)["period"], 1080)
    check("the universe filters the block of 420 down to three elements",
          carrier_block_of(420, universe(2520))[0], (84, 140, 210))
    check("and that filtered block still saves 7",
          carrier_block_saving(420, universe(2520)), 7)
    check("no antichain with lcm 1050 saves more than its block",
          antichain_saving_ceiling(1050)[0], block_saving(1050))
    check("and the witness of that ceiling is the block itself",
          tuple(sorted(antichain_saving_ceiling(1050)[1])), block_of(1050)[0])
    check("the block of 60 in U(840) keeps only a pair",
          carrier_block_of(60, universe(840))[0], (12, 20))
    check("and that pair costs 14 points more than the grid it replaces",
          carrier_block_saving(60, universe(840)), -14)
    check("the carrier universe keeps its first witness",
          exactness_frontier(2000)["witness"], (54, 90, 135))
    frontier = strong_lemma_frontier(200)
    check("the pairs beat the block from 176 primes on",
          (frontier["crossing_primes"], frontier["largest_prime"]), (176, 1049))
    check("and at 175 primes the block still wins, in exact integers",
          frontier["below"]["excess"] < 0, True)
    check("the pairs of 180 are the three products of two primes",
          pair_antichain_of(180), (30, 18, 12))
    check("which save nothing at all",
          pair_antichain_saving(180) < 0, True)
    check("the reciprocal sum of 30 is 31/30",
          reciprocal_prime_sum(30), Fraction(31, 30))
    check("and yet its block saves nothing, because 31 against 30 buys the "
          "same thirty points", block_saving(30), 0)
    check("858 is one plus its own reciprocal",
          reciprocal_prime_sum(858), 1 + Fraction(1, 858))
    check("so it is in the band, where no proof reaches",
          sign_regime(858)["regime"], "band")
    check("60 is not: its reciprocal sum is 31/30 and twice an odd number "
          "is not what it is", sign_regime(60)["regime"], "block-saves")
    check("36 is proved, with no band to worry about",
          sign_regime(36)["regime"], "proven-zero")
    check("the heaviest antichain of 36 is its block",
          worst_antichain_weight(36), (30, (18, 12)))
    check("the heaviest antichain of 858 weighs one more than 858 itself",
          worst_antichain_weight(858)[0], 859)
    check("which is exactly the bound the lemma gives, so it is sharp",
          worst_antichain_weight(858)[0], sum(858 // p for p in primes_of(858)))
    check("and it still saves nothing",
          antichain_saving_ceiling(858)[0], 0)
    census_300 = sign_theorem_census(300)
    check("below 300 the band has exactly one member",
          [row["period"] for row in census_300["band_rows"]], [30])
    check("and the lemma holds on every period it covers",
          census_300["lemma_violations"], [])
    check("and so does the theorem behind it",
          census_300["sign_violations"], [])
    check("the price at 1080 is 135, not 139",
          price_exact(1080, (54, 90, 135)), 135)
    try:
        anchors = verify_file(shipped_manifest())
    except ManifestError as exc:
        check("the shipped manifest reads as a manifest", f"ManifestError: {exc}",
              "readable")
    else:
        divergent = [str(c) for c in anchors if not c.ok]
        check(f"the shipped manifest reproduces its {len(anchors)} checks",
              divergent or "all reproduce", "all reproduce")
    acta = "TECHOS = {2: F(1, 2)}\n"
    check("a table is read without importing the acta",
          parse_tables(acta, ["TECHOS"]), {"TECHOS": {2: Fraction(1, 2)}})
    try:
        parse_tables("TECHOS = {2: 0.5}\n", ["TECHOS"])
    except SourceError:
        check("a float in a quoted table is refused", "SourceError", "SourceError")
    else:
        check("a float in a quoted table is refused", "accepted", "SourceError")
    try:
        parse_tables("TECHOS = dict(TECHOS)\n", ["TECHOS"])
    except SourceError:
        check("a table is never computed at read time", "SourceError", "SourceError")
    else:
        check("a table is never computed at read time", "accepted", "SourceError")
    try:
        mask_table(30, "the self-test")
    except MaskTableTooLarge as exc:
        check("a table above the budget is refused by name", exc.rc, 5)
    else:
        check("a table above the budget is refused by name", "allocated",
              "MaskTableTooLarge")
    return RC_OK if ok else RC_ABOVE_CEILING


def _print_designs(period: int, as_json: bool) -> int:
    data = coverage_census(period)
    if as_json:
        print(json.dumps(dict(data), indent=2, sort_keys=True))
        return RC_OK
    print(f"designs of period {period}: {data['designs']}")
    print(f"  covering   {data['covering']}")
    print(f"  cheapest   {data['cheapest']} points")
    return RC_OK


class _Parser(argparse.ArgumentParser):
    """An argument error is a return code, not an exit.

    ``argparse`` calls ``sys.exit(2)`` from ``error``, which would end the
    process with a code this module does not document.  Raising instead keeps
    the caller in control of ``RC_BAD_REQUEST``.
    """

    def error(self, message: str) -> None:  # type: ignore[override]
        raise ValueError(message)


def build_parser() -> argparse.ArgumentParser:
    parser = _Parser(
        prog="exact-lattice-quantiles",
        description="Exact phase margin for multi-rate periodic sampling.")
    parser.add_argument("--json", action="store_true", help="machine readable output")
    parser.add_argument("--self-test", action="store_true",
                        help="check the claims the package makes and exit")
    sub = parser.add_subparsers(dest="command")

    p_ceiling = sub.add_parser("ceiling", help="the exact ceiling of one period")
    p_ceiling.add_argument("period", type=int)

    p_solve = sub.add_parser("solve", help="cheapest design reaching a margin")
    p_solve.add_argument("period", type=int)
    p_solve.add_argument("--max-cost", type=int, required=True)
    p_solve.add_argument("--probes", type=int, required=True)
    p_solve.add_argument("--min-margin", type=_fraction, default=None)

    p_drift = sub.add_parser(
        "drift", help="margin kept against common drift, for one schedule")
    p_drift.add_argument("period", type=int)
    p_drift.add_argument("--max-cost", type=int, required=True)
    p_drift.add_argument("--probes", type=int, required=True)
    p_drift.add_argument("--min-margin", type=_fraction, default=None)

    p_ladder = sub.add_parser("ladder", help="cost -> best margin for one period")
    p_ladder.add_argument("period", type=int)

    p_census = sub.add_parser("census", help="how many periods beat, hit or miss 1/3")
    p_census.add_argument("low", type=int, nargs="?", default=2)
    p_census.add_argument("high", type=int, nargs="?", default=400)

    p_designs = sub.add_parser("designs", help="census of every design of one period")
    p_designs.add_argument("period", type=int)

    p_price = sub.add_parser(
        "price", help="closed-form price against the exact minimum, both universes")
    p_price.add_argument("period", type=int)

    p_prices = sub.add_parser("prices", help="the same census over a window of periods")
    p_prices.add_argument("low", type=int, nargs="?", default=2)
    p_prices.add_argument("high", type=int, nargs="?", default=400)

    p_frontier = sub.add_parser(
        "frontier", help="the first period the closed form misprices")
    p_frontier.add_argument("--max-period", type=int, default=8000)
    p_frontier.add_argument(
        "--universe", choices=("reduced", "all"), default="reduced",
        help="reduced carriers U(T) (searched) or every divisor D(T) (decided"
             " by the block law)")

    p_blocks = sub.add_parser(
        "blocks", help="the divisibility-minimal blocks that save anything")
    p_blocks.add_argument("--limit", type=int, default=2000)

    p_dominance = sub.add_parser(
        "dominance",
        help="check the lemma the law rests on: no antichain saves more than "
             "the block of its own least common multiple")
    p_dominance.add_argument("--max-period", type=int, default=2000)

    p_strong = sub.add_parser(
        "strong-lemma",
        help="locate where the ranking reading of the lemma fails: the block is "
             "not always the best-saving antichain")
    p_strong.add_argument("--max-primes", type=int, default=400)

    p_sign = sub.add_parser(
        "sign-theorem",
        help="prove the sign reading of the lemma where the reciprocal prime "
             "sum allows it, and name the band where it does not")
    p_sign.add_argument("--max-period", type=int, default=2000)

    p_verify = sub.add_parser(
        "verify-acta",
        help="recompute an acta's anchors from a manifest, naming every field"
             " that diverges")
    p_verify.add_argument("manifest", help="path to the manifest JSON")
    p_verify.add_argument("--source", default="",
                          help="the live acta: measure its blob, its lines and the"
                               " tables the manifest cites (read by AST, never imported)")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except ValueError as exc:
        print(f"BAD_REQUEST: {exc}", file=sys.stderr)
        return RC_BAD_REQUEST
    if args.self_test:
        return _self_test()
    if args.command is None:
        parser.print_help()
        return RC_BAD_REQUEST
    try:
        if args.command == "ceiling":
            return _print_ceiling(args.period, args.json)
        if args.command == "solve":
            return _print_solve(args)
        if args.command == "drift":
            return _print_drift(args)
        if args.command == "ladder":
            return _print_ladder(args.period, args.json)
        if args.command == "census":
            return _print_census(args.low, args.high, args.json)
        if args.command == "designs":
            return _print_designs(args.period, args.json)
        if args.command == "price":
            return _print_price(args.period, args.json)
        if args.command == "prices":
            return _print_price_sweep(args.low, args.high, args.json)
        if args.command == "frontier":
            return _print_frontier(args.max_period,
                                   args.universe == "reduced", args.json)
        if args.command == "blocks":
            return _print_blocks(args.limit, args.json)
        if args.command == "dominance":
            return _print_dominance(args.max_period, args.json)
        if args.command == "strong-lemma":
            return _print_strong_lemma(args.max_primes, args.json)
        if args.command == "sign-theorem":
            return _print_sign_theorem(args.max_period, args.json)
        if args.command == "verify-acta":
            return _print_verify_acta(args.manifest, args.json, args.source)
    except MaskTableTooLarge as exc:
        print(f"BEYOND_TABLE_BUDGET: {exc}", file=sys.stderr)
        return RC_BEYOND_TABLE_BUDGET
    except ValueError as exc:
        print(f"BAD_REQUEST: {exc}", file=sys.stderr)
        return RC_BAD_REQUEST
    parser.print_help()
    return RC_BAD_REQUEST


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
