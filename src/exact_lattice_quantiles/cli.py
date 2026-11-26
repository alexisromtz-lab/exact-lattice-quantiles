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
from itertools import chain, combinations
from math import lcm

from .ceiling import (InfeasibleAboveInstanceCeiling, ceiling_by_product,
                      ceiling_by_scan, ceiling_census, instance_ceiling)
from .acta import (ManifestError, SourceError, inspect_source, load_manifest,
                   parse_tables, render, shipped_manifest, verify_file,
                   verify_manifest)
from .coverage import (coverage_census, coverage_kernel, design_law, drift_law,
                       drift_radius, drift_tolerance, guaranteed_margin,
                       margin_under_drift, minimum_cost_design,
                       witness_margin)
from .arithmetic import primes_of
from .lattice import (TABLE_BUDGET, MaskTableTooLarge, cell_widths,
                      grids, mask_table, period_of, universe)
from .pricing import (antichain_saving_ceiling, arity_ladder_census, band_slack,
                      band_window,
                      block_catalogue,
                      block_advantage, block_defect, block_dominance_census,
                      block_fee, block_fee_closed, block_free,
                      block_of, block_saving,
                      carrier_block_of, carrier_block_saving, cofactor_block,
                      excess_parity, exactness_frontier, extremal_census,
                      fee_frontier_by_blocks, fee_frontier_census,
                      heaviest_antichain,
                      odd_band_census,
                      pair_antichain_of,
                      pair_antichain_saving, price_exact,
                      price_form_census, price_form_sweep, prime_excess,
                      reciprocal_filter_census, reciprocal_prime_sum,
                      setup_fee_frontier, sign_regime, sign_theorem_census,
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


def _print_period_of(durations: list[int], as_json: bool) -> int:
    """The period a set of service lengths generates, and its exact ceiling.

    A system whose lengths are ``d1 ... dn`` repeats only every
    ``lcm(d1 ... dn)``, and its cadences are the divisors of that period above
    one.  This is the one step between a schedule somebody already runs and the
    lattice the rest of the commands price.
    """
    if not durations or any(d < 1 for d in durations):
        print("period: one or more positive durations are required",
              file=sys.stderr)
        return RC_BAD_REQUEST
    period = period_of(durations)
    c = instance_ceiling(period)
    per_period = {d: period // d for d in durations}
    if as_json:
        print(json.dumps({
            "durations": list(durations),
            "period": period,
            "grids": len(grids(period)),
            "universe": len(c.universe),
            "margin": str(c.margin),
            "critical_divisor": c.critical_divisor,
            "radical_bound_hit": c.radical_bound_hit,
            "ceiling_in_units": str(c.margin * period),
            "services_per_period": {str(k): v for k, v in per_period.items()},
        }, indent=2, sort_keys=True))
        return RC_OK
    print(f"durations         {', '.join(str(d) for d in durations)}")
    print(f"period            {period}  (lcm of the durations)")
    print(f"grids             {len(grids(period))} divisors of {period} above one")
    print(f"universe U(T)     {len(c.universe)} carriers of the even grids")
    print(f"ceiling nu(T)     {c.margin}  = {c.margin * period} of the period")
    print(f"critical divisor  {c.critical_divisor}")
    print(f"radical bound hit {c.radical_bound_hit}")
    print("services per period  "
          + ", ".join(f"{d}: {n}" for d, n in per_period.items()))
    return RC_OK


def _schedule_or_rc(period: int, max_cost: int, probes: int,
                    required: Fraction, as_json: bool
                    ) -> tuple[object, int | None]:
    """The schedule for a margin requirement, or the exit code of no schedule.

    Three named impossibilities -- above the instance ceiling, over the budget,
    above the table budget -- are reported here once instead of in every
    command that starts from a requirement, so a new command inherits the
    codes instead of re-deriving them.
    """
    try:
        return solve_lattice_schedule(period, max_cost, probes, required), None
    except InfeasibleAboveInstanceCeiling as exc:
        if as_json:
            print(json.dumps({
                "verdict": "INFEASIBLE_ABOVE_INSTANCE_CEILING",
                "period": exc.period,
                "needed": str(exc.needed),
                "ceiling": str(exc.ceiling),
            }, indent=2, sort_keys=True))
        else:
            print(f"INFEASIBLE_ABOVE_INSTANCE_CEILING: {exc}")
        return None, RC_ABOVE_CEILING
    except InfeasibleWithinBudget as exc:
        if as_json:
            print(json.dumps({
                "verdict": "OVER_BUDGET",
                "period": exc.period,
                "needed": str(exc.needed),
                "cheapest": exc.cheapest,
                "max_cost": exc.max_cost,
            }, indent=2, sort_keys=True))
        else:
            print(f"OVER_BUDGET: {exc}")
        return None, RC_OVER_BUDGET
    except MaskTableTooLarge as exc:
        if as_json:
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
        return None, RC_BEYOND_TABLE_BUDGET


def _print_solve(args: argparse.Namespace) -> int:
    schedule, rc = _schedule_or_rc(
        args.period, args.max_cost, args.probes,
        Fraction(0) if args.min_margin is None else args.min_margin, args.json)
    if rc is not None:
        return rc
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


def _print_plan(args: argparse.Namespace) -> int:
    """The purchase a service level is: keep a margin *under* a common drift.

    ``solve`` takes a margin and reports what the sections hold with no drift.
    A service level is the other shape -- a margin that has to survive a drift
    -- and its price is the price of ``keep + delta``: the drift is spent out
    of the witness margin, so it is one purchase and not two searches.  The
    answer is checked against :func:`guaranteed_margin` before it is printed;
    if the two ever disagree the command reports ``NO_MEDIDO`` instead of a
    number, because a service level read from a stale law is not a reading.
    """
    delta = Fraction(0) if args.delta is None else Fraction(args.delta)
    keep = Fraction(args.keep)
    if keep < 0:
        raise ValueError(f"a margin to keep cannot be negative, got {keep}")
    if delta < 0:
        raise ValueError(f"a drift cannot be negative, got {delta}")
    required = keep + delta
    schedule, rc = _schedule_or_rc(args.period, args.max_cost, args.probes,
                                   required, args.json)
    if rc is not None:
        return rc
    U = universe(schedule.period)
    kept = None
    if schedule.probes and schedule.witness_margin is not None:
        kept = guaranteed_margin(schedule.design, schedule.probes, U, delta)
    short = schedule.verdict != OK or kept is None or kept < keep
    if args.json:
        print(json.dumps({
            "verdict": schedule.verdict,
            "period": schedule.period,
            "keep": str(keep),
            "delta": str(delta),
            "required": str(required),
            "design": list(schedule.design),
            "cost": schedule.cost,
            "max_cost": schedule.max_cost,
            "probes": [str(q) for q in schedule.probes],
            "probes_needed": schedule.probes_needed,
            "witness_margin": None if schedule.witness_margin is None
            else str(schedule.witness_margin),
            "kept_margin": None if kept is None else str(kept),
            "meets_service_level": not short,
            "ceiling": str(schedule.ceiling.margin),
        }, indent=2, sort_keys=True))
        return RC_NO_MEDIDO if short else RC_OK
    print(f"period           {schedule.period}")
    print(f"service level    keep {keep} under a drift of {delta}"
          f"  ->  buy margin {required}")
    print(f"design           {list(schedule.design)}")
    print(f"cost             {schedule.cost} points (budget {schedule.max_cost})")
    print(f"sections ({len(schedule.probes)})  {[str(q) for q in schedule.probes]}")
    print(f"witness margin   {schedule.witness_margin}")
    if kept is not None:
        print(f"under drift      {delta}  ->  keeps margin {kept}")
    print(f"verdict          {schedule.verdict}")
    if short:
        held = "no sections realise it" if kept is None else f"holds {kept}"
        print(f"shortfall        {held} against a service level of {keep}"
              f"  (ceiling {schedule.ceiling.margin})")
    return RC_NO_MEDIDO if short else RC_OK


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
    kept_by_divisor: dict[int, Fraction] = {}
    kept = Fraction(0)
    if args.delta is not None and Q:
        kept_by_divisor = margin_under_drift(schedule.design, Q, U, args.delta)
        kept = min(kept_by_divisor.values())
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
            "delta": None if args.delta is None else str(args.delta),
            "guaranteed_margin": None if args.delta is None else str(kept),
            "divisor_guaranteed_margin": None if args.delta is None else
            {str(d): str(g) for d, g in sorted(kept_by_divisor.items())},
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
    for share, tol in curve:
        print(f"  {str(share):<15} {tol}")
    if args.delta is not None:
        print(f"under drift      {args.delta}  ->  keeps margin {kept}")
        print(f"budget           one currency: keeping {kept} under that drift "
              f"prices the same as keeping {kept + args.delta} with none")
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
    print(line("law, periods with w >= 2", str(census["law_periods"]),
               "excess and slack share their parity"))
    print(line("  violations", str(len(census["law_violations"])),
               "none: the law is an identity, not a census"))
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


def _print_extremal_gap(max_period: int, as_json: bool) -> int:
    """The extremal reading: is the block not merely allowed, but the maximum?

    ``RC_ABOVE_CEILING`` when an antichain outweighs the block -- the package's
    own claim refuted by its own census -- and ``RC_NO_MEDIDO`` when the horizon
    cannot reach the band, because the band is the only region whose warrant is a
    measurement and a census that never enters it would be reporting the proved
    half as if it were the whole.
    """
    census = extremal_census(max_period)
    broken = bool(census["violations"])
    band = census["band_rows"]
    if as_json:
        print(json.dumps(census, indent=2, sort_keys=True))
    else:
        regions = census["regions"]

        def line(label: str, value: str, note: str = "") -> str:
            return f"  {label:<32}{value:>6}  {note}"

        print("extremal reading: the heaviest antichain is the drop-one block")
        print(line(f"periods 2..{max_period}", str(census["periods"]),
                   "with at least two distinct primes"))
        for name in ("proven-zero", "block-saves", "band"):
            row = regions.get(name)
            if row is None:
                continue
            print(line(name, str(row["rows"]),
                       f"gap 0 on {row['gap_zero']}, warrant {row['warrant']}"))
        print(line("the block is the maximum",
                   str(sum(r["block_attains"] for r in regions.values())),
                   "of those periods, on every one"))
        if band:
            print(line("band members",
                       " ".join(str(r["period"]) for r in band),
                       "gap 0, one unit of saving away"))
        tightest = census["tightest"]
        if tightest is not None:
            print(line("closest to a saving", str(tightest["period"]),
                       f"weight {tightest['weight']} on a bound of "
                       f"{tightest['bound']}, {tightest['miss']} points short"))
        for row in census["violations"]:
            print(f"  VIOLATION  {row['period']}: the antichain "
                  f"{tuple(row['witness'])} weighs {row['weight']} against the "
                  f"block's {row['bound']}")
        print(f"  theorem  {census['theorem']}")
        print(f"  frontier {census['frontier']}")
        print(f"  verdict  {census['verdict']}")
    if broken:
        return RC_ABOVE_CEILING
    if not band:
        print(f"NO_MEDIDO: no band member below {max_period}, so the one region "
              f"whose warrant is a measurement was not measured", file=sys.stderr)
        return RC_NO_MEDIDO
    return RC_OK


def _print_setup_fee(period: int, reduced: bool, max_setup: int,
                     as_json: bool) -> int:
    """The closed form's unstated hypothesis: no fixed cost per grid.

    ``RC_NO_MEDIDO`` when no fee inside the horizon breaks the formula, because
    a frontier beyond the horizon is not an absent frontier -- the caller asked
    where the formula stops being the price and the answer was not reached.
    """
    row = setup_fee_frontier(period, reduced=reduced, max_setup=max_setup)
    if as_json:
        payload = {k: v for k, v in row.items() if k != "rows"}
        payload["rows"] = [r.as_dict() for r in row["rows"]]
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        def line(label: str, value: str, note: str = "") -> str:
            return f"  {label:<22}{value:<14}{note}"

        print("setup fee: the closed form assumes a grid costs grid // 2 and nothing else")
        print(line("period", str(period),
                   f"{row['universe']} universe, |U| {row['universe_size']}, "
                   f"{row['consumers']} requests"))
        print(line("fee", "discrepancies", "pairs  witness"))
        for r in row["rows"]:
            pairs = str(r.pair_failures)
            witness = (f"{tuple(r.witness)} by {r.witness_gap}"
                       if r.witness else "-")
            print(line(f"  {r.setup}", str(r.discrepancies),
                       f"{pairs:>5}  {witness}"))
        frontier = row["frontier"]
        print(line("frontier", "none" if frontier is None else str(frontier),
                   "first fee per grid that breaks the closed form"))
        pairs_at = row["pair_frontier"]
        print(line("pairs",
                   "none" if pairs_at is None else str(pairs_at),
                   "first fee at which a pair is what fails"))
        print(f"  theorem  {row['theorem']}")
        print(f"  verdict  {row['verdict']}")
    if row["frontier"] is None:
        print(f"NO_MEDIDO: no fee up to {max_setup} breaks the closed form on "
              f"period {period} ({row['universe']}); the frontier is beyond the "
              f"horizon, not absent", file=sys.stderr)
        return RC_NO_MEDIDO
    return RC_OK


def _print_block_fee(parts: tuple[int, ...], as_json: bool) -> int:
    """The fee of one block, read off its defect instead of searched for.

    ``RC_BAD_REQUEST`` when the argument is not a block at all: one grid has no
    defect and no fee, and answering ``0`` would read as "free".
    """
    if len(parts) < 2:
        print("NO BLOCK: a block needs at least two grids; one grid is never "
              "cheaper than itself", file=sys.stderr)
        return RC_BAD_REQUEST
    lcm_all = lcm(*parts)
    defect = block_defect(parts)
    odd_parts = sum(1 for p in parts if p % 2)
    twice_b = odd_parts - defect if lcm_all % 2 == 0 else len(parts) - 1 - defect
    fee = block_fee_closed(parts)
    if as_json:
        print(json.dumps({"parts": list(parts), "lcm": lcm_all,
                          "sum": sum(parts), "defect": defect,
                          "odd_parts": odd_parts, "twice_b": twice_b,
                          "fee": fee, "fee_by_advantage": block_fee(parts)},
                         indent=2, sort_keys=True))
    else:
        def line(label: str, value: str, note: str = "") -> str:
            return f"  {label:<24}{value:<16}{note}"

        print(f"block fee: {tuple(parts)}")
        print(line("lcm", str(lcm_all), "one grid covering every part"))
        print(line("sum", str(sum(parts)), "one grid per part"))
        print(line("defect", str(defect), "sum - lcm: zero is a perfect block"))
        print(line("odd parts", str(odd_parts),
                   "k - defect, the lcm is even" if lcm_all % 2 == 0 else
                   f"{len(parts)} - 1 - defect, the lcm is odd"))
        print(line("twice B", str(twice_b),
                   "twice the advantage at a zero fee"))
        print(line("fee", str(fee),
                   f"floor(twice B / {2 * (len(parts) - 1)}) + 1"))
        print(line("by advantage", str(block_fee(parts)),
                   "the same fee, inverted the slow way"))
        verdict = ("free: the block already beats its parts with no fee at all"
                   if fee == 0 else
                   f"breaks the closed form from a fee of {fee} up")
        print(f"  verdict  {verdict}")
    return RC_OK


def _print_fee_census(max_period: int, reduced: bool, arity: int,
                      verify_universe: int, as_json: bool) -> int:
    """The fee frontier is a law: two parts cost 2, and a fee of 1 needs a triple.

    ``RC_ABOVE_CEILING`` when the law misprices a period -- the package's own
    census refuting a claim the package makes. ``RC_NO_MEDIDO`` when no period
    was verified, because a law checked nowhere is not a law measured.
    """
    row = fee_frontier_census(max_period=max_period, reduced=reduced,
                              arity=arity, verify_universe=verify_universe)
    if as_json:
        payload = {k: v for k, v in row.items() if k != "rows"}
        payload["rows"] = list(row["rows"])
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        def line(label: str, value: str, note: str = "") -> str:
            return f"  {label:<24}{value:<16}{note}"

        print("fee census: the fee that breaks the closed form, period by period")
        print(line("horizon", str(max_period),
                   f"{row['universe']} universe, blocks of at most {arity} parts"))
        print(line("periods", str(row["periods"]),
                   f"{row['verified']} verified against every antichain "
                   f"(|U| <= {verify_universe})"))
        print(line("mismatches", str(len(row["mismatches"])),
                   "periods where no block of that size attains the frontier"))
        for m in row["mismatches"][:10]:
            print(line(f"  period {m['period']}", str(m['exhaustive']),
                       f"law says {m['law']} by {m['exhaustive_witness']}"))
        print(line("chains", str(len(row["chains"])),
                   "no antichain at all: no fee breaks these"))
        print(line("at zero", str(len(row["at_zero"])),
                   "none in this horizon" if not row["at_zero"] else
                   f"already broken at no fee {row['at_zero'][:7]}"))
        print(line("at one", str(len(row["at_one"])),
                   "break first at a fee of 1"))
        largest = row["largest"]
        print(line("widest frontier",
                   "none" if largest is None else str(largest["frontier"]),
                   "every period here is a chain" if largest is None else
                   f"period {largest['period']} on {largest['witness']}"))
        sizes = row["verified_witness_sizes"]
        print(line("witness sizes",
                   ", ".join(f"{k}: {v}" for k, v in sorted(sizes.items())),
                   "on the verified periods"))
        if row["arity"] >= 3:
            print(line("pairs miss", str(len(row["pair_misses"])),
                       "a triple beats every pair"))
            print(line("pair over", str(len(row["pair_over"])),
                       "of those, above a fee of one: the law refuted if not 0"))
        print(f"  pair law {row['pair_law']}")
        print(f"  verdict  {row['verdict']}")
    if row["mismatches"]:
        print(f"REFUTED: the block law misprices {len(row['mismatches'])} of "
              f"{row['periods']} periods up to {max_period}",
              file=sys.stderr)
        return RC_ABOVE_CEILING
    if row["pair_over"]:
        print(f"REFUTED: {len(row['pair_over'])} periods reach a fee of two or "
              f"more only through a triple, so a pair does not suffice",
              file=sys.stderr)
        return RC_ABOVE_CEILING
    if row["verified"] == 0:
        print(f"NO_MEDIDO: no period up to {max_period} had at most "
              f"{verify_universe} grids, so the law was checked nowhere",
              file=sys.stderr)
        return RC_NO_MEDIDO
    return RC_OK


def _print_arity_ladder(max_period: int, reduced: bool, arities: str,
                        exhaustive_until: int, as_json: bool) -> int:
    """How many parts the frontier needs, and where the free rung opens.

    ``RC_ABOVE_CEILING`` when a rung above the top one still earns -- not
    possible for a fee, since zero is the floor, so seeing it refutes the ladder
    -- and ``RC_NO_MEDIDO`` when no period was verified against every antichain.
    A horizon with no free block is *not* an error: it is the measured floor of
    that universe, and the verdict says so in words.
    """
    rungs = tuple(int(part) for part in arities.replace(",", " ").split())
    row = arity_ladder_census(max_period=max_period, reduced=reduced,
                              arities=rungs or (2, 3, 4, 5),
                              exhaustive_until=exhaustive_until)
    if as_json:
        payload = {k: v for k, v in row.items() if k != "exhaustive"}
        payload["exhaustive"] = {k: v for k, v in row["exhaustive"].items()
                                 if k != "rows"}
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        def line(label: str, value: str, note: str = "") -> str:
            return f"  {label:<22}{value:<18}{note}"

        print("arity ladder: how many parts the fee frontier needs")
        print(line("horizon", str(row["max_period"]),
                   f"{row['universe']} universe, rungs {row['arities']}"))
        for step, periods in row["earnings"].items():
            first = row["first"].get(step)
            print(line(f"  {step}", str(len(periods)),
                       "never earns" if first is None else
                       f"first {first['period']} on {first['witness']} "
                       f"({first['from_fee']} -> {first['to_fee']}, "
                       f"{first['universe_size']} grids)"))
        free = row["first_free"]
        print(line("first free block",
                   "none in horizon" if free is None else str(free["period"]),
                   "every period here keeps a fee of at least one" if free is None
                   else f"{free['witness']}, already free with {free['parts']} "
                        f"part(s)"))
        print(line("periods at zero", str(len(row["at_floor"])),
                   "none in this horizon" if not row["at_floor"]
                   else f"{row['at_floor'][:7]}"))
        print(line("verified", str(row["exhaustive"]["verified"]),
                   f"against every antichain (|U| <= {exhaustive_until})"))
        print(line("mismatches", str(len(row["exhaustive"]["mismatches"])),
                   "the ladder against the exhaustive frontier"))
        print(f"  verdict  {row['verdict']}")
    if row["exhaustive"]["mismatches"]:
        print("REFUTED: the ladder is not the frontier on "
              f"{len(row['exhaustive']['mismatches'])} period(s)",
              file=sys.stderr)
        return RC_ABOVE_CEILING
    if row["exhaustive"]["verified"] == 0:
        print(f"NO_MEDIDO: no period up to {max_period} had at most "
              f"{exhaustive_until} grids, so the ladder was checked nowhere",
              file=sys.stderr)
        return RC_NO_MEDIDO
    return RC_OK


def _print_odd_band(max_omega: int, max_prime: int, slack: int | None,
                    as_json: bool) -> int:
    """The odd half of the band: a window through Giuga's equation, measured.

    ``RC_OK`` when the horizon was walked and nothing is in the band,
    ``RC_ABOVE_CEILING`` when something is -- the package's own claim refuted by
    its own census -- and ``RC_NO_MEDIDO`` when the horizon cannot reach a band
    at all, because reporting "none" from a search that could not cross one
    would be the dangerous green.
    """
    census = odd_band_census(max_omega, max_prime, slack)
    if as_json:
        print(json.dumps(census, indent=2, sort_keys=True))
    else:
        def line(label: str, value: str, note: str = "") -> str:
            return f"  {label:<34}{value:>8}  {note}"

        print("odd band: 1 <= period * (sum 1/p - 1) <= omega - 1")
        print(line("odd primes available", str(census["odd_primes_available"]),
                   f"at most {max_prime}"))
        print(line("prefixes walked", str(census["prefixes"]),
                   "minimal feasible sets only"))
        print(line("threshold", str(census["threshold"]), "omega - 1 by law"))
        floor = census["floor"]
        if floor["omega"] is None:
            print(line("feasibility floor", "not reached", floor["below"]))
        else:
            print(line("feasibility floor", f"omega = {floor['omega']}",
                       f"product {floor['product']} clears the window by "
                       f"{floor['excess']}"))
        print(line("members in the band", str(len(census["candidates"]))))
        for member in census["candidates"][:8]:
            print(f"  {' '.join(str(p) for p in member['primes']):<34}"
                  f"{member['excess']:>8}  excess against slack "
                  f"{member['slack']}")
        margin = census["margin"]
        nearest = margin["nearest"]
        if nearest is None:
            print(line("nearest miss", "none", margin["verdict"]))
        else:
            print(line("nearest miss", str(nearest["excess"]),
                       f"clears the band by {nearest['distance']} at "
                       f"omega = {nearest['omega']}"))
            print(f"  {' '.join(str(p) for p in nearest['primes']):<34}"
                  f"{nearest['excess']:>8}  the closest approach from above")
            for key, row in margin["by_omega"].items():
                print(f"  omega {key:<3}{'':<26}{row['excess']:>8}  by "
                      f"{row['distance']} over {row['threshold']}")
        parity = census["parity"]
        print(line("parity of the excess", "consistent" if parity["consistent"]
                   else "BROKEN", parity["rule"]))
        print(f"  corner: {parity['corner']}")
        print(f"  boundary: {census['boundary']}")
        print(f"  verdict: {census['verdict']}")
    if census["status"] == "below-floor":
        return RC_NO_MEDIDO
    return RC_ABOVE_CEILING if census["candidates"] else RC_OK


def _print_universes(max_period: int, as_json: bool) -> int:
    """The reciprocal sum against both universes: which one decides a period.

    ``RC_OK`` when the census ran.  A period below one that breaks over some
    universe would mean the proof does not reach the census, so that exits
    ``RC_ABOVE_CEILING`` -- the package's own claim refuted by its own table.
    """
    census = reciprocal_filter_census(max_period)
    if as_json:
        print(json.dumps(census, indent=2, sort_keys=True))
    else:
        def line(label: str, value: str, note: str = "") -> str:
            return f"  {label:<34}{value:>8}  {note}"

        print(f"reciprocal sum against universe, even periods to {max_period}")
        print(line("periods", str(census["periods"])))
        print(line("sum 1/p at most one", str(census["proved"]),
                   "no block saves in either universe"))
        print(line("above one", str(census["band"]),
                   "the proof says nothing here"))
        print(line("  clean over U(T)", str(census["band_clean_reduced"]),
                   "the filter folds the band"))
        print(line("  breaking over D(T)", str(census["band_breaking_divisors"]),
                   "same period, same sum, other universe"))
        widest = census["widest"]
        if widest is not None:
            print(line("widest gap", f"T={widest['period']}",
                       f"excess {widest['excess']}  U(T) {widest['reduced']}"
                       f"  D(T) {widest['divisors']}"))
        print(f"  verdict: {census['verdict']}")
    if census["proved_without_break"] < census["proved"]:
        return RC_ABOVE_CEILING
    return RC_OK


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
    clinic = period_of((15, 20, 30, 45, 60))
    check("a clinic of 15/20/30/45/60 runs over period 180", clinic, 180)
    check("its ceiling is the ceiling of 180, not of 60",
          instance_ceiling(clinic).margin, Fraction(1, 3))
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
    check("drift spends the margin one-for-one",
          guaranteed_margin(G, Q, U, Fraction(1, 60)),
          s.witness_margin - Fraction(1, 60))
    check("and the inverse round-trips",
          min(drift_tolerance(G, Q, U,
                              guaranteed_margin(G, Q, U, Fraction(1, 60))).values()),
          Fraction(1, 60))
    check("a margin past the drift is not kept",
          min(drift_tolerance(G, Q, U,
                              guaranteed_margin(G, Q, U, Fraction(1, 60))
                              + Fraction(1, 1000)).values()) < Fraction(1, 60), True)
    check("no drift keeps the whole witness margin",
          guaranteed_margin(G, Q, U), s.witness_margin)
    check("drift wider than the margin keeps none",
          guaranteed_margin(G, Q, U, s.witness_margin * 2), Fraction(0))
    service, drift = Fraction(1, 30), Fraction(1, 60)
    plan = solve_lattice_schedule(180, max_cost=64, k_probes=5,
                                  min_margin=service + drift)
    check("a service level is bought as the margin plus the drift",
          plan.design, (12, 18, 20, 45))
    check("and under that drift it still keeps the service level",
          guaranteed_margin(plan.design, plan.probes, U, drift) >= service, True)
    try:
        solve_lattice_schedule(180, max_cost=10**6, k_probes=8,
                               min_margin=instance_ceiling(180).margin + drift)
    except InfeasibleAboveInstanceCeiling as exc:
        check("a service level above the ceiling minus the drift is refused",
              exc.ceiling, Fraction(1, 3))
    else:
        check("a service level above the ceiling minus the drift is refused",
              "no exception", "exception")
    check("block of 270 is {54, 90, 135}", block_of(270)[0], (54, 90, 135))
    check("block of 270 saves 4", block_saving(270), 4)
    check("block of 30 saves nothing", block_saving(30), 0)
    check("every divisor universe misprices first at 60",
          exactness_frontier(2000, reduced=False)["period"], 60)
    check("every carrier universe misprices first at 1080",
          exactness_frontier(2000)["period"], 1080)
    g2520, u2520 = grids(2520), universe(2520)
    w2520 = cell_widths(g2520, u2520)
    kept, _ = coverage_kernel(g2520, u2520, w2520, Fraction(1, 6))
    check("the kernel of 2520 at 1/6 is 20 of 32 divisors", len(kept), 20)
    check("and its cheapest design costs 888",
          minimum_cost_design(g2520, u2520, w2520, Fraction(1, 6), 2520)[0], 888)
    check("a period the table used to refuse now reaches its ceiling",
          minimum_cost_design(grids(1440), universe(1440),
                              cell_widths(grids(1440), universe(1440)),
                              Fraction(4, 15), 1440)[0], 735)
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
    heavy = heaviest_antichain(30)
    check("so the extremal reading of the lemma holds on the band too",
          (heavy["bound"], heavy["weight"], heavy["gap"]), (31, 31, 0))
    check("and the maximum is the block, not merely a number that ties",
          heavy["block_attains"], True)
    check("with the warrant named: the band is measured, never proved",
          heavy["warrant"], "measured")
    check("a period the lemma does prove carries a different warrant",
          heaviest_antichain(36)["warrant"], "theorem")
    checked = extremal_census(300)
    check("no antichain outweighs the block below 300",
          checked["violations"], [])
    check("and the band's only member below 300 is 30",
          [row["period"] for row in checked["band_rows"]], [30])
    check("the band misses a saving by one unit, the least a miss can be",
          (checked["band_rows"][0]["miss"],
           checked["band_rows"][0]["tightest_miss"]), (2, True))
    check("the miss is the distance from the top of the window",
          (heavy["miss"], heavy["slack"] - heavy["excess"] + 2), (2, 2))
    check("so tightness is m == slack and nothing else",
          (heavy["excess"], heavy["slack"], heavy["tightest_miss"]),
          (1, 1, True))
    check("and a miss is always even, because the weights share a parity",
          (heavy["weight"] % 2, (30 + heavy["slack"]) % 2, heavy["miss"] % 2),
          (1, 1, 0))
    check("a block that ties at a zero fee wins by its saved fees", 
          (block_advantage((6, 10, 15)), block_advantage((6, 10, 15), 1)),
          (0, 2))
    check("an incomparable pair never wins at a zero fee",
          max(block_advantage((a, b)) for a in (4, 6, 12, 20, 30, 60)
              for b in (4, 6, 12, 20, 30, 60)
              if a < b and b % a and a % b), -1)
    fee = setup_fee_frontier(180, max_setup=2)
    check("the closed form is the price until a fixed cost is charged",
          [r.discrepancies for r in fee["rows"]], [0, 32, 146])
    check("the first break is a triple, because folding saves k - 1 fees",
          (fee["frontier"], fee["witness"]), (1, (6, 10, 15)))
    check("and no pair explains that break: pairs need a fee of two",
          (fee["rows"][1].pair_failures, fee["pair_frontier"]), (0, 2))
    try:
        heaviest_antichain(7)
    except ValueError:
        check("a prime period has no antichain to weigh", "ValueError",
              "ValueError")
    census_300 = sign_theorem_census(300)
    check("below 300 the band has exactly one member",
          [row["period"] for row in census_300["band_rows"]], [30])
    check("and the lemma holds on every period it covers",
          census_300["lemma_violations"], [])
    check("and so does the theorem behind it",
          census_300["sign_violations"], [])
    check("the price at 1080 is 135, not 139",
          price_exact(1080, (54, 90, 135)), 135)
    check("the block's saving is the excess minus the slack, halved",
          block_saving(180), (prime_excess(180) - band_slack(180)) // 2)
    check("an odd period carries slack omega - 1, not one",
          (band_slack(45), band_slack(60), band_slack(36)), (1, 0, 0))
    check("the window above {3, 5} closes on 2 and lands on 30",
          (band_window(15, -7, 2), prime_excess(30)), ((2, 2), 1))
    check("so every period is characterised, the odd one as a window",
          [sign_regime(p)["characterisation"] for p in (45, 30, 36)],
          ["giuga-window", "giuga-equation", "no-member-possible"])
    check("no odd set of eight primes reaches a reciprocal sum of one",
          odd_band_census(8, 150)["status"], "below-floor")
    floor = odd_band_census(9, 30)["floor"]
    check("and the first one that does is nine primes, far above the window",
          (floor["omega"], floor["excess"]), (9, 108169298))
    walked = odd_band_census(10, 100)
    check("no odd period in the band with ten primes at most 100",
          (walked["candidates"], walked["prefixes"] > 0), ([], True))
    check("the same walk finds members the moment the threshold is lifted",
          len(odd_band_census(9, 30, slack=10**9)["candidates"]) > 0, True)
    check("the excess and the slack always share their parity",
          {(prime_excess(n) - band_slack(n)) % 2
           for n in range(6, 600) if len(primes_of(n)) >= 2}, {0})
    check("so 858 carries an odd excess and 36 an even one",
          (excess_parity(858), excess_parity(36)), (1, 0))
    check("and an odd period walks only the excesses its omega allows",
          {excess_parity(n) == (len(primes_of(n)) - 1) % 2
           for n in range(3, 600, 2) if len(primes_of(n)) >= 2}, {True})
    nearest = walked["margin"]["nearest"]
    check("the walk reports its nearest miss, not just its silence",
          nearest["excess"] > nearest["threshold"], True)
    check("and that miss sits where the congruence says it must",
          walked["parity"]["consistent"], True)
    check("a pair needs a fee of two and (2, 3) attains it",
          (block_fee((2, 3)), block_fee((6, 10, 15)),
           block_fee((4, 8))), (2, 1, 0))
    check("the fee is half the defect, divided by one less part -- not a search",
          sum(block_fee(parts) != block_fee_closed(parts)
              for period in range(2, 181, 2)
              for parts in chain.from_iterable(
                  combinations(universe(period), size) for size in (2, 3))), 0)
    check("a block one short of perfect pays a fee of one",
          (block_defect((6, 10, 15)), block_defect((6, 14, 21))), (1, -1))
    law = fee_frontier_census(max_period=200, verify_universe=10)
    check("the law holds on every period whose antichains were walked",
          (law["mismatches"],
           law["verified"] == sum(1 for r in law["rows"]
                                  if r["universe_size"] <= 10)), ([], True))
    check("a triple never beats every pair once the fee reaches two",
          law["pair_over"], ())
    check("and where a triple does beat them all, the fee is at most one",
          len(law["pair_misses"]) == len(law["at_zero"]) + len(law["at_one"]),
          True)
    check("a block whose parts are dearer than their lcm pays nothing at all",
          (block_free((12, 20, 30)), block_fee((12, 20, 30)),
           block_free((6, 10, 15)), block_fee((6, 10, 15))),
          (True, 0, False, 1))
    check("and no pair can ever get there: (2, 3) is the best two parts do",
          (block_free((2, 3)), block_fee((2, 3))), (False, 2))
    ladder = arity_ladder_census(max_period=500)
    check("the ladder earns a fee of one from a triple and of two from a pair",
          (ladder["earnings"]["2->3"][:2], ladder["first"]["2->3"]["witness"]),
          ((30, 42), (6, 10, 15)))
    check("the fourth part is free, and no fifth part buys anything",
          (ladder["first"]["3->4"]["witness"], ladder["earnings"]["4->5"]),
          ((42, 66, 154, 231), ()))
    check("the free rung opens at 60, and three parts are enough to reach it",
          (ladder["first_free"]["period"], ladder["first_free"]["parts"]),
          (60, 3))
    check("the ladder is the frontier wherever it could be checked",
          (ladder["exhaustive"]["mismatches"],
           ladder["exhaustive"]["verified"] > 0), ((), True))
    check("over the carriers no period here reaches zero, and the census says so",
          arity_ladder_census(max_period=500,
                              reduced=True)["first_free"], None)
    check("Giuga's equation is one defect short of free",
          (block_defect(cofactor_block((2, 3, 7, 41))),
           block_fee(cofactor_block((2, 3, 7, 41))),
           block_free(cofactor_block((2, 3, 7, 11)))), (1, 1, True))
    check("and a chain answers with no frontier, not with a zero",
          (fee_frontier_by_blocks(4)["chain"],
           fee_frontier_by_blocks(6)["frontier"],
           fee_frontier_by_blocks(180)["frontier"]), (True, 2, 1))
    universes = reciprocal_filter_census(400)
    check("every period below one is clean in both universes, not just one",
          (universes["proved"], universes["proved_without_break"]),
          (187, 187))
    check("the same period breaks over every divisor and not over the carriers",
          [row["divisors"] for row in universes["rows"] if row["period"] == 60],
          [1])
    check("and one period of the band breaks in neither universe",
          [row["period"] for row in universes["rows"]
           if row["above_one"]
           and row["reduced"] == 0 and row["divisors"] == 0], [30])
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

    p_period = sub.add_parser("period",
                              help="the period a set of service lengths generates")
    p_period.add_argument("durations", type=int, nargs="+")

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
    p_drift.add_argument(
        "--delta", type=_fraction, default=None,
        help="common drift to absorb; prints the margin the schedule still keeps")

    p_plan = sub.add_parser(
        "plan", help="cheapest design that keeps a margin under a common drift")
    p_plan.add_argument("period", type=int)
    p_plan.add_argument("--max-cost", type=int, required=True)
    p_plan.add_argument("--probes", type=int, required=True)
    p_plan.add_argument("--keep", type=_fraction, required=True,
                        help="margin the sections must still hold under the drift")
    p_plan.add_argument("--delta", type=_fraction, default=None,
                        help="common drift the schedule has to absorb")

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

    p_odd = sub.add_parser(
        "odd-band",
        help="walk every odd period the band can reach below a horizon: an odd "
             "member is a window through Giuga's equation")
    p_odd.add_argument("--max-omega", type=int, default=11,
                       help="most odd primes in a candidate period")
    p_odd.add_argument("--max-prime", type=int, default=150,
                       help="largest odd prime a candidate may use")
    p_odd.add_argument("--slack", type=int, default=None,
                       help="replace the per-member threshold omega - 1, which"
                            " is how the walk is shown not to be blind")

    p_extremal = sub.add_parser(
        "extremal-gap",
        help="read the lemma as an extremum: is the drop-one block the heaviest "
             "antichain, with the warrant of its region attached?")
    p_extremal.add_argument("--limit", type=int, default=2000)

    p_setup = sub.add_parser(
        "setup-fee",
        help="find the first fixed cost per grid at which the closed form stops "
             "being the price, and the request that breaks it")
    p_setup.add_argument("--period", type=int, default=180)
    p_setup.add_argument(
        "--universe", choices=("reduced", "divisors"), default="reduced",
        help="reduced carriers U(T), where the formula is exact at a zero fee, "
             "or every divisor D(T), where it is already broken")
    p_setup.add_argument("--max-setup", type=int, default=4,
                         help="largest fee per grid the walk reaches")

    p_block_fee = sub.add_parser(
        "block-fee",
        help="the fee of one block, read off its defect: half the distance "
             "between its lcm and its sum, over one less part")
    p_block_fee.add_argument("parts", nargs="+", type=int,
                             help="the grids folded into the block")

    p_fee = sub.add_parser(
        "fee-census",
        help="the fee at which the closed form breaks is always attained by a "
             "block of two or three parts")
    p_fee.add_argument("--max-period", type=int, default=2000)
    p_fee.add_argument(
        "--universe", choices=("reduced", "divisors"), default="reduced",
        help="reduced carriers U(T), where the formula is exact at a zero fee, "
             "or every divisor D(T), where it is already broken")
    p_fee.add_argument("--arity", type=int, default=3,
                       help="largest block the law is allowed to use")
    p_fee.add_argument("--verify-universe", type=int, default=14,
                       help="verify against every antichain when |U| is at most"
                            " this (the cost is exponential in |U|)")

    p_ladder = sub.add_parser(
        "arity-ladder",
        help="how many parts the fee frontier needs: the rungs 2 -> 3 -> 4, "
             "and where the free rung first opens")
    p_ladder.add_argument("--max-period", type=int, default=1000)
    p_ladder.add_argument(
        "--universe", choices=("reduced", "divisors"), default="divisors",
        help="every divisor D(T), where a free block appears at 60, or the "
             "carriers U(T), whose frontier never reached zero up to 1000")
    p_ladder.add_argument("--arities", default="2,3,4,5",
                          help="the rungs to climb, in order (2,3,4,5)")
    p_ladder.add_argument("--exhaustive-until", type=int, default=15,
                          help="check the ladder against every antichain when "
                               "|U| is at most this; 15 is what reaches 462")

    p_univ = sub.add_parser(
        "universes",
        help="census the reciprocal sum against both universes, period by "
             "period: the proof covers one column, the filter decides both")
    p_univ.add_argument("--max-period", type=int, default=400)

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
        if args.command == "period":
            return _print_period_of(args.durations, args.json)
        if args.command == "ceiling":
            return _print_ceiling(args.period, args.json)
        if args.command == "solve":
            return _print_solve(args)
        if args.command == "plan":
            return _print_plan(args)
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
        if args.command == "extremal-gap":
            return _print_extremal_gap(args.limit, args.json)
        if args.command == "setup-fee":
            return _print_setup_fee(args.period, args.universe == "reduced",
                                    args.max_setup, args.json)
        if args.command == "block-fee":
            return _print_block_fee(tuple(args.parts), args.json)
        if args.command == "fee-census":
            return _print_fee_census(args.max_period,
                                     args.universe == "reduced", args.arity,
                                     args.verify_universe, args.json)
        if args.command == "verify-acta":
            return _print_verify_acta(args.manifest, args.json, args.source)
        if args.command == "odd-band":
            return _print_odd_band(args.max_omega, args.max_prime,
                                   args.slack, args.json)
        if args.command == "arity-ladder":
            return _print_arity_ladder(args.max_period,
                                       args.universe == "reduced",
                                       args.arities, args.exhaustive_until,
                                       args.json)
        if args.command == "universes":
            return _print_universes(args.max_period, args.json)
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
