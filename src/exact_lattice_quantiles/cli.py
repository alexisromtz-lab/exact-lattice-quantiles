"""Command line: the same numbers the library computes, printable and testable.

Exit codes are part of the contract:

``0``  the request was answered
``1``  the answer is no: the margin requested is above the instance ceiling --
       impossible at any price, and the ceiling is printed so the caller learns
       the number -- or an ``acta`` anchor does not reproduce
``2``  the margin requested exists but costs more points than the budget
``3``  the request itself is malformed
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction

from .ceiling import (InfeasibleAboveInstanceCeiling, ceiling_by_product,
                      ceiling_by_scan, ceiling_census, instance_ceiling)
from .acta import ManifestError, render, shipped_manifest, verify_file
from .coverage import (coverage_census, drift_law, drift_radius,
                       drift_tolerance, witness_margin)
from .lattice import universe
from .pricing import (block_catalogue, block_of, block_saving,
                      exactness_frontier, price_exact,
                      price_form_census, price_form_sweep)
from .solve import InfeasibleWithinBudget, ladder, solve_lattice_schedule

RC_OK = 0
RC_ABOVE_CEILING = 1
RC_OVER_BUDGET = 2
RC_BAD_REQUEST = 3


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
    if args.json:
        print(json.dumps(schedule.as_dict(), indent=2, sort_keys=True))
        return RC_OK
    print(f"period           {schedule.period}")
    print(f"design           {list(schedule.design)}")
    print(f"cost             {schedule.cost} points (budget {schedule.max_cost})")
    print(f"design margin    {schedule.design_margin}")
    print(f"sections ({len(schedule.probes)})  {[str(q) for q in schedule.probes]}")
    print(f"witness margin   {schedule.witness_margin}")
    print(f"robust margin    {schedule.robust_margin}")
    print(f"probes needed    {schedule.probes_needed}")
    print(f"ceiling          {schedule.ceiling.margin} "
          f"(critical divisor {schedule.ceiling.critical_divisor})")
    print(f"verdict          {schedule.verdict}")
    return RC_OK


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
    U = universe(schedule.period)
    Q = schedule.probes
    curve = []
    for share in (Fraction(1), Fraction(3, 4), Fraction(1, 2), Fraction(1, 4), Fraction(0)):
        kept = schedule.witness_margin * share
        tol = min(drift_tolerance(schedule.design, Q, U, kept).values()) if Q else Fraction(0)
        curve.append((kept, tol))
    coverage = min(drift_radius(schedule.design, Q, U).values()) if Q else Fraction(0)
    edge = (min(drift_law(schedule.design, Q, U).values(),
                key=lambda law: law.critical_margin) if Q else None)
    if args.json:
        print(json.dumps({
            "period": schedule.period,
            "design": list(schedule.design),
            "cost": schedule.cost,
            "probes": [str(q) for q in Q],
            "witness_margin": str(schedule.witness_margin),
            "coverage_drift": str(coverage),
            "critical_margin": str(edge.critical_margin) if edge else None,
            "cliff_at_critical_margin": str(edge.cliff) if edge else None,
            "kept_margin_to_drift": [[str(v), str(t)] for v, t in curve],
        }, indent=2, sort_keys=True))
        return RC_OK
    print(f"period           {schedule.period}")
    print(f"design           {list(schedule.design)}")
    print(f"cost             {schedule.cost} points (budget {schedule.max_cost})")
    print(f"sections ({len(Q)})  {[str(q) for q in Q]}")
    print(f"witness margin   {schedule.witness_margin}")
    print(f"coverage drift   {coverage}  (captured, margin not kept)")
    if edge is not None:
        print(f"linear up to     {edge.critical_margin}  "
              f"(one-for-one price, then a drop of {edge.cliff})")
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


def _print_verify_acta(manifest: str, as_json: bool) -> int:
    """Recompute the anchors of a manifest: the boundary is a gate, not a memory."""
    try:
        checks = verify_file(manifest)
    except ManifestError as exc:
        print(f"BAD_REQUEST: {exc}", file=sys.stderr)
        return RC_BAD_REQUEST
    divergent = [check for check in checks if not check.ok]
    anchors = len({check.anchor for check in checks})
    if as_json:
        print(json.dumps({
            "manifest": str(manifest),
            "anchors": anchors,
            "checks": len(checks),
            "divergent": len(divergent),
            "results": [{"anchor": c.anchor, "kind": c.kind, "field": c.field,
                         "manifest": render(c.expected),
                         "package": render(c.measured), "ok": c.ok}
                        for c in checks],
        }, indent=2, sort_keys=True))
    else:
        for check in checks:
            print(check)
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
    check("block of 270 is {54, 90, 135}", block_of(270)[0], (54, 90, 135))
    check("block of 270 saves 4", block_saving(270), 4)
    check("block of 30 saves nothing", block_saving(30), 0)
    check("every divisor universe misprices first at 60",
          exactness_frontier(2000, reduced=False)["period"], 60)
    check("every carrier universe misprices first at 1080",
          exactness_frontier(2000)["period"], 1080)
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

    p_verify = sub.add_parser(
        "verify-acta",
        help="recompute an acta's anchors from a manifest, naming every field"
             " that diverges")
    p_verify.add_argument("manifest", help="path to the manifest JSON")
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
        if args.command == "verify-acta":
            return _print_verify_acta(args.manifest, args.json)
    except ValueError as exc:
        print(f"BAD_REQUEST: {exc}", file=sys.stderr)
        return RC_BAD_REQUEST
    parser.print_help()
    return RC_BAD_REQUEST


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
