"""Command line: the same numbers the library computes, printable and testable.

Exit codes are part of the contract:

``0``  the request was answered
``1``  the margin requested is above the instance ceiling -- impossible at any
       price, and the ceiling is printed so the caller learns the number
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
from .coverage import coverage_census
from .pricing import price_form_census
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
        if args.command == "ladder":
            return _print_ladder(args.period, args.json)
        if args.command == "census":
            return _print_census(args.low, args.high, args.json)
        if args.command == "designs":
            return _print_designs(args.period, args.json)
        if args.command == "price":
            return _print_price(args.period, args.json)
    except ValueError as exc:
        print(f"BAD_REQUEST: {exc}", file=sys.stderr)
        return RC_BAD_REQUEST
    parser.print_help()
    return RC_BAD_REQUEST


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
