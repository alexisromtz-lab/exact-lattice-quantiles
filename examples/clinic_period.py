"""A clinic with heterogeneous service lengths is a period.

A chair that offers 15, 20, 30, 45 and 60 minute services repeats its booking
grammar every ``lcm(15, 20, 30, 45, 60) = 180`` minutes, and the candidate
cadences of that grammar are the divisors of 180 above one.  So the lattice of
``T = 180`` is the lattice of that chair, and every number below is computed
here in exact rational arithmetic -- none of them is quoted.

The interpretation is ours, the arithmetic is not: this is a *reading*, not a
scheduler.  What it reads is one specific question: a delay that hits every
service at once -- a late opening, a cleaning cycle running long -- is a
common-mode drift ``delta``, and the number worth knowing is how much margin
survives it rather than the average slot length.

Run it with ``PYTHONPATH=src python examples/clinic_period.py``.
"""

from fractions import Fraction
from math import lcm

from exact_lattice_quantiles import (grids, guaranteed_margin, instance_ceiling,
                                     ladder, minimum_cost_design,
                                     solve_lattice_schedule, universe, cell_widths)

DURATIONS = (15, 20, 30, 45, 60)


def minutes(fraction: Fraction, period: int) -> str:
    """The same quantity in the unit the clinic reads it in."""
    return f"{fraction} = {fraction * period} minute(s)"


def main() -> int:
    period = lcm(*DURATIONS)
    U = universe(period)
    ceiling = instance_ceiling(period)
    cheapest, design = minimum_cost_design(grids(period), U,
                                           cell_widths(grids(period), U),
                                           Fraction(0), period)
    at_ceiling = solve_lattice_schedule(period, max_cost=10 ** 6, k_probes=4,
                                        min_margin=ceiling.margin)
    keep, drift = Fraction(1, 30), Fraction(1, 60)
    service = solve_lattice_schedule(period, max_cost=10 ** 6, k_probes=5,
                                     min_margin=keep + drift)
    kept = guaranteed_margin(service.design, service.probes, U, drift)

    print(f"durations            {', '.join(str(d) for d in DURATIONS)} minutes")
    print(f"period               {period} minutes"
          f"  (lcm of the durations, {len(grids(period))} candidate cadences)")
    print(f"instance ceiling     {minutes(ceiling.margin, period)}"
          f"  (critical divisor {ceiling.critical_divisor})")
    print(f"cheapest design      {cheapest} points: cadences {list(design)}")
    print(f"at the ceiling       {at_ceiling.cost} points:"
          f" cadences {list(at_ceiling.design)}"
          f"  ->  holds {minutes(at_ceiling.witness_margin, period)}")
    print(f"ladder of period {period}")
    for cost, margin in sorted(ladder(period).items()):
        print(f"  {cost:5d} points -> margin {minutes(margin, period)}")
    print(f"service level        keep {minutes(keep, period)}"
          f" under a drift of {minutes(drift, period)}")
    print(f"  bought as          margin {keep + drift}, {service.cost} points,"
          f" {len(service.probes)} shared sections")
    print(f"  still keeps        {minutes(kept, period)}")
    print(f"  ceiling check      keep + drift = {keep + drift}"
          f" <= {ceiling.margin}: {keep + drift <= ceiling.margin}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
