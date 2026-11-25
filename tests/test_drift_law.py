"""The drift tolerance is linear in the margin, and then it falls off a cliff.

`drift_tolerance` buys margin with tolerance at a price of one for one, but
only up to a per-divisor edge.  These tests pin both halves: the line below the
edge, and the drop at it -- which is the part a linear extrapolation promises
and the span geometry does not pay.
"""

import unittest
from fractions import Fraction

from exact_lattice_quantiles import (DriftLaw, drift_law, drift_radius,
                                     drift_tolerance, universe)
from exact_lattice_quantiles.solve import solve_lattice_schedule

PERIODS = list(range(2, 122, 2))


def schedule(period: int, probes: int = 4):
    return solve_lattice_schedule(period, period, k_probes=probes,
                                  min_margin=Fraction(0))


class DriftLawTest(unittest.TestCase):
    def test_the_radius_is_the_tolerance_at_zero_margin(self):
        for period in PERIODS:
            s = schedule(period)
            radius = drift_radius(s.design, s.probes, universe(period))
            for d, law in drift_law(s.design, s.probes, universe(period)).items():
                self.assertEqual(law.radius, radius[d], f"period {period} divisor {d}")

    def test_the_price_is_one_for_one_below_the_critical_margin(self):
        """Halfway to the edge, the line has not bent yet."""
        for period in PERIODS:
            s = schedule(period)
            U = universe(period)
            for d, law in drift_law(s.design, s.probes, U).items():
                if law.critical_margin == 0 or law.radius == 0:
                    continue
                below = law.critical_margin / 2
                self.assertEqual(
                    drift_tolerance(s.design, s.probes, U, below)[d],
                    law.radius - below,
                    f"period {period} divisor {d} off the line at {below}",
                )

    def test_the_critical_margin_is_a_cliff_and_not_a_slope(self):
        """Past the edge the tolerance is under the line, never on it."""
        for period in PERIODS:
            s = schedule(period)
            U = universe(period)
            for d, law in drift_law(s.design, s.probes, U).items():
                if law.critical_margin == 0:
                    continue
                if law.critical_margin >= law.radius:
                    self.assertEqual(law.cliff, 0, f"period {period} divisor {d}")
                    continue
                self.assertGreater(law.cliff, 0, f"period {period} divisor {d}")
                self.assertEqual(
                    drift_tolerance(s.design, s.probes, U, law.critical_margin)[d],
                    law.radius - law.critical_margin - law.cliff,
                    f"period {period} divisor {d} cliff misreported",
                )

    def test_the_named_anchor_of_t180_at_eight_probes(self):
        # The published anchor is the pair (design 12,18,20,45 at cost 47,
        # witness margin 1/6).  The per-divisor radii below are a property of
        # the sections the solver returns, and several section sets reach 1/6
        # on this design; the radii moved up when the candidate set stopped
        # offering only band right edges, so they are pinned to the sections
        # actually returned.  The design's own radius is the minimum over
        # divisors and did not move: d = 2 breaks at 1/6.
        s = solve_lattice_schedule(180, 47, k_probes=8, min_margin=Fraction(0))
        self.assertEqual(s.design, (12, 18, 20, 45))
        self.assertEqual(s.witness_margin, Fraction(1, 6))
        U = universe(180)
        law = drift_law(s.design, s.probes, U)
        self.assertEqual(law[2], DriftLaw(Fraction(1, 6), Fraction(1, 6), Fraction(0)))
        self.assertEqual(law[6], DriftLaw(Fraction(5, 6), Fraction(1, 90), Fraction(4, 9)))
        self.assertEqual(law[9], DriftLaw(Fraction(7, 6), Fraction(1, 30), Fraction(29, 45)))
        self.assertEqual(min(v.radius for v in law.values()), Fraction(1, 6))


if __name__ == "__main__":
    unittest.main()
