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
        s = solve_lattice_schedule(180, 47, k_probes=8, min_margin=Fraction(0))
        law = drift_law(s.design, s.probes, universe(180))
        self.assertEqual(law[2], DriftLaw(Fraction(1, 6), Fraction(1, 6), Fraction(0)))
        self.assertEqual(law[6], DriftLaw(Fraction(13, 18), Fraction(1, 90), Fraction(1, 5)))
        self.assertEqual(law[9], DriftLaw(Fraction(77, 90), Fraction(11, 90), Fraction(1, 45)))


if __name__ == "__main__":
    unittest.main()
