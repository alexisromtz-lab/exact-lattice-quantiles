"""A design has no cliff: the line runs all the way to zero tolerance.

`drift_law` reports, per divisor, a critical margin followed by a drop.  Those
drops never reach a design, because taking the minimum over divisors favours
the divisor that holds the radius, and that divisor's own regime is the longest
one.  These tests pin the contrast: the divisor breaks, the design does not.
"""

import unittest
from fractions import Fraction

from exact_lattice_quantiles import (design_law, drift_law, drift_radius,
                                     drift_tolerance, universe)
from exact_lattice_quantiles.coverage import _reach_from_zero, _shrinks, _spans
from exact_lattice_quantiles.solve import (InfeasibleWithinBudget,
                                           solve_lattice_schedule)


def design_tolerance(G, Q, U, margin):
    return min(drift_tolerance(G, Q, U, margin).values())


def assert_linear_over_the_range(case, G, Q, U):
    """Check the identity on every interval between consecutive knots.

    The tolerance is a minimum of piecewise-linear functions of slope -1 with
    downward jumps at the knots, so one interior point per interval settles the
    whole interval.
    """
    radius = design_law(G, Q, U).radius
    knots = sorted({2 * knot for d in U for knot in _shrinks(_spans(G, Q, d))
                    if 0 < 2 * knot < radius})
    bounds = [Fraction(0)] + knots + [radius]
    for low, high in zip(bounds, bounds[1:]):
        middle = (low + high) / 2
        case.assertEqual(design_tolerance(G, Q, U, middle), radius - middle,
                         f"margin {middle} of {radius}")


class DesignLinearityTest(unittest.TestCase):
    def test_the_cheap_anchor_is_linear_to_zero(self):
        s = solve_lattice_schedule(180, max_cost=64, k_probes=4,
                                   min_margin=Fraction(0))
        U = universe(180)
        self.assertEqual(design_law(s.design, s.probes, U).radius, Fraction(1, 15))
        assert_linear_over_the_range(self, s.design, s.probes, U)

    def test_the_wide_anchor_is_linear_to_zero(self):
        s = solve_lattice_schedule(180, max_cost=64, k_probes=4,
                                   min_margin=Fraction(1, 3))
        U = universe(180)
        shape = design_law(s.design, s.probes, U)
        self.assertEqual(shape.radius, Fraction(1, 3))
        assert_linear_over_the_range(self, s.design, s.probes, U)

    def test_a_divisor_break_does_not_reach_the_design(self):
        s = solve_lattice_schedule(20, max_cost=64, k_probes=4,
                                   min_margin=Fraction(0))
        U = universe(20)
        G, Q = s.design, s.probes
        shape = design_law(G, Q, U)
        radii = drift_radius(G, Q, U)
        laws = drift_law(G, Q, U)
        self.assertEqual(shape.radius, min(radii.values()))
        self.assertLess(shape.first_divisor_break, shape.radius)
        self.assertTrue(any(law.cliff > 0 for d, law in laws.items()
                            if radii[d] > shape.radius))
        past = (shape.first_divisor_break + shape.radius) / 2
        self.assertEqual(design_tolerance(G, Q, U, past), shape.radius - past)

    def test_first_divisor_break_is_the_narrowest_divisor(self):
        for period in range(2, 122, 2):
            s = solve_lattice_schedule(period, period, k_probes=4,
                                       min_margin=Fraction(0))
            U = universe(period)
            shape = design_law(s.design, s.probes, U)
            laws = drift_law(s.design, s.probes, U)
            self.assertEqual(shape.first_divisor_break,
                             min(law.critical_margin for law in laws.values()),
                             f"period {period}")

    def test_linearity_survives_a_sweep_of_budgets_and_probe_counts(self):
        measured = 0
        for period in (6, 20, 30, 60, 180):
            for cost in (16, 47, 64, 96):
                for probes in (2, 4, 6, 8):
                    try:
                        s = solve_lattice_schedule(period, max_cost=cost,
                                                   k_probes=probes,
                                                   min_margin=Fraction(0))
                    except InfeasibleWithinBudget:
                        continue
                    if not s.probes:
                        continue
                    U = universe(period)
                    with self.subTest(period=period, cost=cost, probes=probes):
                        assert_linear_over_the_range(self, s.design, s.probes, U)
                    measured += 1
        self.assertGreater(measured, 40)

    def test_the_undoubled_view_is_the_same_measurement(self):
        s = solve_lattice_schedule(180, max_cost=64, k_probes=4,
                                   min_margin=Fraction(1, 3))
        U = universe(180)
        G, Q = s.design, s.probes
        half = design_law(G, Q, U, doubled=False)
        doubled = design_law(G, Q, U)
        self.assertEqual(half.radius * 2, doubled.radius)
        self.assertEqual(half.first_divisor_break * 2, doubled.first_divisor_break)
        for d, law in drift_law(G, Q, U, doubled=False).items():
            self.assertEqual(law.critical_margin * 2,
                             drift_law(G, Q, U)[d].critical_margin, f"divisor {d}")
            self.assertEqual(
                drift_tolerance(G, Q, U, law.critical_margin, doubled=False)[d],
                law.radius - law.critical_margin - law.cliff, f"divisor {d}")


class DesignReachSanityTest(unittest.TestCase):
    """The knot set must contain every place the reach can bend or jump."""

    def test_the_reach_is_zero_at_the_radius(self):
        s = solve_lattice_schedule(180, max_cost=64, k_probes=4,
                                   min_margin=Fraction(1, 3))
        U = universe(180)
        h = design_law(s.design, s.probes, U).radius / 2
        for d in U:
            self.assertEqual(
                _reach_from_zero([(a + h, b - h) for a, b in _spans(s.design, s.probes, d)]),
                Fraction(0), f"divisor {d}")


if __name__ == "__main__":
    unittest.main()
