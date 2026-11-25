"""One call: the two budgets, the two exceptions, and the named shortfall."""

from __future__ import annotations

import unittest
from fractions import Fraction

from exact_lattice_quantiles import (OK, PROBES_INSUFFICIENT,
                                     InfeasibleAboveInstanceCeiling,
                                     InfeasibleWithinBudget, ladder,
                                     solve_lattice_schedule)


class TestTheThreeOutcomes(unittest.TestCase):
    def test_the_ceiling_is_reached_at_64_points_with_four_sections(self):
        s = solve_lattice_schedule(180, max_cost=64, k_probes=4,
                                   min_margin=Fraction(1, 3))
        self.assertEqual(s.design, (4, 6, 10, 12, 15, 18, 20, 45))
        self.assertEqual(s.cost, 64)
        self.assertEqual(s.witness_margin, Fraction(1, 3))
        self.assertEqual(s.probes_needed, 4)
        self.assertTrue(s.meets_requirement)

    def test_cheap_design_short_of_sections_is_named_not_silent(self):
        s = solve_lattice_schedule(180, max_cost=47, k_probes=4,
                                   min_margin=Fraction(1, 6))
        self.assertEqual(s.design_margin, Fraction(1, 6))
        self.assertEqual(s.probes_needed, 8)
        self.assertEqual(s.verdict, PROBES_INSUFFICIENT)
        self.assertFalse(s.meets_requirement)
        self.assertEqual(s.witness_margin, Fraction(1, 15))

    def test_the_same_design_with_eight_sections_is_ok(self):
        s = solve_lattice_schedule(180, max_cost=47, k_probes=8,
                                   min_margin=Fraction(1, 6))
        self.assertEqual(s.verdict, OK)
        self.assertEqual(s.witness_margin, Fraction(1, 6))

    def test_zero_margin_is_four_sections_and_not_zero(self):
        # Zero required margin used to be published as ``probes_needed == 0``
        # and ``verdict == OK``: the count came from a shortcut that never asked
        # whether a section holds its divisor, so "no section is needed" was
        # reported for a design that needs four.  Zero margin buys coverage, and
        # coverage is exactly what four sections cost.
        for k in (1, 2, 3):
            s = solve_lattice_schedule(180, max_cost=47, k_probes=k)
            self.assertEqual((s.cost, s.design), (47, (12, 18, 20, 45)))
            self.assertEqual(s.probes_needed, 4)
            self.assertEqual(s.verdict, PROBES_INSUFFICIENT)
            self.assertFalse(s.meets_requirement)
            # Nothing is held, so there is no margin to report -- which is not
            # the same statement as a margin of zero.
            self.assertIsNone(s.witness_margin)
            self.assertIsNone(s.robust_margin)
            self.assertEqual(len(s.uncovered), 12)

        s = solve_lattice_schedule(180, max_cost=47, k_probes=4)
        self.assertEqual(s.verdict, OK)
        self.assertEqual(s.probes_needed, 4)
        self.assertEqual(s.uncovered, ())
        self.assertEqual(s.witness_margin, Fraction(1, 15))

    def test_no_schedule_is_met_without_covering_its_divisors(self):
        # The dangerous direction: a schedule that does not hold some divisor
        # must never read as meeting its requirement, whatever the margin asked.
        for period in (12, 30, 60, 180):
            for k in (1, 2, 3):
                s = solve_lattice_schedule(period, max_cost=10 ** 6, k_probes=k)
                if s.uncovered:
                    self.assertFalse(s.meets_requirement, (period, k))
                    self.assertIsNone(s.witness_margin, (period, k))
                else:
                    self.assertTrue(s.meets_requirement, (period, k))


    def test_a_cover_that_is_paid_for_is_always_placed(self):
        # The belt in ``solve_lattice_schedule``: ``best_probes`` maximises the
        # *margin*, and at margin zero a subset that misses a divisor is
        # indistinguishable from one that holds it with no slack.  Where the
        # budget pays for a cover, the sections placed must therefore be a
        # cover -- measured over the cheap design of every even period up to 60,
        # at the count the design needs and one more (224 combinations up to
        # 120 agree).
        for period in range(2, 61, 2):
            need = solve_lattice_schedule(period, max_cost=10 ** 6,
                                          k_probes=64).probes_needed
            self.assertIsNotNone(need, period)
            for k in (need, need + 1):
                s = solve_lattice_schedule(period, max_cost=10 ** 6, k_probes=k)
                self.assertEqual(s.uncovered, (), (period, k))
                self.assertEqual(s.verdict, OK, (period, k))
                self.assertIsNotNone(s.witness_margin, (period, k))


class TestTwoInfeasibilities(unittest.TestCase):
    def test_above_the_ceiling_is_not_a_budget_problem(self):
        with self.assertRaises(InfeasibleAboveInstanceCeiling) as caught:
            solve_lattice_schedule(180, max_cost=10 ** 9, k_probes=32,
                                   min_margin=Fraction(2, 5))
        self.assertEqual(caught.exception.ceiling, Fraction(1, 3))
        self.assertEqual(caught.exception.rc, 1)

    def test_below_the_cheapest_design_is_a_budget_problem(self):
        with self.assertRaises(InfeasibleWithinBudget) as caught:
            solve_lattice_schedule(180, max_cost=46, k_probes=8,
                                   min_margin=Fraction(1, 6))
        self.assertEqual(caught.exception.cheapest, 47)
        self.assertEqual(caught.exception.rc, 2)

    def test_the_two_exceptions_are_not_interchangeable(self):
        self.assertFalse(issubclass(InfeasibleAboveInstanceCeiling,
                                    InfeasibleWithinBudget))
        self.assertNotEqual(InfeasibleAboveInstanceCeiling.rc,
                            InfeasibleWithinBudget.rc)

    def test_exactly_at_the_ceiling_is_solvable(self):
        # The ceiling is attained by a real design, not just bounded: for
        # T = 360 it costs 141 points and four sections.
        s = solve_lattice_schedule(360, max_cost=141, k_probes=4,
                                   min_margin=Fraction(4, 15))
        self.assertEqual(s.ceiling.margin, Fraction(4, 15))
        self.assertEqual(s.cost, 141)
        self.assertEqual(s.witness_margin, Fraction(4, 15))
        self.assertTrue(s.meets_requirement)

    def test_a_ceiling_below_one_third_still_costs_more(self):
        # 4/15 < 1/3, and yet T = 360 needs 141 points against the 64 of
        # T = 180: a lower ceiling is not a cheaper instance.
        self.assertEqual(ladder(360)[141], Fraction(4, 15))
        self.assertEqual(ladder(180)[64], Fraction(1, 3))


class TestSerialisation(unittest.TestCase):
    def test_no_float_ever_escapes(self):
        import json

        s = solve_lattice_schedule(180, max_cost=47, k_probes=4,
                                   min_margin=Fraction(1, 6))
        payload = s.as_dict()
        json.dumps(payload)  # would raise on a Fraction, which is the point
        for key in ("design_margin", "witness_margin", "robust_margin"):
            self.assertIsInstance(payload[key], str, key)
        for q in payload["probes"]:
            self.assertIsInstance(q, str)
        self.assertIsInstance(payload["ceiling"]["margin"], str)

    def test_sections_round_trip(self):
        s = solve_lattice_schedule(180, max_cost=47, k_probes=8,
                                   min_margin=Fraction(1, 6))
        self.assertEqual([Fraction(q) for q in s.as_dict()["probes"]], list(s.probes))


class TestLadder(unittest.TestCase):
    def test_ladder_is_the_same_curve_at_every_entry(self):
        curve = ladder(180)
        for cost, margin in curve.items():
            s = solve_lattice_schedule(180, max_cost=cost, k_probes=8,
                                       min_margin=margin)
            self.assertEqual(s.cost, cost, cost)
            self.assertGreaterEqual(s.design_margin, margin, cost)


if __name__ == "__main__":
    unittest.main()
