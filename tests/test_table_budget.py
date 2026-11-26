"""The table budget: a limit of the implementation, refused by name before the
allocation instead of by the kernel after it.

Two refusals with two meanings live in this package.  An odd period has no
answer at all; a period above the budget has an answer and no room to compute
it.  The tests below pin the second one to its arithmetic and to the sites that
used to allocate past it.
"""

from __future__ import annotations

import unittest
from fractions import Fraction

from exact_lattice_quantiles import (TABLE_BUDGET, MaskTableTooLarge, grids,
                                     mask_table, universe)
from exact_lattice_quantiles.coverage import (coverage_census, margin_ladder,
                                              minimum_probes,
                                              minimum_cost_design)
from exact_lattice_quantiles.lattice import cell_widths
from exact_lattice_quantiles.pricing import price_exact, price_form_census


class TestTheBudgetIsAnInequality(unittest.TestCase):
    def test_the_boundary_is_where_the_budget_says(self):
        self.assertEqual(mask_table(29, "probe"), 1 << 29)
        with self.assertRaises(MaskTableTooLarge):
            mask_table(30, "probe")

    def test_width_counts_slots_and_not_only_rows(self):
        # Cells are what the loop touches, so four rows at 2^27 is the same
        # reserve as one row at 2^29: accepted at the boundary, refused above.
        self.assertEqual(mask_table(27, "probe", width=4), 1 << 27)
        self.assertEqual(mask_table(28, "probe", width=2), 1 << 28)
        with self.assertRaises(MaskTableTooLarge):
            mask_table(28, "probe", width=4)
        with self.assertRaises(MaskTableTooLarge):
            mask_table(29, "probe", width=2)

    def test_a_width_below_one_is_a_malformed_request(self):
        # Not a budget question: a table with no rows is not a table.
        with self.assertRaises(ValueError) as caught:
            mask_table(1, "probe", width=0)
        self.assertNotIsInstance(caught.exception, MaskTableTooLarge)

    def test_the_budget_is_the_number_the_messages_quote(self):
        self.assertEqual(TABLE_BUDGET, 1 << 29)


class TestTheRefusalNamesTheInstance(unittest.TestCase):
    def test_attributes_and_rc(self):
        with self.assertRaises(MaskTableTooLarge) as caught:
            mask_table(27, "a DP", period=1440, width=36)
        exc = caught.exception
        self.assertEqual(exc.rc, 5)
        self.assertEqual(exc.period, 1440)
        self.assertEqual(exc.count, 27)
        self.assertEqual(exc.width, 36)
        self.assertEqual(exc.cells, 36 << 27)
        self.assertEqual(exc.site, "a DP")
        self.assertIn("period 1440", str(exc))
        self.assertIn(str(TABLE_BUDGET), str(exc))

    def test_a_caller_catching_value_error_keeps_working(self):
        # The package's refusals are ValueErrors first and typed verdicts
        # second, so an existing ``except ValueError`` still catches this one.
        self.assertTrue(issubclass(MaskTableTooLarge, ValueError))


class TestTheSitesRefuseBeforeAllocating(unittest.TestCase):
    def test_the_design_dp_answers_what_it_used_to_refuse(self):
        # The period the table used to name is now computed: its kernel has 18
        # divisors, so the DP asks for 2^18 masks and not for 2^36.
        gr, U = grids(1440), universe(1440)
        widths = cell_widths(gr, U)
        cost, design = minimum_cost_design(gr, U, widths, Fraction(1, 6), 1440)
        self.assertEqual(cost, 590)
        self.assertIsNotNone(design)

    def test_the_ladder_of_that_period_reaches_its_ceiling(self):
        gr, U = grids(1440), universe(1440)
        widths = cell_widths(gr, U)
        ladder = margin_ladder(gr, U, widths, 1440)
        self.assertEqual(ladder[max(ladder)], Fraction(4, 15))
        self.assertEqual(max(ladder), 735)

    def test_the_set_cover_dp_refuses_on_its_own_table(self):
        # The one site whose budget is not the calling period's design table:
        # it asks for two rows of 2^|U| and must refuse on that own count.
        with self.assertRaises(MaskTableTooLarge) as caught:
            minimum_probes(grids(2520), universe(2520), Fraction(1, 6))
        self.assertEqual(caught.exception.site, "the set-cover DP")
        self.assertEqual(caught.exception.width, 2)

    def test_the_census_names_the_period(self):
        with self.assertRaises(MaskTableTooLarge) as caught:
            coverage_census(720)
        self.assertEqual(caught.exception.period, 720)
        self.assertEqual(caught.exception.site, "the design census")

    def test_the_price_census_names_the_period(self):
        with self.assertRaises(MaskTableTooLarge) as caught:
            price_form_census(2520)
        self.assertEqual(caught.exception.period, 2520)

    def test_a_wide_kernel_is_refused_by_name(self):
        # Thirty divisors that no grid covers two at a time: nothing implies
        # anything, the kernel is the whole universe, and the table would be
        # twelve bytes wide over 2^30 states.
        n = 30
        gr = tuple(2 * (i + 1) for i in range(n))
        widths = {(M, d): Fraction(1, 2) if M == d else Fraction(0)
                  for M in gr for d in gr}
        with self.assertRaises(MaskTableTooLarge) as caught:
            minimum_cost_design(gr, gr, widths, Fraction(1, 6), 60)
        self.assertEqual(caught.exception.site, "the cheapest-design DP")
        self.assertEqual(caught.exception.period, 60)
        self.assertEqual(caught.exception.count, n)

    def test_more_grids_than_one_mask_word_still_answers(self):
        # Seventy grids over three divisors: every grid covers everything, so
        # the kernel is one divisor and the cheapest grid wins.  The grid mask
        # no longer fits in one word, which is the other branch of the same DP.
        gr = tuple(2 * (i + 1) for i in range(70))
        widths = {(M, d): Fraction(1, 2) for M in gr for d in (2, 4, 6)}
        self.assertEqual(minimum_cost_design(gr, (2, 4, 6), widths,
                                             Fraction(1, 6), 8), (1, (2,)))


class TestTheInstancesThatReturnedStillReturn(unittest.TestCase):
    def test_the_360_point_is_still_112(self):
        # The control the budget had to preserve: the largest class below the
        # boundary is computed, not refused.
        gr, U = grids(360), universe(360)
        widths = cell_widths(gr, U)
        cost, design = minimum_cost_design(gr, U, widths, Fraction(1, 6), 360)
        self.assertEqual(cost, 112)
        self.assertIsNotNone(design)

    def test_the_budget_refuses_the_table_not_the_request(self):
        # 2520 has 32 divisors and 6 maximals: the census over every request
        # is refused, and the price of the whole universe -- its frontier -- is
        # still exact and immediate.  A budget that refused both would be
        # hiding an answer it can compute.
        with self.assertRaises(MaskTableTooLarge):
            price_form_census(2520)
        self.assertEqual(price_exact(2520, universe(2520)), 649)


if __name__ == "__main__":
    unittest.main()
