"""The price of a request: exact minimum, and the closed form that only works
on one of the two universes."""

from __future__ import annotations

import unittest
from fractions import Fraction

from exact_lattice_quantiles import (grids, maximal_divisors, price_exact,
                                     price_form_census, price_one_per_maximal,
                                     universe)


class TestClosedFormIsAnUpperBound(unittest.TestCase):
    def test_the_formula_never_goes_below_the_minimum(self):
        for period in (60, 120, 180):
            U = grids(period)
            census = price_form_census(period, reduced=False)
            self.assertGreater(census.discrepancies, 0, period)
            self.assertFalse(census.adequate, period)
            self.assertLessEqual(price_exact(period, U, reduced=False),
                                 price_one_per_maximal(U))

    def test_a_prime_period_has_one_grid_so_the_formula_is_exact(self):
        # 7 has the single divisor 7 above one, so there is nothing to combine.
        self.assertEqual(price_form_census(7, reduced=False).discrepancies, 0)
        self.assertEqual(price_one_per_maximal((7,)), 3)
        self.assertEqual(price_exact(7, (7,), reduced=False), 3)

    def test_maximal_divisors_drop_the_multiples(self):
        self.assertEqual(maximal_divisors((2, 4, 6, 8)), (6, 8))
        self.assertEqual(maximal_divisors((12, 18, 20, 45)), (12, 18, 20, 45))
        self.assertEqual(maximal_divisors(()), ())


class TestTheTwoUniversesOf180(unittest.TestCase):
    """P58 measured these two rows; the package reproduces them from scratch."""

    def test_reduced_universe_the_formula_is_the_price(self):
        census = price_form_census(180)
        self.assertTrue(census.adequate)
        self.assertEqual(census.universe_size, 12)
        self.assertEqual(census.consumers, 4095)
        self.assertEqual(census.discrepancies, 0)
        self.assertEqual(census.cheapest_cover, 47)
        self.assertEqual(census.one_grid, 90)

    def test_every_divisor_the_formula_overcharges(self):
        census = price_form_census(180, reduced=False)
        self.assertFalse(census.adequate)
        self.assertEqual(census.universe_size, 17)
        self.assertEqual(census.consumers, 131071)
        self.assertEqual(census.discrepancies, 9984)
        self.assertEqual(census.share, Fraction(9984, 131071))
        self.assertEqual(census.frontiers, 173)
        self.assertEqual(census.saving_frontiers, 10)
        self.assertEqual(census.cheapest_cover, 90)
        self.assertEqual(census.one_grid, 90)

    def test_a_triple_is_bought_whole_for_less_than_its_parts(self):
        self.assertEqual(price_one_per_maximal((12, 20, 30)), 31)
        self.assertEqual(price_exact(180, (12, 20, 30), reduced=False), 30)
        self.assertEqual(price_one_per_maximal((36, 60, 90)), 93)
        self.assertEqual(price_exact(180, (36, 60, 90), reduced=False), 90)

    def test_a_pair_never_saves(self):
        # lcm(a, b) >= a + b for incomparable a, b, so no pair can be cheaper.
        U = grids(180)
        for i, a in enumerate(U):
            for b in U[i + 1:]:
                if a % b == 0 or b % a == 0:
                    continue
                self.assertEqual(price_one_per_maximal((a, b)),
                                 price_exact(180, (a, b), reduced=False), (a, b))

    def test_the_formula_never_crosses_the_ceiling_of_one_grid(self):
        census = price_form_census(180, reduced=False)
        self.assertEqual(census.cheapest_cover, census.one_grid)
        self.assertLessEqual(price_exact(180, universe(180)), census.one_grid)

    def test_the_dict_form_names_the_universe(self):
        self.assertEqual(price_form_census(180).as_dict()["universe"], "reduced")
        self.assertEqual(price_form_census(180, reduced=False).as_dict()["universe"],
                         "divisors")


if __name__ == "__main__":
    unittest.main()
