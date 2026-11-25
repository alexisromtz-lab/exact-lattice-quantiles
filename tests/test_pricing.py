"""The price of a request: exact minimum, and the closed form that only works
on one of the two universes."""

from __future__ import annotations

import unittest
from fractions import Fraction

from exact_lattice_quantiles import (best_block_saving, exactness_frontier,
                                     grids, maximal_divisors, price_exact,
                                     price_form_census, price_form_sweep,
                                     price_one_per_maximal, universe, widths)


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


class TestTheSweepOverEvenPeriods(unittest.TestCase):
    def test_the_formula_is_the_price_on_every_reduced_universe_to_400(self):
        data = price_form_sweep(2, 400)
        self.assertEqual(data["periods"], 200)
        self.assertEqual(data["adequate"], 200)
        self.assertEqual(data["inadequate"], 0)
        self.assertEqual(data["skipped"], [])

    def test_several_grids_are_still_needed_almost_everywhere(self):
        # The formula is the optimum, but never because one grid would do:
        # in 192 of the 200 periods the cheapest cover beats every single grid.
        data = price_form_sweep(2, 400)
        self.assertEqual(data["cheaper_than_one_grid"], 192)
        self.assertEqual(data["costliest"], (160, 384))

    def test_the_sweep_publishes_what_it_could_not_afford(self):
        data = price_form_sweep(360, 360, reduced=False, max_masks=1 << 4)
        self.assertEqual(data["periods"], 0)
        self.assertEqual(data["skipped"], [360])


class TestTheBoundaryOfTheSweep(unittest.TestCase):
    """200/200 to 400 is a bounded census, not a law.  It breaks at 1080."""

    def test_a_block_is_an_antichain_and_the_naive_predicate_double_counts(self):
        # The 17 divisors of 180: the naive "any subset cheaper than its parts"
        # reports a deficit of 180 with the whole universe as witness.  Four of
        # the five in the naive witness divide 45, so their cost is charged
        # twice; the honest witness is {36, 60, 90}, which is what the DP finds.
        divisors_of_180 = tuple(d for d in range(2, 181) if 180 % d == 0)
        saving, witness, bound, cover = best_block_saving(divisors_of_180)
        self.assertEqual((saving, witness), (3, (36, 60, 90)))
        self.assertEqual((bound, cover), (93, 90))
        self.assertEqual(price_one_per_maximal(witness), bound)

    def test_the_reduced_universe_of_180_has_no_block_at_all(self):
        self.assertEqual(best_block_saving(universe(180)), (0, (), 0, 0))

    def test_the_first_mispriced_universe_is_1080(self):
        data = exactness_frontier()
        self.assertEqual(data["period"], 1080)
        self.assertEqual(data["witness"], (54, 90, 135))
        self.assertEqual((data["bound"], data["cover"], data["grid"]), (139, 135, 270))
        self.assertEqual((data["breaks"], data["scanned"]), (16, 4000))

    def test_the_witness_is_a_cover_the_library_itself_confirms(self):
        # The certificate is the design {270}, not a DP: at 1080 the DP is 2^23
        # masks.  The grid 270 must genuinely cover all three by its own widths.
        w = widths(1080)
        for d in (54, 90, 135):
            self.assertGreater(w[(270, d)], Fraction(0), d)

    def test_divisibility_by_30_is_necessary_and_not_sufficient(self):
        # 840 has 12 and 20 in U(840) and lacks 30, so the triple is not there
        # and no block exists; 265 of the 4000 periods hold a multiple of 30.
        U840 = universe(840)
        self.assertIn(12, U840)
        self.assertIn(20, U840)
        self.assertNotIn(30, U840)
        self.assertEqual(best_block_saving(U840)[0], 0)


if __name__ == "__main__":
    unittest.main()
