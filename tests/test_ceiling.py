"""The ceiling: three routes to the same number, and the census that keeps
``1/3`` from passing itself off as a law."""

from __future__ import annotations

import unittest
from fractions import Fraction

from exact_lattice_quantiles import (OddPeriodUnsupported,
                                     ceiling_by_product, ceiling_by_scan,
                                     ceiling_census, instance_ceiling, mass,
                                     prime_product, universe)
from exact_lattice_quantiles.ceiling import InfeasibleAboveInstanceCeiling
from exact_lattice_quantiles.solve import solve_lattice_schedule


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


def _phi(d: int) -> int:
    return sum(1 for i in range(1, d + 1) if _gcd(i, d) == 1)


class TestMass(unittest.TestCase):
    def test_two_is_one_class_not_half_of_one(self):
        # phi(2) / 2 = 1/2 is not a number of classes: the only i with
        # 1 <= i <= 1 and gcd(i, 2) == 1 is i = 1.
        self.assertEqual(mass(2), 2)

    def test_mass_is_totient_above_two(self):
        for d in (3, 4, 5, 6, 9, 10, 12, 15, 18, 20, 45, 90, 210):
            self.assertEqual(mass(d), _phi(d), d)


class TestThreeRoutes(unittest.TestCase):
    def test_routes_agree_on_every_period_to_400(self):
        for T in range(4, 401, 2):
            closed = instance_ceiling(T).margin
            self.assertEqual(closed, ceiling_by_product(T), T)
            self.assertEqual(closed, ceiling_by_scan(T), T)

    def test_routes_agree_out_of_sample(self):
        for T in (720, 840, 1260, 2520):
            closed = instance_ceiling(T).margin
            self.assertEqual(closed, ceiling_by_product(T), T)
            self.assertEqual(closed, ceiling_by_scan(T), T)

    def test_out_of_sample_values(self):
        self.assertEqual(instance_ceiling(720).margin, Fraction(4, 15))
        self.assertEqual(instance_ceiling(720).critical_divisor, 90)
        for T in (840, 1260, 2520):
            self.assertEqual(instance_ceiling(T).margin, Fraction(8, 35), T)
            self.assertEqual(instance_ceiling(T).critical_divisor, 210, T)

    def test_two_is_the_one_degenerate_period(self):
        # The product route divides by a mass ratio, and at d = 2 there is no
        # ratio: mass(2) = 2 while phi(2) = 1.  T = 2 is the only period where
        # the two routes disagree, and the disagreement is the measurement that
        # locates the convention.
        self.assertEqual(instance_ceiling(2).margin, Fraction(1))
        self.assertEqual(ceiling_by_product(2), Fraction(1))
        self.assertEqual(prime_product(2), Fraction(1, 2))
        disagreeing = [T for T in range(2, 401, 2)
                       if ceiling_by_product(T) != instance_ceiling(T).margin]
        self.assertEqual(disagreeing, [])


class TestCeilingOf180(unittest.TestCase):
    def test_ceiling_is_one_third_tied_by_six(self):
        c = instance_ceiling(180)
        self.assertEqual(c.margin, Fraction(1, 3))
        self.assertEqual(c.critical_divisor, 6)
        self.assertEqual(c.universe,
                         (2, 3, 4, 5, 6, 9, 10, 12, 15, 18, 20, 45))

    def test_one_third_is_the_product_route_not_the_bare_prime_product(self):
        # The divisor that ties the ceiling is 6, which carries no factor 5, so
        # F = {5} and 1/3 = P(180) * 5/4 = 4/15 * 5/4.
        c = instance_ceiling(180)
        self.assertEqual(c.prime_product, Fraction(4, 15))
        self.assertEqual(c.dropped, (5,))
        self.assertFalse(c.radical_bound_hit)
        self.assertEqual(c.prime_product * Fraction(5, 4), c.margin)

    def test_percent_is_exact(self):
        self.assertEqual(instance_ceiling(180).percent, Fraction(100, 3))


class TestCensus(unittest.TestCase):
    def test_167_above_30_equal_3_below(self):
        data = ceiling_census(2, 400)
        self.assertEqual(data["periods"], 200)
        self.assertEqual(data["above"], 167)
        self.assertEqual(data["equal"], 30)
        self.assertEqual(data["below"], 3)

    def test_the_three_that_fall_below(self):
        third = Fraction(1, 3)
        below = {T: instance_ceiling(T).margin for T in range(2, 401, 2)
                 if instance_ceiling(T).margin < third}
        self.assertEqual(below, {252: Fraction(2, 7), 360: Fraction(4, 15),
                                 396: Fraction(10, 33)})

    def test_one_third_is_not_a_universal_bound(self):
        self.assertGreater(instance_ceiling(90).margin, Fraction(1, 3))
        self.assertEqual(ceiling_census(2, 400)["lowest"], (Fraction(4, 15), 360))


class TestCeilingIsPerInstance(unittest.TestCase):
    def test_same_radical_three_ceilings(self):
        # rad(T) = 30 for all six; the exponents decide, through U(T).
        got = {T: instance_ceiling(T).margin for T in (60, 90, 120, 180, 240, 360)}
        self.assertEqual(got, {
            60: Fraction(1, 3),
            90: Fraction(8, 15),
            120: Fraction(1, 3),
            180: Fraction(1, 3),
            240: Fraction(1, 3),
            360: Fraction(4, 15),
        })

    def test_powers_of_two_do_not_collapse(self):
        for T in (4, 8, 16, 32, 64, 128, 256):
            self.assertEqual(instance_ceiling(T).margin, Fraction(1, 2), T)
            self.assertEqual(prime_product(T), Fraction(1, 2), T)


class TestDroppedPrimes(unittest.TestCase):
    def test_census_of_dropped_sets_to_400(self):
        counts = {0: 0, 1: 0}
        for T in range(2, 401, 2):
            n = len(instance_ceiling(T).dropped)
            self.assertLessEqual(n, 1, T)
            counts[n] = counts.get(n, 0) + 1
        self.assertEqual(counts, {0: 84, 1: 116})

    def test_two_dropped_primes_are_not_excluded(self):
        # |F| <= 1 is a regularity of one range, not a law of the problem.
        two = [T for T in range(4, 2001, 2) if len(instance_ceiling(T).dropped) == 2]
        self.assertEqual(two, [420, 660, 780, 924, 1020, 1092, 1140, 1380,
                               1428, 1596, 1716, 1740, 1848, 1860, 1932])
        for T in two:
            self.assertEqual(instance_ceiling(T).margin, Fraction(1, 3), T)


class TestTheDomain(unittest.TestCase):
    def test_the_universe_is_empty_exactly_when_the_period_is_odd(self):
        # Every divisor of an odd T is odd, and only an even grid carries a
        # divisor: so U(T) is empty for every odd T and non-empty for every
        # even one.  400 periods, no search involved.
        for T in range(1, 401):
            self.assertEqual(bool(universe(T)), T % 2 == 0, T)

    def test_an_odd_period_is_refused_by_name(self):
        # ``T = 1`` is refused earlier and for another reason: it has no
        # phase of two classes at all.
        for T in (5, 7, 9, 11, 15, 21, 45, 77, 105):
            with self.assertRaises(OddPeriodUnsupported, msg=T) as caught:
                instance_ceiling(T)
            self.assertEqual(caught.exception.period, T)
            self.assertIn(str(T), str(caught.exception))
            self.assertIn("even", str(caught.exception))

    def test_the_refusal_is_still_a_value_error(self):
        # Every route raised ``ValueError`` before -- except that it raised
        # CPython's "min() iterable argument is empty", which names the
        # arithmetic.  A caller that caught the old one still catches this.
        with self.assertRaises(ValueError):
            instance_ceiling(9)

    def test_every_route_refuses_the_same_way(self):
        for route in (instance_ceiling, ceiling_by_product, ceiling_by_scan):
            with self.assertRaises(OddPeriodUnsupported, msg=route.__name__):
                route(15)
        with self.assertRaises(OddPeriodUnsupported):
            solve_lattice_schedule(15, max_cost=64, k_probes=4, min_margin=0)

    def test_the_even_neighbour_is_untouched(self):
        self.assertEqual(instance_ceiling(180).margin, Fraction(1, 3))
        self.assertEqual(universe(180)[-1], 45)


class TestInfeasible(unittest.TestCase):
    def test_the_exception_carries_the_ceiling(self):
        with self.assertRaises(InfeasibleAboveInstanceCeiling) as caught:
            raise_if_above(180, Fraction(2, 5))
        self.assertEqual(caught.exception.ceiling, Fraction(1, 3))
        self.assertEqual(caught.exception.needed, Fraction(2, 5))
        self.assertEqual(caught.exception.rc, 1)

    def test_exactly_at_the_ceiling_is_not_infeasible(self):
        raise_if_above(180, Fraction(1, 3))


def raise_if_above(period: int, margin: Fraction) -> None:
    c = instance_ceiling(period)
    if margin > c.margin:
        raise InfeasibleAboveInstanceCeiling(period, margin, c.margin)


if __name__ == "__main__":
    unittest.main()
