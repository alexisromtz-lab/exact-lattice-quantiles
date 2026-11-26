"""The unified law is a congruence, and the odd band reports its distance.

``block_saving(n) == (prime_excess(n) - band_slack(n)) // 2`` used to be a
census: it was checked on a population and published as measured.  It is in fact
an identity, because the two sides have the same parity for every ``n`` with at
least two distinct primes, so the floor in that division never loses a unit and
the law needs no population at all.  The second half of the file pins the other
consequence: an empty band is not a silence but a distance, carried with the
primes that produce it.
"""

from __future__ import annotations

import unittest
from itertools import combinations

from exact_lattice_quantiles.arithmetic import primes_of
from exact_lattice_quantiles.pricing import (band_slack, block_saving,
                                             excess_parity, odd_band_census,
                                             odd_members, prime_excess,
                                             sign_theorem_census)

HORIZON = 2000


def _excess(primes) -> int:
    """``sum(L / p) - L`` from the primes, without waiting for a factorisation."""
    product = 1
    for prime in primes:
        product *= prime
    return sum(product // p for p in primes) - product


def _lawful() -> list[int]:
    return [n for n in range(2, HORIZON + 1) if len(primes_of(n)) >= 2]


class TheCongruence(unittest.TestCase):
    def test_the_excess_and_the_slack_always_share_their_parity(self):
        parities = {(prime_excess(n) - band_slack(n)) % 2 for n in _lawful()}
        self.assertEqual(parities, {0})

    def test_so_the_law_needs_no_floor_to_be_exact(self):
        disagreements = [n for n in _lawful()
                         if block_saving(n) != (prime_excess(n) - band_slack(n)) // 2]
        self.assertEqual(disagreements, [])

    def test_the_parity_is_read_off_the_members_and_not_recomputed(self):
        for n in _lawful():
            self.assertEqual(excess_parity(n), prime_excess(n) % 2, n)
            self.assertEqual(excess_parity(n), (odd_members(n) - n % 2) % 2, n)

    def test_an_odd_period_walks_the_excesses_its_omega_allows(self):
        for n in _lawful():
            if n % 2:
                self.assertEqual(excess_parity(n),
                                 (len(primes_of(n)) - 1) % 2, n)

    def test_the_giuga_corner_therefore_needs_an_even_omega(self):
        # An odd Giuga number would sit at excess 1 on an odd period, and the
        # congruence leaves that residue to an even number of prime factors.
        odd_periods = [n for n in _lawful() if n % 2]
        self.assertTrue(all(excess_parity(n) == (len(primes_of(n)) - 1) % 2
                            for n in odd_periods))
        self.assertTrue(all(len(primes_of(n)) % 2 == 0
                            for n in odd_periods if excess_parity(n) == 1))

    def test_a_prime_period_has_no_parity_to_report(self):
        for n in (1, 2, 4, 8, 49, 121):
            self.assertLess(len(primes_of(n)), 2, n)
            with self.assertRaises(ValueError):
                excess_parity(n)

    def test_the_census_reports_the_congruence_as_an_identity(self):
        census = sign_theorem_census(HORIZON)
        self.assertEqual(census["law_violations"], [])
        self.assertEqual(census["law_periods"], len(_lawful()))
        self.assertIn("congruence", census["law_verdict"])


class TheMargin(unittest.TestCase):
    """An empty band with a witness is a measurement; without one it is a hunch."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.census = odd_band_census(10, 100)
        cls.nearest = cls.census["margin"]["nearest"]

    def test_the_walk_reports_the_nearest_miss(self):
        self.assertEqual(self.census["candidates"], [])
        self.assertEqual(self.nearest["omega"], 9)
        self.assertEqual(tuple(self.nearest["primes"]),
                         (3, 5, 7, 11, 13, 17, 19, 43, 47))
        self.assertEqual(self.nearest["threshold"], 8)

    def test_the_witness_is_recomputed_from_its_primes(self):
        self.assertEqual(_excess(self.nearest["primes"]),
                         self.nearest["excess"])
        self.assertEqual(self.nearest["distance"],
                         self.nearest["excess"] - self.nearest["threshold"])

    def test_the_miss_sits_above_the_band_by_more_than_its_own_slack(self):
        self.assertGreater(self.nearest["excess"], self.nearest["threshold"])
        self.assertGreater(self.nearest["excess"], 10 ** 4)

    def test_every_omega_reports_its_own_nearest_miss(self):
        rows = self.census["margin"]["by_omega"]
        self.assertEqual(sorted(rows), ["10", "9"])
        self.assertEqual(rows["9"]["excess"], self.nearest["excess"])
        for row in rows.values():
            self.assertGreater(row["distance"], 0)
            self.assertEqual(row["distance"], row["excess"] - row["threshold"])

    def test_the_floor_minimises_the_product_and_the_margin_the_excess(self):
        floor = self.census["floor"]
        self.assertEqual(floor["excess"], 108169298)
        self.assertLess(floor["product"], _product_of(self.nearest["primes"]))
        self.assertLess(self.nearest["excess"], floor["excess"])

    def test_the_parity_holds_on_every_miss_it_reports(self):
        self.assertTrue(self.census["parity"]["consistent"])
        for row in self.census["margin"]["by_omega"].values():
            self.assertEqual((row["excess"] - (row["omega"] - 1)) % 2, 0)

    def test_the_nearest_miss_agrees_with_exhaustive_enumeration(self):
        # Pruned walks are only worth their speed if they land on the same
        # optimum the unpruned product of combinations does, for every omega.
        odd = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
        census = odd_band_census(len(odd) - 3, max(odd))
        self.assertEqual(census["candidates"], [])
        for omega in (9, 10, 11):
            best = None
            for primes in combinations(odd, omega):
                total = _excess(primes)
                if total > 0 and (best is None or total < best[0]):
                    best = (total, primes)
            row = census["margin"]["by_omega"][str(omega)]
            self.assertEqual((row["excess"], tuple(row["primes"])), best)


def _product_of(primes) -> int:
    product = 1
    for prime in primes:
        product *= prime
    return product


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
