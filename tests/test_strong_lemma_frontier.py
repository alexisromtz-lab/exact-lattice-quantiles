"""Where the strong reading of the lemma fails, and why the law does not care.

The lemma has two readings.  The **sign** reading -- no antichain saves where no
block of the same least common multiple saves -- is the one the law uses, and
the census has not found a counterexample.  The **ranking** reading is stronger:
the drop-one block is the antichain that saves *most*.  That one is false, and
this file pins where.

The witness is a primorial.  Divide every saving by ``period // 2`` and the
block pays ``e1 - 1`` while the pair antichain pays ``e2 - 1``, where ``e1`` is
the reciprocal sum of the prime set and ``e2`` the sum over its pairs.  Both
scale with ``period``, so the comparison is about primes alone -- and ``e2``,
growing like the square of ``e1``, overtakes it at 176 primes, the largest
being 1049.
"""

import unittest
from math import lcm, prod

from exact_lattice_quantiles.pricing import (antichain_saving_ceiling,
                                             block_saving, pair_antichain_of,
                                             pair_antichain_saving,
                                             primes_of, strong_lemma_frontier)

CROSSING = 176
LARGEST_PRIME = 1049


class TheRankingReadingIsFalse(unittest.TestCase):
    def test_the_pairs_beat_the_block_from_176_primes_on(self):
        frontier = strong_lemma_frontier(200)
        self.assertEqual(frontier["crossing_primes"], CROSSING)
        self.assertEqual(frontier["largest_prime"], LARGEST_PRIME)
        self.assertEqual(frontier["verdict"],
                         "the drop-one block is not the best-saving antichain")

    def test_at_175_primes_the_block_still_wins(self):
        below = strong_lemma_frontier(200)["below"]
        self.assertEqual(below["primes"], CROSSING - 1)
        self.assertEqual(below["largest_prime"], 1039)
        self.assertGreater(below["block_saving"], 0)
        self.assertLess(below["excess"], 0)

    def test_the_margin_of_the_defeat_is_a_fraction_of_a_percent(self):
        crossing = strong_lemma_frontier(200)["crossing"]
        self.assertGreater(crossing["excess"], 0)
        # A defeat by a hair: under one part in a thousand of the saving.
        self.assertLess(crossing["excess"] * 1000, crossing["block_saving"])
        self.assertEqual(crossing["digits"], 440)

    def test_the_integer_rule_says_the_same_as_the_reciprocals(self):
        crossing = strong_lemma_frontier(200)["crossing"]
        period = prod(self._primes_up_to(CROSSING))
        block = sum(period // p // 2 for p in primes_of(period)) - period // 2
        pairs = (sum(period // (a * b) // 2
                     for a in self._primes_up_to(CROSSING)
                     for b in self._primes_up_to(CROSSING) if a < b)
                 - period // 2)
        self.assertEqual(block, crossing["block_saving"])
        self.assertEqual(pairs, crossing["pair_saving"])
        self.assertEqual(pairs - block, crossing["excess"])

    def _primes_up_to(self, count: int) -> list[int]:
        found: list[int] = []
        candidate = 2
        while len(found) < count:
            if all(candidate % p for p in found if p * p <= candidate):
                found.append(candidate)
            candidate += 1
        return found

    def test_a_short_search_returns_no_crossing_rather_than_a_verdict(self):
        frontier = strong_lemma_frontier(100)
        self.assertIsNone(frontier["crossing_primes"])
        self.assertIn("below 100 primes", frontier["verdict"])
        self.assertIsNone(frontier["crossing"])

    def test_two_primes_cannot_be_told_apart(self):
        with self.assertRaises(ValueError):
            strong_lemma_frontier(2)


class ThePairAntichain(unittest.TestCase):
    def test_it_is_an_antichain_that_spans_the_period(self):
        for period in (180, 210, 2310, 2520):
            grids = pair_antichain_of(period)
            self.assertEqual(len(grids), len(primes_of(period)) * (len(primes_of(period)) - 1) // 2)
            self.assertTrue(all(period % d == 0 for d in grids))
            self.assertEqual(lcm(*grids), period)
            for i, left in enumerate(grids):
                for right in grids[i + 1:]:
                    self.assertNotEqual(left % right, 0, f"{period}: {right} | {left}")
                    self.assertNotEqual(right % left, 0, f"{period}: {left} | {right}")

    def test_it_is_the_products_of_two_primes_of_180(self):
        self.assertEqual(pair_antichain_of(180), (30, 18, 12))

    def test_three_primes_leave_the_family_in_the_red(self):
        self.assertEqual(pair_antichain_saving(180), -60)
        self.assertEqual(pair_antichain_saving(210), -56)

    def test_the_pairs_never_beat_the_ceiling_once_a_block_saves(self):
        for period in (60, 90, 120, 180, 420, 840):
            if block_saving(period) <= 0:
                continue
            self.assertLessEqual(pair_antichain_saving(period),
                                 antichain_saving_ceiling(period)[0],
                                 f"{period}")


if __name__ == "__main__":
    unittest.main()
