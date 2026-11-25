"""The block law: when a set of grids is cheaper bought whole.

``{12, 20, 30}`` costs ``31`` points bought separately and ``30`` as the single
grid ``60`` -- its least common multiple.  So the closed form that buys one grid
per maximal element is not the price.  This module pins the law that decides it
without a search, over the full divisor universe, and pins the ceiling of that
law: reducing a period to its carriers deletes witnesses, and the two universes
break at periods eighteen times apart.
"""

from __future__ import annotations

import itertools
import time
import unittest

from exact_lattice_quantiles import (block_catalogue, block_of, block_saving,
                                     exactness_frontier, frontier_by_blocks,
                                     price_exact, universe)
from exact_lattice_quantiles.arithmetic import divisors
from exact_lattice_quantiles.pricing import _cheapest, best_block_saving


class TestOneBlock(unittest.TestCase):
    def test_the_block_drops_one_prime(self):
        self.assertEqual(block_of(270), ((54, 90, 135), 270))
        self.assertEqual(block_of(60), ((12, 20, 30), 60))

    def test_the_block_is_an_antichain(self):
        for n in (60, 210, 420, 2310):
            block, _ = block_of(n)
            for a, b in itertools.permutations(block, 2):
                self.assertNotEqual(a % b, 0, f"{a} divides {b} in the block of {n}")

    def test_the_grid_that_covers_a_block_is_its_lcm(self):
        for n in range(2, 500):
            block, whole = block_of(n)
            if not block:
                continue
            lcm = 1
            for d in block:
                lcm = lcm * d // __import__("math").gcd(lcm, d)
            self.assertEqual(whole, lcm)

    def test_the_saving_is_a_sum_of_reciprocal_primes(self):
        self.assertEqual(block_saving(60), 1)
        self.assertEqual(block_saving(270), 4)
        self.assertEqual(block_saving(210), 18)

    def test_floors_are_not_noise(self):
        """``sum(1/p)`` above one is not enough; the points have to pay."""
        for n in (30, 858, 1722):
            self.assertGreater(sum(1 / p for p in _primes(n)), 1.0)
            self.assertEqual(block_saving(n), 0, f"{n} should not save")

    def test_a_block_can_be_cheaper_apart(self):
        self.assertEqual(block_of(6), ((2, 3), 6))
        self.assertLess(block_saving(6), 0)

    def test_primes_leave_no_block(self):
        self.assertEqual(block_of(2), ((), 1))
        self.assertEqual(block_saving(2), 0)
        self.assertEqual(block_of(4), ((2,), 2))


class TestTheLawOverEveryDivisor(unittest.TestCase):
    def test_law_matches_the_exhaustive_search(self):
        """No search is needed: ``exists N | T`` decides it, witness included."""
        for period in range(2, 1001, 2):
            divisor_universe = tuple(d for d in divisors(period) if d >= 2)
            measured, _, _, _ = best_block_saving(divisor_universe)
            predicted = max(block_saving(d) for d in divisor_universe)
            self.assertEqual(
                measured > 0, predicted > 0,
                f"period {period}: search says {measured}, law says {predicted}")
            self.assertEqual(measured, predicted, f"period {period}")

    def test_frontier_over_every_divisor_is_sixty(self):
        data = frontier_by_blocks(2000)
        self.assertEqual(data["period"], 60)
        self.assertEqual(data["witness"], (12, 20, 30))
        self.assertEqual((data["bound"], data["cover"]), (31, 30))
        self.assertEqual((data["breaks"], data["scanned"]), (82, 1000))

    def test_multiple_of_sixty_is_sufficient_and_not_necessary(self):
        """``60 | T`` always breaks; ``49`` of the ``82`` breaks to 2000 don't.

        The sufficient condition is a divisor of the period, so it can only read
        the divisors: ``90`` breaks through ``N = 90``, which is no multiple of
        ``60``.  A rule stated on multiples of ``60`` would miss six breaks in
        ten.
        """
        breaks = [T for T in range(2, 2001, 2)
                  if any(block_saving(d) > 0 for d in divisors(T) if d >= 2)]
        off = [T for T in breaks if T % 60]
        self.assertEqual(len(breaks), 82)
        self.assertEqual(off[:5], [90, 150, 210, 270, 330])
        self.assertEqual(len(off), 49)
        self.assertTrue(all(T in breaks for T in range(60, 2001, 60)))


class TestTheCatalogue(unittest.TestCase):
    def test_minimal_blocks(self):
        rows = block_catalogue(250)
        self.assertEqual([n for n, _, _ in rows], [60, 90, 150, 210])
        self.assertEqual(rows[3], (210, 18, (30, 42, 70, 105)))

    def test_every_minimal_row_really_saves(self):
        for n, saving, block in block_catalogue(2000):
            self.assertGreater(saving, 0)
            self.assertEqual(saving, block_saving(n))
            self.assertEqual(block, block_of(n)[0])
            for d in divisors(n):
                if 2 <= d < n:
                    self.assertLessEqual(
                        block_saving(d), 0, f"{n} is not minimal: {d} already saves")


class TestTwoUniverses(unittest.TestCase):
    def test_the_reduced_universe_is_eighteen_times_later(self):
        every = exactness_frontier(2000, reduced=False)
        reduced = exactness_frontier(2000, reduced=True)
        self.assertEqual((every["period"], reduced["period"]), (60, 1080))
        self.assertEqual(reduced["witness"], (54, 90, 135))

    def test_reduction_deletes_witnesses(self):
        """``840`` carries ``12`` and ``20`` but not ``30``: no block of ``60``."""
        carriers = universe(840)
        self.assertIn(12, carriers)
        self.assertIn(20, carriers)
        self.assertNotIn(30, carriers)
        self.assertEqual(exactness_frontier(2)["period"], None)


class TestTheAntichainOnly(unittest.TestCase):
    """A request costs what its maximal elements cost, so the DP can shrink."""

    def test_identity_on_the_reduced_universe(self):
        U, best, _ = _cheapest(180, True)
        for mask in range(1, 1 << len(U)):
            need = tuple(U[i] for i in range(len(U)) if mask >> i & 1)
            self.assertEqual(price_exact(180, need), best[mask])

    def test_identity_on_a_sample_of_every_divisor(self):
        U, best, _ = _cheapest(180, False)
        step = 97  # deterministic, and spread over the 2^17 masks
        for mask in range(1, 1 << len(U), step):
            need = tuple(U[i] for i in range(len(U)) if mask >> i & 1)
            self.assertEqual(price_exact(180, need, reduced=False), best[mask])

    def test_the_price_survives_where_the_full_table_died(self):
        """``2^23`` and ``2^32`` masks are not a place a table can run."""
        started = time.perf_counter()
        self.assertEqual(price_exact(1080, (54, 90, 135)), 135)
        self.assertEqual(price_exact(2520, (84, 140, 210)), 210)
        self.assertLess(time.perf_counter() - started, 1.0)

    def test_a_carrier_outside_the_universe_is_refused(self):
        with self.assertRaises(ValueError):
            price_exact(1080, (54, 90, 136))
        with self.assertRaises(ValueError):
            price_exact(180, (36, 60, 90))  # legal over every divisor, not over U
        self.assertEqual(price_exact(180, (36, 60, 90), reduced=False), 180 // 2)


def _primes(n: int) -> tuple[int, ...]:
    out, p = [], 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return tuple(out)


if __name__ == "__main__":
    unittest.main()
