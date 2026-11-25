"""One law, two universes: the block law with the universe as a filter.

The law over every divisor says a period is mispriced exactly when one of its
divisors has a drop-one block that saves.  That sentence survives the move to
the reduced universe ``U(T)`` only after one change: the block is *kept* where
the universe carries it, and the cover is recomputed over what survives.

``U(2520)`` holds ``84``, ``140`` and ``210`` and not ``60``, though all four
are the block of ``420``.  The kept block costs ``84 + 140 + 210 = 434`` apart
against ``420`` as one grid, so the saving is ``7`` points and ``2520`` breaks.
``U(840)`` holds ``12`` and ``20`` and not ``30``: the block of ``60`` survives
there with two elements, and two incomparable grids never save -- which is why
``840`` is a multiple of ``60`` and still does not break.

This module also pins the limit of the law: **not every saving set is a block.**
In ``D(1050)`` the antichain ``{42, 150, 350, 525}`` saves ``8`` points and its
reciprocals are ``1/2 + 1/3 + 1/7 + 1/25 > 1``, with ``42`` no ``1050 / p``.
What survives is the maximum: the block of ``1050`` still takes ``92``.
"""

from __future__ import annotations

import itertools
import unittest

from exact_lattice_quantiles import (best_block_saving, block_of, block_saving,
                                     carrier_block_of, carrier_block_saving,
                                     exactness_frontier, frontier_by_blocks,
                                     universe)
from exact_lattice_quantiles.arithmetic import divisors


def law_saving(period: int, reduced: bool) -> int:
    """The law's answer for one period: best kept block over its divisors."""
    carried = set(universe(period)) if reduced else set(divisors(period))
    return max((carrier_block_saving(n, carried)
                for n in divisors(period) if n >= 2), default=0)


class TestTheFilterIsTheLaw(unittest.TestCase):
    def test_the_block_of_420_loses_exactly_60_in_2520(self):
        self.assertEqual(block_of(420), ((60, 84, 140, 210), 420))
        self.assertEqual(carrier_block_of(420, set(universe(2520))),
                         ((84, 140, 210), 420))

    def test_the_kept_block_saves_where_the_whole_one_cannot_be_carried(self):
        kept = carrier_block_of(420, set(universe(2520)))
        self.assertEqual(sum(d // 2 for d in kept[0]) - kept[1] // 2, 7)
        reported = exactness_frontier(8000, reduced=True)
        self.assertEqual(reported["period"], 1080)

    def test_840_is_a_multiple_of_60_that_does_not_break(self):
        carriers = set(universe(840))
        self.assertEqual(carrier_block_of(60, carriers), ((12, 20), 60))
        self.assertLess(carrier_block_saving(60, carriers), 0)
        self.assertEqual(carrier_block_saving(420, carriers), -35)

    def test_a_kept_pair_never_saves(self):
        """Two incomparable grids cost at least as much whole as apart."""
        for n in range(2, 1001):
            block, _ = block_of(n)
            for pair in itertools.combinations(block, 2):
                self.assertLessEqual(carrier_block_saving(n, set(pair)), 0,
                                     f"{pair} of the block of {n}")

    def test_a_prime_keeps_nothing(self):
        self.assertEqual(carrier_block_of(2, set(universe(180))), ((), 1))
        self.assertEqual(carrier_block_saving(2, set(universe(180))), 0)


class TestTheLawOverBothUniverses(unittest.TestCase):
    def test_law_matches_the_exhaustive_search(self):
        """The law names the value, not only the fact, in both universes."""
        for period in range(2, 2001, 2):
            for reduced in (False, True):
                carried = (universe(period) if reduced
                           else tuple(d for d in divisors(period) if d >= 2))
                measured = best_block_saving(carried)[0]
                self.assertEqual(
                    measured, law_saving(period, reduced),
                    f"period {period}, reduced={reduced}")

    def test_both_frontiers_and_their_counts(self):
        every = frontier_by_blocks(8000)
        reduced = frontier_by_blocks(8000, reduced=True)
        self.assertEqual((every["period"], every["breaks"]), (60, 336))
        self.assertEqual((reduced["period"], reduced["breaks"]), (1080, 16))
        self.assertEqual((reduced["witness"], reduced["grid"]), ((54, 90, 135), 270))
        self.assertEqual((every["bound"], every["cover"]), (31, 30))
        self.assertEqual(reduced["law"], "carrier blocks")

    def test_reduced_breaks_are_filtered_blocks_of_a_divisor(self):
        """``2520`` breaks through ``420``, which is not even in ``U(2520)``."""
        carriers = set(universe(2520))
        self.assertNotIn(420, carriers)
        self.assertEqual(carrier_block_of(420, carriers)[0], (84, 140, 210))
        rows = [T for T in range(2, 8001, 2) if law_saving(T, True) > 0]
        self.assertEqual(rows[:5], [1080, 2160, 2520, 3240, 3960])


class TestTheMaximumIsABlock(unittest.TestCase):
    def test_not_every_saver_is_a_block(self):
        """A saving antichain that is no ``lcm / prime`` -- and still not the max."""
        witness = (42, 150, 350, 525)
        saving = sum(d // 2 for d in witness) - 1050 // 2
        self.assertEqual(saving, 8)
        self.assertNotIn(42, block_of(1050)[0])
        self.assertEqual(block_saving(1050), 92)

    def test_the_search_maximum_is_a_block_on_every_break(self):
        for period in range(2, 1501, 2):
            divisor_universe = tuple(d for d in divisors(period) if d >= 2)
            saving, witness, _, _ = best_block_saving(divisor_universe)
            if saving <= 0:
                continue
            blocks = {block_of(n)[0] for n in divisor_universe}
            self.assertIn(tuple(witness), blocks, f"period {period}")


if __name__ == "__main__":
    unittest.main()
