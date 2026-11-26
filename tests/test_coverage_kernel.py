"""The kernel of the coverage matrix: the cheapest-design DP only has to tell
apart the divisors whose fate some other divisor's fate does not already fix.

The claim is an equivalence, so it is pinned against an enumeration of every
design rather than against a second DP: ``2^|grids|`` is small enough on
``T = 180`` and knows nothing about how the DP is written.
"""

from __future__ import annotations

import unittest
from fractions import Fraction

from exact_lattice_quantiles import (design_margin, grids, instance_ceiling,
                                     universe)
from exact_lattice_quantiles.coverage import (coverage_kernel,
                                              minimum_cost_design)
from exact_lattice_quantiles.lattice import cell_widths


def _widths(period: int) -> dict[tuple[int, int], Fraction]:
    return cell_widths(grids(period), universe(period))


def _brute_force(period: int, threshold: Fraction) -> int | None:
    """Cheapest cost over every design, by enumeration instead of by DP."""
    gr, U = grids(period), universe(period)
    widths = _widths(period)
    reachable = [sum(1 << i for i, d in enumerate(U)
                     if widths.get((M, d), Fraction(0)) >= threshold)
                 for M in gr]
    full = (1 << len(U)) - 1
    best = None
    covers = [0] * (1 << len(gr))
    costs = [0] * (1 << len(gr))
    for bits in range(1, 1 << len(gr)):
        low = bits & -bits
        i = low.bit_length() - 1
        rest = bits ^ low
        covers[bits] = covers[rest] | reachable[i]
        costs[bits] = costs[rest] + gr[i] // 2
        if covers[bits] == full and (best is None or costs[bits] < best):
            best = costs[bits]
    return best


class TestTheKernelIsExact(unittest.TestCase):
    def test_every_level_of_180_agrees_with_the_enumeration(self):
        gr, U = grids(180), universe(180)
        widths = _widths(180)
        for v in sorted({v for v in widths.values()}):
            self.assertEqual(minimum_cost_design(gr, U, widths, v, 180)[0],
                             _brute_force(180, v), v)

    def test_a_level_above_the_ceiling_answers_none_on_both_sides(self):
        gr, U = grids(180), universe(180)
        widths = _widths(180)
        above = instance_ceiling(180).margin + Fraction(1, 100)
        self.assertEqual(minimum_cost_design(gr, U, widths, above, 180),
                         (None, None))
        self.assertIsNone(_brute_force(180, above))


class TestTheKernelIsStrictlySmaller(unittest.TestCase):
    def test_the_measured_sizes_and_their_costs(self):
        # The reference points of the package, and the periods the old table
        # refused: 2520 at 1/6 would have asked for 2^32 masks.
        measured = {180: (12, 5, 47), 360: (17, 9, 112), 720: (22, 14, 271),
                    1260: (23, 12, 431), 1680: (27, 19, 652), 2520: (32, 20, 888)}
        for period, (divisors, kernel, cost) in measured.items():
            gr, U = grids(period), universe(period)
            widths = _widths(period)
            kept, _ = coverage_kernel(gr, U, widths, Fraction(1, 6))
            self.assertEqual((len(U), len(kept)), (divisors, kernel), period)
            self.assertEqual(minimum_cost_design(gr, U, widths, Fraction(1, 6),
                                                 period)[0], cost, period)

    def test_the_kernel_is_a_property_of_the_threshold(self):
        # It is not a constant of the period: near the ceiling almost nothing
        # is implied any more, which is why the guard is on the kernel and not
        # on an instance.
        gr, U = grids(2520), universe(2520)
        widths = _widths(2520)
        low, _ = coverage_kernel(gr, U, widths, Fraction(1, 126))
        high, _ = coverage_kernel(gr, U, widths, instance_ceiling(2520).margin)
        self.assertLess(len(low), len(high))
        self.assertEqual(len(high), 23)


class TestADroppedDivisorIsImplied(unittest.TestCase):
    def test_every_dropped_divisor_is_contained_in_a_kept_one(self):
        for period, v in ((2520, Fraction(1, 6)), (1680, Fraction(1, 7)),
                          (360, Fraction(1, 3)), (180, Fraction(2, 9))):
            gr, U = grids(period), universe(period)
            widths = _widths(period)
            kept, columns = coverage_kernel(gr, U, widths, v)
            for d in U:
                if d in kept:
                    continue
                self.assertTrue(
                    any(columns[e] & ~columns[d] == 0 for e in kept),
                    (period, v, d))

    def test_the_design_that_covers_the_kernel_covers_the_universe(self):
        for period, v in ((2520, Fraction(1, 5)), (1680, Fraction(1, 6)),
                          (720, Fraction(1, 4))):
            gr, U = grids(period), universe(period)
            widths = _widths(period)
            cost, design = minimum_cost_design(gr, U, widths, v, period)
            self.assertIsNotNone(design)
            self.assertGreaterEqual(design_margin(design, U, widths), v)
            self.assertEqual(sum(M // 2 for M in design), cost)

    def test_margin_zero_is_coverage_and_not_every_pair(self):
        # At threshold 0 a pair that does not exist would still read as "wide
        # enough" if the columns were built with a default, and every grid
        # would appear to cover every divisor: the cheapest design would shrink
        # to a single grid.  Margin zero is the coverage price, 47 for 180.
        gr, U = grids(180), universe(180)
        widths = _widths(180)
        cost, design = minimum_cost_design(gr, U, widths, Fraction(0), 180)
        self.assertEqual((cost, design), (47, (12, 18, 20, 45)))
        self.assertEqual(set(design) & set(U), set(design))
        for d in U:
            self.assertTrue(any((M, d) in widths for M in design), d)

    def test_above_the_ceiling_the_kernel_is_one_empty_column(self):
        gr, U = grids(2520), universe(2520)
        widths = _widths(2520)
        above = instance_ceiling(2520).margin + Fraction(1, 100)
        kept, columns = coverage_kernel(gr, U, widths, above)
        self.assertEqual(len(kept), 1)
        self.assertEqual(columns[kept[0]], 0)


if __name__ == "__main__":
    unittest.main()
