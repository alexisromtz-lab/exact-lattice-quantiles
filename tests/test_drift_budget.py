"""Drift is not a second budget: it is spent from the margin.

`drift_tolerance` reads a margin and answers the drift it tolerates.
`guaranteed_margin` reads the drift and answers the margin that survives it.
These tests pin the exchange between the two, and pin the one place it does not
hold: a design whose coverage radius sits above its witness margin carries
reach it cannot spend, because its tolerance breaks at the witness margin
instead of paying the line all the way down.

The supremum is exact and is *approached*, not held, at a break: every margin
strictly below the answer tolerates the drift, the answer itself may not.
"""

import unittest
from fractions import Fraction as Fraction

from exact_lattice_quantiles import (cell_widths, cover_masks, design_law,
                                     drift_radius, drift_tolerance, grids,
                                     guaranteed_margin, margin_under_drift,
                                     minimum_cost_design, solve_lattice_schedule,
                                     universe, witness_margin)
from exact_lattice_quantiles.coverage import all_designs


def design_tolerance(G, Q, U, margin):
    return min(drift_tolerance(G, Q, U, margin).values())


class OneBudgetTest(unittest.TestCase):
    def test_no_drift_keeps_the_whole_radius(self):
        s = solve_lattice_schedule(180, max_cost=64, k_probes=4, min_margin=Fraction(0))
        U = universe(180)
        self.assertEqual(guaranteed_margin(s.design, s.probes, U),
                         design_law(s.design, s.probes, U).radius)

    def test_a_cheapest_design_spends_its_margin_one_for_one(self):
        U = universe(180)
        for k in (2, 4):
            for delta in (Fraction(1, 120), Fraction(1, 60), Fraction(1, 30)):
                s = solve_lattice_schedule(180, max_cost=64, k_probes=k, min_margin=delta)
                if not s.probes or s.witness_margin is None:
                    continue
                kept = guaranteed_margin(s.design, s.probes, U, delta)
                self.assertEqual(kept, s.witness_margin - delta,
                                 f"k={k} delta={delta}")
                self.assertEqual(kept, design_law(s.design, s.probes, U).radius - delta,
                                 f"k={k} delta={delta}: the radius is the witness margin")

    def test_the_answer_round_trips_to_the_drift(self):
        U = universe(180)
        for delta in (Fraction(1, 60), Fraction(1, 30), Fraction(1, 20)):
            s = solve_lattice_schedule(180, max_cost=64, k_probes=4, min_margin=delta)
            kept = guaranteed_margin(s.design, s.probes, U, delta)
            self.assertEqual(design_tolerance(s.design, s.probes, U, kept), delta)

    def test_the_guarantee_is_strict_at_the_answer(self):
        U = universe(180)
        delta = Fraction(1, 30)
        s = solve_lattice_schedule(180, max_cost=64, k_probes=4, min_margin=delta)
        kept = guaranteed_margin(s.design, s.probes, U, delta)
        self.assertGreaterEqual(design_tolerance(s.design, s.probes, U, kept - Fraction(1, 500)),
                                delta)
        self.assertLess(design_tolerance(s.design, s.probes, U, kept + Fraction(1, 500)), delta)

    def test_the_drift_is_not_a_second_budget(self):
        """The witness margin is one budget, split between margin and drift."""
        U = universe(180)
        delta = Fraction(1, 60)
        s = solve_lattice_schedule(180, max_cost=64, k_probes=4,
                                   min_margin=Fraction(1, 30))
        self.assertEqual(s.witness_margin, Fraction(1, 15))
        self.assertEqual(guaranteed_margin(s.design, s.probes, U, delta), Fraction(1, 20))
        self.assertEqual(guaranteed_margin(s.design, s.probes, U, delta) + delta,
                         s.witness_margin)

    def test_a_census_of_cheapest_designs_keeps_the_law(self):
        for period in (12, 20, 30, 60):
            U = universe(period)
            gr = grids(period)
            cost, _ = minimum_cost_design(gr, U, cell_widths(gr, U), Fraction(0))
            for k in (2, 4):
                for delta in (Fraction(0), Fraction(1, 60), Fraction(1, 30)):
                    s = solve_lattice_schedule(period, max_cost=cost, k_probes=k,
                                               min_margin=delta)
                    if not s.probes:
                        continue
                    self.assertEqual(guaranteed_margin(s.design, s.probes, U, delta),
                                     max(Fraction(0), s.witness_margin - delta),
                                     f"period={period} k={k} delta={delta}")

    def test_per_divisor_readings_agree_with_the_design(self):
        U = universe(180)
        delta = Fraction(1, 60)
        s = solve_lattice_schedule(180, max_cost=64, k_probes=4, min_margin=delta)
        per_divisor = margin_under_drift(s.design, s.probes, U, delta)
        self.assertEqual(min(per_divisor.values()),
                         guaranteed_margin(s.design, s.probes, U, delta))

    def test_a_negative_drift_is_refused(self):
        U = universe(12)
        with self.assertRaises(ValueError):
            margin_under_drift((4, 6), (Fraction(1, 6),), U, Fraction(-1, 4))


class ReachIsNotSpendableTest(unittest.TestCase):
    """The gap between coverage and margin cannot be traded for drift.

    ``(4, 6)`` at ``T = 12`` covers a radius of ``1/3`` while its witness margin
    is ``1/6``: the extra reach exists only because a drift may push a section
    into a block it never captured.  Past the witness margin the tolerance is
    zero, so the gap buys nothing.
    """

    U = universe(12)
    G = (4, 6)
    Q = (Fraction(1, 6), Fraction(1, 4), Fraction(1, 2), Fraction(7, 12), Fraction(11, 12))

    def test_the_radius_exceeds_the_witness_margin(self):
        self.assertEqual(min(drift_radius(self.G, self.Q, self.U).values()), Fraction(1, 3))
        self.assertEqual(min(witness_margin(self.G, self.Q, self.U).values()), Fraction(1, 6))

    def test_the_line_pays_until_the_witness_margin(self):
        self.assertEqual(design_tolerance(self.G, self.Q, self.U, Fraction(1, 12)),
                         Fraction(1, 4))

    def test_and_stops_there(self):
        self.assertEqual(design_tolerance(self.G, self.Q, self.U, Fraction(1, 6)),
                         Fraction(0))

    def test_so_the_gap_is_not_a_budget(self):
        """The exchange rate is the witness margin, never the radius."""
        for delta in (Fraction(1, 100), Fraction(1, 12), Fraction(1, 8)):
            self.assertEqual(guaranteed_margin(self.G, self.Q, self.U, delta),
                             Fraction(1, 6), f"delta={delta}")


class TheRadiusNeverFallsBelowTheMarginTest(unittest.TestCase):
    """``radius >= witness_margin`` for every design, and the gap is inert.

    Not a census but a theorem: ``drift_radius`` is the largest margin whose
    tolerance is still non-negative, so it cannot sit below the margin the
    design already holds.  Swept rather than argued as well, because a design
    sweep is what caught the pieces this law is assembled from.

    The sweep also shows the gap is not rare -- 49 of its 186 designs reach
    further than they hold -- and that none of that reach buys drift: past the
    witness margin every one of them has tolerance zero.
    """

    DELTAS = (Fraction(1, 100), Fraction(1, 12), Fraction(1, 3))

    def designs(self):
        for period in (12, 18, 20):
            U = universe(period)
            gr = grids(period)
            masks = cover_masks(gr, U, cell_widths(gr, U), Fraction(0))
            full = (1 << len(U)) - 1
            for G in all_designs(gr):
                reached = 0
                for M in G:
                    reached |= masks[M]
                if reached != full:
                    continue
                for k in (2, 3, 5):
                    Q = tuple(sorted({Fraction((j * 7 + 3) % period, period)
                                      for j in range(k)}))
                    if Q:
                        yield period, U, G, Q

    def test_the_radius_is_never_below_the_witness_margin(self):
        readings = 0
        for period, U, G, Q in self.designs():
            w = min(witness_margin(G, Q, U).values())
            r = min(drift_radius(G, Q, U).values())
            self.assertGreaterEqual(r, w, f"T={period} G={G} Q={Q}")
            readings += 1
        self.assertGreater(readings, 150)

    def test_the_tolerance_dies_at_the_witness_margin(self):
        for period, U, G, Q in self.designs():
            w = min(witness_margin(G, Q, U).values())
            self.assertEqual(min(drift_tolerance(G, Q, U, w).values()),
                             Fraction(0), f"T={period} G={G} Q={Q}")

    def test_so_no_drift_keeps_more_than_the_witness_margin(self):
        for period, U, G, Q in self.designs():
            w = min(witness_margin(G, Q, U).values())
            for delta in self.DELTAS:
                self.assertLessEqual(guaranteed_margin(G, Q, U, delta), w,
                                     f"T={period} G={G} Q={Q} delta={delta}")

    def test_the_gap_is_real_and_unspendable(self):
        """A design that reaches further than it holds, measured once."""
        gap = 0
        for period, U, G, Q in self.designs():
            w = min(witness_margin(G, Q, U).values())
            if min(drift_radius(G, Q, U).values()) > w:
                gap += 1
        self.assertGreater(gap, 0, "the sweep stopped exercising the gap")


if __name__ == "__main__":
    unittest.main()
