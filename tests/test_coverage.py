"""Coverage duality: the same margin counted in grids, in sections, and in the
distance to the ceiling."""

from __future__ import annotations

import unittest
from fractions import Fraction
from itertools import combinations

from exact_lattice_quantiles import (best_probes, design_margin, drift_radius,
                                     drift_tolerance, grids, ladder,
                                     margin_ladder, minimum_cost_design,
                                     minimum_probes, universe, witness_margin)
from exact_lattice_quantiles.coverage import coverage_census
from exact_lattice_quantiles.lattice import (cell_intervals, cell_width,
                                             cell_widths)

PAIRS = ((12, 18, 20, 45), (4, 6, 10, 12, 15, 18, 20, 45))
COSTAR, COST64 = 47, 64


def _widths(period: int):
    return cell_widths(grids(period), universe(period))


class TestWidthIsAMassRatio(unittest.TestCase):
    def test_width_is_mass_over_grid(self):
        from exact_lattice_quantiles import mass

        for period in (60, 180):
            for (M, d), w in _widths(period).items():
                self.assertEqual(w, Fraction(mass(d), M), (M, d))

    def test_widest_grid_of_a_divisor_is_itself(self):
        for period in (60, 180):
            for d in universe(period):
                widest = max((cell_width(M, d), M) for M in grids(period)
                             if M % d == 0)
                self.assertEqual(widest[0], cell_width(d, d), d)


class TestLadderOf180(unittest.TestCase):
    def test_the_four_quoted_points(self):
        self.assertEqual(ladder(180), {
            COSTAR: Fraction(1, 6),
            52: Fraction(8, 45),
            59: Fraction(1, 5),
            COST64: Fraction(1, 3),
        })

    def test_ladder_is_monotone_in_cost(self):
        curve = ladder(180)
        costs = sorted(curve)
        for small, big in zip(costs, costs[1:]):
            self.assertLessEqual(curve[small], curve[big], (small, big))

    def test_ladder_never_passes_the_ceiling(self):
        from exact_lattice_quantiles import instance_ceiling

        ceiling = instance_ceiling(180).margin
        self.assertEqual(max(ladder(180).values()), ceiling)


class TestDesigns(unittest.TestCase):
    def test_minimum_design_of_the_cheapest_point(self):
        cost, design = minimum_cost_design(
            grids(180), universe(180), _widths(180), Fraction(1, 6))
        self.assertEqual((cost, design), (47, (12, 18, 20, 45)))

    def test_the_64_point_design_is_the_ceiling(self):
        cost, design = minimum_cost_design(
            grids(180), universe(180), _widths(180), Fraction(1, 3))
        self.assertEqual((cost, design), (64, (4, 6, 10, 12, 15, 18, 20, 45)))

    def test_dp_agrees_with_brute_force(self):
        # The DP is the instrument; the enumeration is the control.  Without
        # this check a wrong parent pointer reads as a cheaper design.
        for period in (12, 24, 60):
            gr, U, w = grids(period), universe(period), _widths(period)
            for v in sorted(set(w.values())):
                dp_cost, _ = minimum_cost_design(gr, U, w, v)
                brute = None
                for r in range(1, len(gr) + 1):
                    for design in combinations(gr, r):
                        if design_margin(design, U, w) >= v:
                            cost = sum(M // 2 for M in design)
                            if brute is None or cost < brute:
                                brute = cost
                self.assertEqual(dp_cost, brute, (period, v))

    def test_uncovered_design_has_margin_zero(self):
        w = _widths(180)
        self.assertEqual(design_margin((2,), universe(180), w), Fraction(0))


class TestProbes(unittest.TestCase):
    def test_probes_needed_for_the_two_quoted_ceilings(self):
        U = universe(180)
        self.assertEqual(minimum_probes(PAIRS[0], U, Fraction(1, 6))[0], 8)
        self.assertEqual(minimum_probes(PAIRS[1], U, Fraction(1, 3))[0], 4)

    def test_k_min_is_not_the_design_cost(self):
        # 47 points of grids need 8 sections; 64 points need 4.  The two
        # budgets run in opposite directions, which is the third currency.
        self.assertGreater(minimum_probes(PAIRS[0], universe(180), Fraction(1, 6))[0],
                           minimum_probes(PAIRS[1], universe(180), Fraction(1, 3))[0])

    def test_four_sections_on_the_cheap_design(self):
        mu, Q = best_probes(PAIRS[0], universe(180), 4)
        self.assertEqual(mu, Fraction(2, 33))
        self.assertEqual(len(Q), 4)
        self.assertEqual(min(witness_margin(PAIRS[0], Q, universe(180)).values()), mu)

    def test_eight_sections_reach_the_design_margin(self):
        mu, Q = best_probes(PAIRS[0], universe(180), 8)
        self.assertEqual(mu, Fraction(1, 6))
        self.assertEqual(min(witness_margin(PAIRS[0], Q, universe(180)).values()), mu)

    def test_four_sections_on_the_rich_design_reach_the_ceiling(self):
        mu, Q = best_probes(PAIRS[1], universe(180), 4)
        self.assertEqual(mu, Fraction(1, 3))
        self.assertEqual(min(witness_margin(PAIRS[1], Q, universe(180)).values()), mu)

    def test_no_level_holds_a_divisor_without_a_band(self):
        # The best this design can hold is 1/6, so a hair above it leaves some
        # divisor with no band at all: infeasible for every k, not a budget.
        self.assertEqual(minimum_probes(PAIRS[0], universe(180),
                                        Fraction(1, 6) + Fraction(1, 1000))[0], None)


class TestDrift(unittest.TestCase):
    def test_the_common_drift_radius_is_at_least_the_witness_slack(self):
        # Capture is a weaker condition than slack -- a section may leave the
        # block it sat in and be pushed into another one -- so the radius can
        # only be larger.  A test asserting the opposite asserts a falsehood.
        U = universe(180)
        for G in PAIRS:
            for k in (2, 4, 6):
                mu, Q = best_probes(G, U, k)
                if not Q:
                    continue
                drift = min(drift_radius(G, Q, U).values())
                witness = min(witness_margin(G, Q, U).values())
                self.assertGreaterEqual(drift, witness, (G, k))

    def test_the_radius_is_the_component_and_not_two_maxima(self):
        # This design has a section far from every block it captures, and drift
        # can push it in.  The version that only looks at pairs capturing d at
        # delta = 0 -- ``min(max(q - L), max(R - q))`` -- answers 1/6 here; the
        # exact component answers 37/210.
        U = universe(180)
        G = (2, 3, 4, 6, 9, 12, 15, 20, 36, 45, 180)
        _, Q = best_probes(G, U, 4)
        self.assertEqual(min(drift_tolerance(G, Q, U, Fraction(0)).values()),
                         Fraction(37, 210))
        peor = None
        for d in U:
            down = up = Fraction(0)
            for M in G:
                bloques = {x: (L, R) for x, L, R in cell_intervals(M)}
                if d not in bloques:
                    continue
                L, R = bloques[d]
                for q in Q:
                    if L < q <= R:
                        down = max(down, q - L)
                        up = max(up, R - q)
            alcanza = min(down, up) * 2
            peor = alcanza if peor is None else min(peor, alcanza)
        self.assertEqual(peor, Fraction(1, 6))

    def test_a_margin_of_zero_gives_the_coverage_radius(self):
        U = universe(180)
        for G in PAIRS:
            for k in range(1, 9):
                _, Q = best_probes(G, U, k)
                if not Q:
                    continue
                self.assertEqual(
                    min(drift_tolerance(G, Q, U, Fraction(0)).values()),
                    min(drift_radius(G, Q, U).values()), (G, k))

    def test_the_witness_margin_survives_no_drift_at_all(self):
        # Measured, not assumed: at the witness margin the binding divisor sits
        # on the edge of its block, so the tolerance is zero for every design
        # checked here.  This is why the curve exists and why the radius is not
        # a promise about the margin.
        U = universe(180)
        for G in PAIRS + ((2, 3, 4, 6, 9, 12, 15, 20, 36, 45, 180),
                          (4, 6, 12, 18, 20, 36, 45, 90, 180)):
            for k in (2, 4, 6):
                _, Q = best_probes(G, U, k)
                if not Q:
                    continue
                mu = min(witness_margin(G, Q, U).values())
                self.assertEqual(min(drift_tolerance(G, Q, U, mu).values()),
                                 Fraction(0), (G, k))

    def test_the_curve_of_the_ladder_designs_is_linear(self):
        # For the two designs the ladder prints the trade is one for one:
        # every unit of margin kept costs a unit of drift.  Measured over the
        # five shares; a design where it does not hold would fail here.
        U = universe(180)
        for G in PAIRS:
            for k in range(1, 9):
                _, Q = best_probes(G, U, k)
                if not Q:
                    continue
                mu = min(witness_margin(G, Q, U).values())
                for share in (Fraction(0), Fraction(1, 4), Fraction(1, 2),
                              Fraction(3, 4), Fraction(1)):
                    self.assertEqual(
                        min(drift_tolerance(G, Q, U, mu * share).values()),
                        mu * (1 - share), (G, k, share))

    def test_the_witness_set_of_p52_holds_its_radius(self):
        U = universe(180)
        mu, Q = best_probes(PAIRS[0], U, 8)
        self.assertEqual(min(drift_radius(PAIRS[0], Q, U).values()), mu)


class TestCensus(unittest.TestCase):
    def test_designs_of_180(self):
        data = coverage_census(180)
        self.assertEqual(data["designs"], 2 ** len(grids(180)))
        self.assertEqual(data["cheapest"], COSTAR)
        self.assertLess(data["covering"], data["designs"])

    def test_margin_ladder_matches_the_quoted_curve(self):
        gr, U, w = grids(180), universe(180), _widths(180)
        self.assertEqual(margin_ladder(gr, U, w), ladder(180))


if __name__ == "__main__":
    unittest.main()
