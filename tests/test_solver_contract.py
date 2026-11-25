"""The two margins a Schedule reports have to nest.

``design_margin`` is what the design guarantees across the universe;
``witness_margin`` is what the chosen sections actually hold.  A witness above
the design is an over-promise, and it is what the instance ceiling exists to
bound.

Issue #6 recorded the nesting failing at 96 of the 200 even periods in
``2..400``: ``_interval_map`` normalised blocks by ``M // 2`` while
``cell_width`` normalised them by ``M / 2``, and the two agree only on even
grids.  The repair puts both in the period currency, and this file pins the
result rather than the truce: the exempt set is now empty, and it stays empty.
The ceiling itself was never in doubt -- no design beats it at any period
measured here, in either currency.
"""

from __future__ import annotations

import math
import unittest
from fractions import Fraction

from exact_lattice_quantiles import cell_width, seats
from exact_lattice_quantiles.ceiling import instance_ceiling
from exact_lattice_quantiles.solve import solve_lattice_schedule

EVEN_PERIODS = tuple(range(2, 401, 2))

ODD_GRIDS = (15, 21, 27, 33, 45, 63, 75, 99, 105)
EVEN_GRIDS = (12, 18, 20, 30, 42, 60)


def _phase_density(M: int) -> dict[int, Fraction]:
    """``d -> share of the M phases that land in ``d``, counted from the gcd.

    Written from the definition of a class and nothing else: phase ``j`` belongs
    to ``M / gcd(j, M)``.  It never consults the seat map, so it can disagree
    with :func:`cell_width` -- which is the whole point.
    """
    count: dict[int, int] = {}
    for j in range(M):
        d = M // math.gcd(j, M)
        count[d] = count.get(d, 0) + 1
    return {d: Fraction(c, M) for d, c in count.items()}


def _schedule(period: int):
    return solve_lattice_schedule(period, max_cost=10**6, k_probes=10**3,
                                  min_margin=Fraction(0))


class TestSolverContract(unittest.TestCase):

    def test_no_design_beats_the_instance_ceiling(self):
        for T in EVEN_PERIODS:
            s = _schedule(T)
            self.assertLessEqual(s.design_margin, instance_ceiling(T).margin, T)

    def test_drift_radius_never_exceeds_the_witness_it_holds(self):
        for T in EVEN_PERIODS:
            s = _schedule(T)
            self.assertLessEqual(s.robust_margin, s.witness_margin, T)

    def test_the_witness_never_exceeds_the_design(self):
        above = []
        for T in EVEN_PERIODS:
            s = _schedule(T)
            if s.witness_margin > s.design_margin:
                above.append((T, str(s.design_margin), str(s.witness_margin)))
        self.assertEqual(above, [], "a witness above its design is an over-promise")

    def test_the_ceiling_refusal_agrees_with_the_witness(self):
        # T = 6 was the worst ratio of issue #6: the design guaranteed 2/3 while
        # the witness reported 1, so asking for 1 was refused as impossible by a
        # call that had just returned it.  The witness now stops at the ceiling.
        s = _schedule(6)
        self.assertEqual(s.design_margin, Fraction(2, 3))
        self.assertEqual(s.witness_margin, Fraction(2, 3))
        self.assertEqual(s.ceiling.margin, Fraction(2, 3))
        with self.assertRaises(Exception) as caught:
            solve_lattice_schedule(6, max_cost=10**6, k_probes=10**3,
                                   min_margin=Fraction(1))
        self.assertEqual(caught.exception.rc, 1)


class TestTheCurrencyHasAnIndependentReferee(unittest.TestCase):
    """The width of a block, checked against phases instead of against seats.

    Issue #6 was found by two functions of this package disagreeing with each
    other.  That is enough to prove one of them wrong, not to say which.  The
    referee here counts phases from ``M / gcd(j, M)``: on an odd grid it agrees
    with :func:`cell_width` at every divisor, so the ``a / (M // 2)`` reading
    -- the one that stretches an odd grid by ``M / (M - 1)`` and tiles every
    grid to exactly ``1`` -- is wrong, and only a shared-section argument can
    see the damage.
    """

    def test_an_odd_grid_width_matches_the_phase_density(self):
        for M in ODD_GRIDS:
            density = _phase_density(M)
            for d, _, _ in seats(M):
                self.assertEqual(cell_width(M, d), density[d], f"M={M} d={d}")

    def test_except_the_antipode_an_even_grid_also_matches(self):
        # On an even grid exactly one class is its own antipode: the phase
        # M // 2.  The seat convention books it as a whole class, so it is the
        # single divisor where the seat map and the density disagree -- by one
        # class, in both implementations alike.  Pinned so that a later change
        # of convention has to be an announced change of currency.
        for M in EVEN_GRIDS:
            density = _phase_density(M)
            off = [d for d, _, _ in seats(M) if cell_width(M, d) != density[d]]
            self.assertEqual(len(off), 1, f"M={M}: {off}")
            # The odd one out is the antipodal class, and it is booked wider.
            antipodal = M // math.gcd(M // 2, M)
            self.assertEqual(off, [antipodal], f"M={M}")
            self.assertEqual(cell_width(M, antipodal) - density[antipodal],
                             Fraction(1, M))

    def test_an_odd_grid_leaves_exactly_the_carrier_uncovered(self):
        # This is what makes the wrong currency invisible: in ``a / (M // 2)``
        # the blocks tile the section circle to exactly 1 for every grid, so
        # the measure looks self-evidently right.  In the period currency an odd
        # grid leaves exactly 1 / M uncovered -- the phase that is exactly zero,
        # which no section occupies.
        for M in ODD_GRIDS:
            covered = sum(cell_width(M, d) for d, _, _ in seats(M))
            self.assertEqual(covered, 1 - Fraction(1, M), f"M={M}")
        for M in EVEN_GRIDS:
            covered = sum(cell_width(M, d) for d, _, _ in seats(M))
            self.assertEqual(covered, Fraction(1), f"M={M}")


if __name__ == "__main__":
    unittest.main()

