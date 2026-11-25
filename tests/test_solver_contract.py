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

import unittest
from fractions import Fraction

from exact_lattice_quantiles.ceiling import instance_ceiling
from exact_lattice_quantiles.solve import solve_lattice_schedule

EVEN_PERIODS = tuple(range(2, 401, 2))


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


if __name__ == "__main__":
    unittest.main()

