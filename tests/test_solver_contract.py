"""The two margins a Schedule reports have to nest.

``design_margin`` is what the design guarantees across the universe;
``witness_margin`` is what the chosen sections actually hold.  A witness above
the design is an over-promise, and it is what the instance ceiling exists to
bound.

Issue #6 records that the nesting fails today at 96 of the 200 even
periods in ``2..400``.  The set is pinned rather than tolerated: a fix has to
shrink it on purpose, and a regression cannot grow it unnoticed.  The ceiling
itself is not in doubt -- no design beats it at any period measured here.
"""

from __future__ import annotations

import unittest
from fractions import Fraction

from exact_lattice_quantiles.ceiling import instance_ceiling
from exact_lattice_quantiles.solve import solve_lattice_schedule

EVEN_PERIODS = tuple(range(2, 401, 2))

#: Periods where ``witness_margin > design_margin`` as of issue #6.
KNOWN_WITNESS_ABOVE_DESIGN = (
    6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50,
    54, 58, 62, 66, 70, 74, 78, 82, 86, 90, 94, 98,
    102, 106, 110, 114, 118, 122, 126, 130, 134, 138, 142, 146,
    150, 154, 158, 162, 166, 170, 174, 178, 182, 186, 190, 194,
    198, 202, 206, 214, 218, 222, 226, 230, 234, 238, 242, 246,
    250, 254, 258, 262, 266, 270, 274, 278, 282, 286, 290, 294,
    298, 302, 306, 310, 314, 318, 322, 326, 334, 338, 342, 346,
    350, 354, 358, 362, 366, 370, 374, 378, 382, 386, 394, 398,
)


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

    def test_the_witness_is_bounded_by_the_design_except_where_issue_6_says(self):
        measured = []
        for T in EVEN_PERIODS:
            s = _schedule(T)
            if s.witness_margin > s.design_margin:
                measured.append(T)
        measured = tuple(measured)
        self.assertEqual(measured, KNOWN_WITNESS_ABOVE_DESIGN,
                         "issue #6 changed: update the pinned set deliberately")

    def test_the_ceiling_refusal_can_contradict_the_witness_it_returns(self):
        # T = 6 is the worst ratio of issue #6: the design guarantees 2/3, the
        # witness reports 1, and asking for 1 is refused as impossible at any
        # cost.  Pinned so the day the two agree this test has to be rewritten.
        s = _schedule(6)
        self.assertEqual(s.design_margin, Fraction(2, 3))
        self.assertEqual(s.witness_margin, Fraction(1))
        self.assertEqual(s.ceiling.margin, Fraction(2, 3))
        with self.assertRaises(Exception) as caught:
            solve_lattice_schedule(6, max_cost=10**6, k_probes=10**3,
                                   min_margin=Fraction(1))
        self.assertEqual(caught.exception.rc, 1)


if __name__ == "__main__":
    unittest.main()
