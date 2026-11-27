"""The cofactor window is necessary, and the census of what that buys.

:func:`carrier_cofactor_window` is a theorem: whatever carries ``a`` has a
cofactor in a range read off ``a`` alone, so the search for a carrier is a
search inside a window and never outside one.  The natural next sentence --
*in the window, therefore carried* -- is false, and it is false immediately:
``a = 30`` sits in ``[2, 18]`` and no even grid carries it at all, because
``carriers(30)`` is ``(15,)`` and the mass has already passed the middle by the
time it reaches ``30``.

Three things are pinned here.

* The window is never empty, and that is arithmetic and not luck:
  ``sum(phi(d) for d | a) == a`` keeps ``q = 2`` inside every one of them, so
  the emptiness count cannot be the evidence that necessity is doing work.
* The gap -- nonempty window, no cofactor in it carries -- is not empty either:
  ``45`` divisors to ``600``, beginning ``30``, and every one of them is a
  multiple of six at that horizon.
* A reader that trusts the window names the gap as carried, which is the
  mutant: the census counts the disagreement and the first one is ``30``.
"""

from __future__ import annotations

import io
import json
import unittest
from contextlib import redirect_stdout

from exact_lattice_quantiles.cli import RC_OK, main
from exact_lattice_quantiles.lattice import carriers
from exact_lattice_quantiles.pricing import (carrier_cofactor_window,
                                             carrier_cofactors,
                                             carrier_window_census)

HORIZON = 600


class TestTheWindowIsANecessaryCondition(unittest.TestCase):
    """What the theorem actually says, and what it does not."""

    def test_two_is_in_every_window_because_the_mass_sums_to_a(self):
        for a in range(2, 3001):
            floor, ceiling = carrier_cofactor_window(a)
            self.assertLessEqual(floor, 2, a)
            self.assertGreaterEqual(ceiling, 2, a)

    def test_the_census_finds_no_empty_window_and_no_window_without_two(self):
        census = carrier_window_census(HORIZON)
        self.assertEqual(census["empty_windows"], ())
        self.assertEqual(census["windows_without_two"], ())

    def test_thirty_sits_in_its_window_and_no_even_grid_carries_it(self):
        self.assertEqual(carrier_cofactor_window(30), (2, 18))
        self.assertEqual(carriers(30), (15,))
        self.assertEqual(carrier_cofactors(30), ())

    def test_the_window_is_not_a_decision_procedure(self):
        census = carrier_window_census(HORIZON)
        self.assertNotEqual(census["gap"], ())
        self.assertEqual(census["gap"][0], 30)
        self.assertEqual(len(census["gap"]),
                         HORIZON - 1 - census["carried"])


class TestTheGapIsTheEvensAlone(unittest.TestCase):
    """The shape of the failure, held to the horizon it was measured at."""

    def test_every_gap_divisor_is_a_multiple_of_six(self):
        census = carrier_window_census(HORIZON)
        self.assertEqual(census["gap_all_multiples_of_six"], True)
        self.assertEqual([a for a in census["gap"] if a % 6], [])

    def test_the_carried_divisors_need_no_cofactor_above_four(self):
        census = carrier_window_census(HORIZON)
        self.assertEqual(census["largest_least_cofactor"], 4)
        self.assertEqual(census["largest_witness"], 90)

    def test_the_widest_window_is_the_largest_divisor_that_cannot_be_carried(
            self):
        census = carrier_window_census(HORIZON)
        widest = census["widest_window"]
        self.assertEqual(widest["a"], 599)
        self.assertEqual((widest["floor"], widest["ceiling"]), (1, 365))

    def test_the_gap_is_disjoint_from_the_carried_set(self):
        census = carrier_window_census(HORIZON)
        carried = set(range(2, HORIZON + 1)) - set(census["gap"])
        self.assertEqual(len(carried), census["carried"])
        self.assertIn(90, carried)
        self.assertNotIn(30, carried)


class TestTheMutantIsTheBackwardsReading(unittest.TestCase):
    """Trust the window and the answer changes on exactly the gap."""

    def test_a_window_only_reader_carries_the_whole_gap(self):
        census = carrier_window_census(HORIZON)
        by_window = {a for a in range(2, HORIZON + 1)
                     if carrier_cofactor_window(a)[0]
                     <= carrier_cofactor_window(a)[1]}
        truly = set(range(2, HORIZON + 1)) - set(census["gap"])
        self.assertEqual(sorted(by_window - truly), list(census["gap"]))

    def test_the_verdict_names_the_gap_and_not_the_emptiness(self):
        verdict = carrier_window_census(HORIZON)["verdict"]
        self.assertIn("necessary and not sufficient", verdict)
        self.assertIn("multiple of six", verdict)


class TestTheCliCarriesTheCensus(unittest.TestCase):
    """The numbers reach the shell, and the exit code is not a failure code."""

    def test_carrier_law_prints_the_window_block(self):
        out = io.StringIO()
        with redirect_stdout(out):
            code = main(["carrier-law", "--max-period", "120",
                         "--max-divisor", "120"])
        text = out.getvalue()
        self.assertEqual(code, RC_OK)
        self.assertIn("window empty", text)
        self.assertIn("in window, uncarried", text)
        self.assertIn("gap all sixes", text)
        self.assertIn("widest window", text)

    def test_the_json_payload_is_well_formed(self):
        out = io.StringIO()
        with redirect_stdout(out):
            main(["carrier-law", "--max-period", "120", "--max-divisor", "120",
                  "--json"])
        payload = json.loads(out.getvalue())
        self.assertEqual(payload["cofactor_window"]["max_divisor"], 120)
        self.assertIn("gap", payload["cofactor_window"])


if __name__ == "__main__":
    unittest.main()
