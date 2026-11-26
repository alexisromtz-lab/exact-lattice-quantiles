"""The closed form's unstated hypothesis: a grid costs ``grid // 2`` and nothing else.

Every census the package published before this file priced a grid at half its
size.  Charge a fixed cost per grid bought -- a setup, a minimum billing unit --
and the closed form, which buys one grid per maximal divisor and never buys a
block, is no longer the price.  This file pins the identity that decides when::

    advantage(B, setup) = setup * (len(B) - 1) + sum(b // 2) - lcm(B) // 2

Folding ``k`` grids into their lcm refunds ``k - 1`` fees, so the fee **favours
the big block**.  On an incomparable pair ``lcm(a, b) >= a + b`` makes the
advantage at most ``-1``, so no pair can win below a fee of two, and the first
break is a triple whose zero-fee advantage was exactly zero: it tied, and the
two refunded fees turned the tie into a saving.

The last class keeps the boundary honest on both sides: the command returns
``RC_NO_MEDIDO`` when the horizon never reaches a break, and ``frontier = 0``
when the universe was already broken before any fee was charged.
"""

from __future__ import annotations

import contextlib
import io
import json
import unittest
from itertools import combinations
from math import lcm

from exact_lattice_quantiles import (block_advantage, price_exact,
                                     price_form_census, price_one_per_maximal,
                                     setup_fee_frontier)
from exact_lattice_quantiles.cli import main
from exact_lattice_quantiles.pricing import divisors


class TheIdentityThatDecidesIt(unittest.TestCase):
    def test_the_advantage_is_the_closed_form_minus_the_lcm(self):
        for parts in ((6, 10, 15), (12, 20, 30), (4, 6), (30, 42, 70)):
            for setup in (0, 1, 2, 5):
                self.assertEqual(
                    block_advantage(parts, setup),
                    sum(setup + p // 2 for p in parts)
                    - (setup + lcm(*parts) // 2))

    def test_each_folded_grid_refunds_exactly_one_fee(self):
        for parts in ((6, 10, 15), (12, 20, 30), (4, 6)):
            base = block_advantage(parts)
            for setup in (0, 1, 3, 7):
                self.assertEqual(
                    block_advantage(parts, setup),
                    base + setup * (len(parts) - 1))

    def test_the_exact_price_charges_the_fee_once_per_grid(self):
        self.assertEqual(price_exact(180, (6, 10, 15)), 15)
        self.assertEqual(price_exact(180, (6, 10, 15), setup=1), 16)
        self.assertEqual(price_one_per_maximal((6, 10, 15)), 15)
        self.assertEqual(price_one_per_maximal((6, 10, 15), setup=1), 18)

    def test_a_zero_fee_is_the_published_census(self):
        self.assertEqual(price_form_census(180).discrepancies, 0)
        self.assertEqual(price_form_census(180, setup=0).as_dict(),
                         price_form_census(180).as_dict())


class NoPairExplainsTheFirstBreak(unittest.TestCase):
    def test_no_incomparable_pair_wins_below_a_fee_of_two(self):
        grid = [d for d in divisors(180) if d >= 2]
        pairs = [(a, b) for a, b in combinations(grid, 2)
                 if b % a and a % b]
        self.assertTrue(pairs)
        for a, b in pairs:
            self.assertLessEqual(block_advantage((a, b)), -1)
            self.assertLessEqual(block_advantage((a, b), 1), 0)
        self.assertLessEqual(max(block_advantage((a, b), 1) for a, b in pairs),
                             0)

    def test_the_pair_bound_is_attained_and_therefore_not_slack(self):
        self.assertEqual(block_advantage((2, 3)), -1)
        self.assertEqual(block_advantage((2, 3), 1), 0)
        self.assertEqual(block_advantage((2, 3), 2), 1)

    def test_a_pair_can_win_from_a_fee_of_two(self):
        self.assertEqual(block_advantage((4, 6)), -1)
        self.assertEqual(block_advantage((4, 6), 1), 0)
        self.assertEqual(block_advantage((4, 6), 2), 1)

    def test_the_census_agrees_the_first_break_has_no_pair_in_it(self):
        row = setup_fee_frontier(180, max_setup=2)
        self.assertEqual(row["rows"][1].pair_failures, 0)
        self.assertEqual(row["pair_frontier"], 2)
        self.assertEqual(row["frontier"], 1)


class TheFrontierOfAnInstance(unittest.TestCase):
    def test_the_formula_stops_being_the_price_at_a_fee_of_one(self):
        row = setup_fee_frontier(180, max_setup=2)
        self.assertEqual([r.discrepancies for r in row["rows"]], [0, 32, 146])
        self.assertEqual(row["frontier"], 1)
        self.assertEqual(row["witness"], (6, 10, 15))
        self.assertEqual(row["consumers"], 4095)

    def test_the_witness_gap_grows_by_the_refunded_fees(self):
        row = setup_fee_frontier(180, max_setup=4)
        gaps = [r.witness_gap for r in row["rows"][1:]]
        self.assertEqual(gaps, [2, 4, 6, 8])
        self.assertEqual(gaps[1] - gaps[0], len(row["witness"]) - 1)

    def test_a_universe_already_broken_reports_a_zero_frontier(self):
        row = setup_fee_frontier(180, reduced=False, max_setup=1)
        self.assertEqual(row["frontier"], 0)
        self.assertEqual(row["witness"], (36, 60, 90))
        self.assertEqual(row["rows"][0].discrepancies, 9984)

    def test_the_verdict_names_the_witness_and_both_prices(self):
        row = setup_fee_frontier(180, max_setup=1)
        self.assertIn("(6, 10, 15)", row["verdict"])
        self.assertIn("16 exact", row["verdict"])
        self.assertIn("18 by the formula", row["verdict"])

    def test_a_horizon_beyond_the_break_is_not_a_frontier(self):
        row = setup_fee_frontier(180, max_setup=0)
        self.assertIsNone(row["frontier"])
        self.assertIn("NO_MEDIDO", row["verdict"])


class TheCommandLine(unittest.TestCase):
    def run_cli(self, *argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(list(argv))
        return rc, out.getvalue(), err.getvalue()

    def test_the_command_prints_the_frontier_and_exits_zero(self):
        rc, out, _ = self.run_cli("setup-fee", "--max-setup", "1")
        self.assertEqual(rc, 0)
        self.assertIn("frontier", out)
        self.assertIn("(6, 10, 15)", out)

    def test_a_horizon_that_never_breaks_exits_no_medido(self):
        rc, _, err = self.run_cli("setup-fee", "--max-setup", "0")
        self.assertEqual(rc, 4)
        self.assertIn("NO_MEDIDO", err)

    def test_json_publishes_one_row_per_fee(self):
        rc, out, _ = self.run_cli("--json", "setup-fee", "--max-setup", "1")
        self.assertEqual(rc, 0)
        payload = json.loads(out)
        self.assertEqual([r["setup"] for r in payload["rows"]], [0, 1])
        self.assertEqual(payload["frontier"], 1)
        self.assertEqual(payload["witness"], [6, 10, 15])

    def test_the_divisors_universe_is_selectable(self):
        rc, out, _ = self.run_cli("setup-fee", "--universe", "divisors",
                                  "--max-setup", "1")
        self.assertEqual(rc, 0)
        self.assertIn("divisors universe", out)


if __name__ == "__main__":
    unittest.main()
