"""The fee frontier is a law: two parts cost two, and a fee of one needs a triple.

``setup_fee_frontier`` walks the fee upwards and reports where the closed form
first breaks, one period at a time.  ``fee_frontier_census`` answers the same
question without walking: the fee at which a request breaks is ``block_fee`` of
an antichain, so the frontier of a period is a minimum over its antichains -- and
the claim under test is that this minimum is always attained with two or three
parts.  Only an antichain can break the formula, because it prices the maximals
of a request and a part that divides another is never priced apart.

The three states are kept apart: ``frontier = None`` means the universe is a
chain and no fee ever breaks it (a positive answer, not a zero), ``0`` means the
formula already misprices with no fee charged, and ``RC_NO_MEDIDO`` guards the
case where the law was computed but checked nowhere.
"""

import contextlib
import io
import json
import unittest
from unittest import mock

from exact_lattice_quantiles import (block_advantage, block_fee,
                                     fee_frontier_by_blocks,
                                     fee_frontier_census)
from exact_lattice_quantiles.cli import main


class TheFeeOfABlock(unittest.TestCase):
    def test_a_pair_costs_two_and_attains_it(self):
        self.assertEqual(block_advantage((2, 3)), -1)
        self.assertEqual(block_fee((2, 3)), 2)

    def test_a_triple_can_cost_only_one(self):
        self.assertEqual(block_advantage((6, 10, 15)), 0)
        self.assertEqual(block_fee((6, 10, 15)), 1)
        self.assertEqual(block_fee((6, 14, 21)), 1)

    def test_a_chain_block_needs_no_fee_because_it_is_already_cheaper(self):
        self.assertGreater(block_advantage((4, 8)), 0)
        self.assertEqual(block_fee((4, 8)), 0)

    def test_a_block_of_one_grid_is_not_a_block(self):
        with self.assertRaises(ValueError):
            block_fee((6,))


class TheFrontierOfAPeriod(unittest.TestCase):
    def test_a_universe_that_is_a_chain_answers_none_not_zero(self):
        row = fee_frontier_by_blocks(4)
        self.assertTrue(row["chain"])
        self.assertIsNone(row["frontier"])
        self.assertEqual(row["witness"], ())

    def test_the_first_period_that_breaks_does_so_at_two(self):
        self.assertEqual(fee_frontier_by_blocks(6)["frontier"], 2)

    def test_and_the_first_period_needing_only_one_is_1080(self):
        self.assertEqual(fee_frontier_by_blocks(180)["frontier"], 1)
        self.assertEqual(fee_frontier_by_blocks(1080)["frontier"], 0)

    def test_an_arity_of_one_is_refused(self):
        with self.assertRaises(ValueError):
            fee_frontier_by_blocks(180, arity=1)


class TheCensus(unittest.TestCase):
    def setUp(self):
        self.census = fee_frontier_census(max_period=400)

    def test_the_law_agrees_with_every_antichain_where_both_were_computed(self):
        self.assertEqual(self.census["mismatches"], [])
        self.assertGreater(self.census["verified"], 150)

    def test_only_the_powers_of_two_below_400_are_chains(self):
        self.assertEqual(self.census["chains"], (2, 4, 8, 16, 32, 64, 128, 256))
        self.assertEqual(self.census["at_zero"], ())

    def test_the_two_shapes_that_break_at_a_fee_of_one(self):
        shapes = {witness for _, witness in self.census["at_one"]}
        self.assertEqual(shapes, {(6, 10, 15), (6, 14, 21)})

    def test_every_verified_witness_has_two_or_three_parts(self):
        self.assertEqual(sorted(self.census["verified_witness_sizes"]), [2, 3])

    def test_the_verdict_counts_the_periods_it_actually_checked(self):
        self.assertIn(f"all {self.census['verified']} periods",
                      self.census["verdict"])
        self.assertNotIn("NO_MEDIDO", self.census["verdict"])
        self.assertNotIn("REFUTED", self.census["verdict"])

    def test_restricting_the_law_to_pairs_refutes_it(self):
        pairs = fee_frontier_census(max_period=200, arity=2)
        self.assertNotEqual(pairs["mismatches"], [])
        self.assertIn("REFUTED", pairs["verdict"])
        self.assertEqual(
            [m["period"] for m in pairs["mismatches"]][:2], [60, 84])


class ThePairLaw(unittest.TestCase):
    """The third part buys one thing only: a fee of one.

    A pair is where the frontier lands whenever the fee reaches two, so the
    census could enumerate pairs instead of triples -- and the periods where it
    cannot are exactly the ones a triple is needed for, which is what makes the
    two statements one.
    """

    def test_a_triple_beats_every_pair_only_below_a_fee_of_two(self):
        row = fee_frontier_by_blocks(180)
        self.assertEqual((row["frontier"], row["pair_frontier"]), (1, 2))
        self.assertEqual(row["pair_witness"], (2, 3))

    def test_and_where_no_triple_wins_the_pair_is_the_frontier(self):
        row = fee_frontier_by_blocks(1994)
        self.assertEqual(row["frontier"], row["pair_frontier"])

    def test_a_chain_has_neither_frontier_nor_pair(self):
        row = fee_frontier_by_blocks(4)
        self.assertIsNone(row["pair_frontier"])
        self.assertEqual(row["pair_witness"], ())

    def test_the_periods_a_triple_wins_are_the_periods_at_a_fee_of_one(self):
        census = fee_frontier_census(max_period=400)
        self.assertEqual(census["pair_over"], ())
        self.assertEqual(len(census["pair_misses"]),
                         len(census["at_zero"]) + len(census["at_one"]))
        self.assertEqual(census["at_zero"], ())

    def test_the_law_names_the_fee_it_leaves_to_the_third_part(self):
        self.assertIn("the fee of one",
                      fee_frontier_census(max_period=400)["pair_law"])

    def test_at_two_parts_the_pair_law_declares_itself_unmeasured(self):
        row = fee_frontier_census(max_period=400, arity=2)
        self.assertTrue(row["pair_law"].startswith("NOT_APPLICABLE"),
                        row["pair_law"])
        self.assertEqual((row["pair_misses"], row["pair_over"]), ((), ()))


class TheCommandLine(unittest.TestCase):
    def run_cli(self, *argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(list(argv))
        return rc, out.getvalue(), err.getvalue()

    def test_the_command_prints_the_law_and_exits_zero(self):
        rc, out, _ = self.run_cli("fee-census", "--max-period", "200")
        self.assertEqual(rc, 0)
        self.assertIn("mismatches", out)
        self.assertIn("(6, 10, 15)", out)

    def test_the_pair_law_is_printed_beside_the_verdict(self):
        rc, out, _ = self.run_cli("fee-census", "--max-period", "200")
        self.assertEqual(rc, 0)
        self.assertIn("pairs miss", out)
        self.assertIn("pair law", out)

    def test_a_two_part_census_prints_no_pair_rows(self):
        rc, out, _ = self.run_cli("fee-census", "--max-period", "200",
                                  "--arity", "2")
        self.assertEqual(rc, 1)
        self.assertNotIn("pairs miss", out)
        self.assertIn("NOT_APPLICABLE", out)

    def test_a_refuted_law_exits_one(self):
        rc, _, err = self.run_cli("fee-census", "--max-period", "200",
                                  "--arity", "2")
        self.assertEqual(rc, 1)
        self.assertIn("REFUTED", err)

    def test_a_law_checked_nowhere_exits_no_medido(self):
        empty = {"periods": 1, "verified": 0, "mismatches": [], "chains": (),
                 "at_zero": (), "at_one": (), "rows": (),
                 "universe": "U(T)", "arity": 3, "verify_universe": 14,
                 "largest": None, "verified_witness_sizes": {},
                 "pair_misses": (), "pair_over": (), "pair_law": "",
                 "theorem": "", "verdict": "NO_MEDIDO"}
        with mock.patch("exact_lattice_quantiles.cli.fee_frontier_census",
                        return_value=empty):
            rc, _, err = self.run_cli("fee-census")
        self.assertEqual(rc, 4)
        self.assertIn("NO_MEDIDO", err)

    def test_json_publishes_one_row_per_period(self):
        rc, out, _ = self.run_cli("--json", "fee-census", "--max-period", "100")
        self.assertEqual(rc, 0)
        payload = json.loads(out)
        self.assertEqual(len(payload["rows"]), 50)
        self.assertEqual(payload["rows"][1]["period"], 4)
        self.assertTrue(payload["rows"][1]["chain"])

    def test_the_divisors_universe_is_selectable(self):
        rc, out, _ = self.run_cli("fee-census", "--max-period", "100",
                                  "--universe", "divisors")
        self.assertEqual(rc, 0)
        self.assertIn("D(T) universe", out)


if __name__ == "__main__":
    unittest.main()
