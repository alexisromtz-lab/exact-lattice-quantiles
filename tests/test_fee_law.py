"""The fee frontier is a law: two parts cost two, a triple one, and a fourth is free.

``setup_fee_frontier`` walks the fee upwards and reports where the closed form
first breaks, one period at a time.  ``fee_frontier_census`` answers the same
question without walking: the fee at which a request breaks is ``block_fee`` of
an antichain, so the frontier of a period is a minimum over its antichains --
and the first claim under test, that this minimum is always attained with two
or three parts, is the one the package published until v0.19.0 and its own
``TheFourthPart`` broke in v0.20.0.  Only an antichain can break the formula,
because it prices the maximals of a request and a part that divides another is
never priced apart.

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

from itertools import combinations

from exact_lattice_quantiles import (arity_ladder_census, block_advantage,
                                     block_defect, block_fee, block_fee_closed,
                                     block_free, cofactor_block,
                                     fee_frontier_by_blocks,
                                     fee_frontier_census, grids, universe)
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

    def test_the_law_is_the_same_in_both_universes(self):
        census = fee_frontier_census(max_period=200, reduced=False,
                                     verify_universe=16)
        self.assertEqual(census["pair_over"], ())
        self.assertEqual(len(census["pair_misses"]),
                         len(census["at_zero"]) + len(census["at_one"]))

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


class TheFeeIsTheDefect(unittest.TestCase):
    """``block_fee`` is one division on ``sum - lcm``: an identity, not a search."""

    def test_both_routes_price_every_block_the_same_in_both_universes(self):
        checked = 0
        for period in range(2, 201, 2):
            for U in (universe(period), grids(period)):
                for size in range(2, min(3, len(U)) + 1):
                    for parts in combinations(U, size):
                        self.assertEqual(block_fee(parts),
                                         block_fee_closed(parts), (period, parts))
                        checked += 1
        self.assertGreater(checked, 10_000)

    def test_the_blocks_that_are_free_are_exactly_the_ones_that_save(self):
        for period in range(2, 301, 2):
            for parts in combinations(universe(period), 3):
                self.assertEqual(block_fee_closed(parts) == 0,
                                 block_advantage(parts) > 0, (period, parts))

    def test_the_fee_of_one_is_a_block_one_short_of_perfect(self):
        self.assertEqual(block_defect((6, 10, 15)), 1)
        self.assertEqual(block_defect((6, 14, 21)), -1)
        self.assertEqual(block_fee_closed((6, 10, 15)), 1)
        self.assertEqual(block_fee_closed((6, 14, 21)), 1)

    def test_a_block_of_one_grid_is_refused(self):
        for route in (block_defect, block_fee_closed):
            with self.assertRaises(ValueError):
                route((6,))


class TheBlockFeeCommand(unittest.TestCase):
    def run_cli(self, *argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(list(argv))
        return rc, out.getvalue(), err.getvalue()

    def test_the_command_reads_the_fee_off_the_defect(self):
        rc, out, _ = self.run_cli("block-fee", "6", "10", "15")
        self.assertEqual(rc, 0)
        self.assertIn("defect", out)
        self.assertIn("by advantage", out)

    def test_json_carries_both_routes(self):
        rc, out, _ = self.run_cli("--json", "block-fee", "6", "21", "28")
        self.assertEqual(rc, 0)
        payload = json.loads(out)
        self.assertEqual((payload["defect"], payload["fee"]), (-29, 8))
        self.assertEqual(payload["fee_by_advantage"], payload["fee"])

    def test_one_grid_is_refused_rather_than_priced_at_zero(self):
        rc, out, err = self.run_cli("block-fee", "6")
        self.assertEqual(rc, 3)
        self.assertEqual(out, "")
        self.assertIn("NO BLOCK", err)


class TheFourthPart(unittest.TestCase):
    """The ladder of arities: a pair costs two, a triple one, a fourth is free.

    The claim the package published until v0.19.0 -- that the frontier is always
    attained with two or three parts -- is false, and what refutes it is a block
    of four parts that costs nothing where every triple of the same period costs
    one.  The census that missed it did not miss it by arithmetic: its
    verification horizon stopped one grid below the first witness.  So the tests
    here pin the witnesses, pin the horizon that reaches them, and keep the free
    rung's floor -- zero -- separate from the fee that attains it.
    """

    def test_free_means_the_parts_are_dearer_than_their_lcm(self):
        self.assertTrue(block_free((12, 20, 30)))       # 62 > 60, all parts even
        self.assertEqual(block_fee((12, 20, 30)), 0)
        self.assertFalse(block_free((6, 10, 15)))       # 31 > 30 but one part odd
        self.assertEqual(block_fee((6, 10, 15)), 1)
        self.assertTrue(block_free((42, 66, 154, 231)))  # 493 > 462
        self.assertEqual(block_fee((42, 66, 154, 231)), 0)

    def test_no_pair_is_ever_free(self):
        """``lcm`` doubles the larger part, so a pair's defect is never positive.

        This is why the ladder's second rung is not 1: a fee of two is the floor
        for two parts, in both universes, for every period measured.
        """
        checked = 0
        for period in range(6, 301, 6):
            for U in (universe(period), grids(period)):
                for parts in combinations(U, 2):
                    if max(parts) % min(parts) == 0:
                        continue
                    checked += 1
                    self.assertFalse(block_free(parts), (period, parts))
                    self.assertGreaterEqual(block_fee(parts), 2, (period, parts))
        self.assertGreater(checked, 1_000)

    def test_a_free_block_always_has_a_reciprocal_sum_above_one(self):
        """``free => sum(1/complement) > 1``, which a pair cannot satisfy."""
        from fractions import Fraction
        from math import lcm as _lcm
        checked = 0
        for period in range(2, 301, 2):
            for U in (universe(period), grids(period)):
                for size in (3, 4):
                    for parts in combinations(U, size):
                        if any(max(p) % min(p) == 0
                               for p in combinations(parts, 2)):
                            continue
                        checked += 1
                        if not block_free(parts):
                            continue
                        L = 1
                        for part in parts:
                            L = _lcm(L, part)
                        reciprocal = sum(Fraction(1, L // part)
                                         for part in parts)
                        self.assertGreater(reciprocal, 1, (period, parts))
        self.assertGreater(checked, 1_000)

    def test_the_fourth_part_is_free_at_462_and_the_triples_there_cost_one(self):
        three = fee_frontier_by_blocks(462, reduced=False, arity=3)
        four = fee_frontier_by_blocks(462, reduced=False, arity=4)
        self.assertEqual(three["frontier"], 1)
        self.assertEqual(four["frontier"], 0)
        self.assertEqual(four["witness"], (42, 66, 154, 231))
        self.assertEqual(three["universe_size"], 15)

    def test_the_first_free_block_arrives_with_three_parts_and_not_with_four(self):
        """The quadruple is not the first rung to reach zero: the triple is.

        At ``60`` the cofactors of ``2, 3, 5`` are all even, so ``k = 0`` and a
        defect of ``2`` already frees the block -- the same prime set that pays
        ``1`` at ``30``, where ``15`` is odd and the defect is ``1``.
        """
        ladder = arity_ladder_census(max_period=500)
        self.assertEqual(ladder["first_free"]["period"], 60)
        self.assertEqual(ladder["first_free"]["parts"], 3)
        self.assertEqual(ladder["first_free"]["witness"], (12, 20, 30))
        self.assertEqual(fee_frontier_by_blocks(30, reduced=False,
                                                arity=3)["frontier"], 1)

    def test_the_ladder_earns_at_exactly_these_periods_up_to_1000(self):
        ladder = arity_ladder_census(max_period=1000)
        self.assertEqual(ladder["earnings"]["3->4"],
                         (462, 546, 714, 798, 858, 924, 966))
        self.assertEqual(ladder["earnings"]["4->5"], ())
        self.assertEqual(ladder["first"]["2->3"]["witness"], (6, 10, 15))

    def test_the_ladder_is_the_frontier_wherever_it_could_be_checked(self):
        """``exhaustive_until=15`` is what reaches 462; 14 was one grid short."""
        ladder = arity_ladder_census(max_period=500, exhaustive_until=15)
        self.assertEqual(ladder["exhaustive"]["mismatches"], ())
        self.assertGreater(ladder["exhaustive"]["verified"], 100)
        short = arity_ladder_census(max_period=500, exhaustive_until=14)
        self.assertEqual(short["exhaustive"]["mismatches"], ())

    def test_over_the_carriers_no_period_reaches_zero(self):
        """The other universe has a floor of one, and the census says so in words."""
        ladder = arity_ladder_census(max_period=500, reduced=True)
        self.assertIsNone(ladder["first_free"])
        self.assertEqual(ladder["at_floor"], ())
        self.assertIn("NO_PERIOD REACHES ZERO", ladder["verdict"])

    def test_giugas_equation_is_one_defect_short_of_free(self):
        """``delta = 1`` is Giuga, and one defect past it is a fee of zero."""
        giuga = cofactor_block((2, 3, 7, 41))    # 1722 = sum of its cofactors + 1
        past = cofactor_block((2, 3, 7, 11))     # 462, defect 31
        short = cofactor_block((2, 3, 7, 43))    # primary pseudoperfect, defect -1
        self.assertEqual(block_defect(giuga), 1)
        self.assertEqual(block_fee(giuga), 1)
        self.assertFalse(block_free(giuga))
        self.assertEqual(block_defect(past), 31)
        self.assertTrue(block_free(past))
        self.assertEqual(block_defect(short), -1)
        self.assertEqual(block_fee(short), 1)

    def test_a_cofactor_block_needs_two_distinct_primes(self):
        with self.assertRaises(ValueError):
            cofactor_block((3,))
        with self.assertRaises(ValueError):
            cofactor_block((3, 3))
        self.assertEqual(cofactor_block((2, 3, 5)), (6, 10, 15))

    def test_the_command_reports_the_ladder_and_exits_clean(self):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(["arity-ladder", "--max-period", "500"])
        self.assertEqual((rc, err.getvalue()), (0, ""))
        self.assertIn("first free block", out.getvalue())
        self.assertIn("(42, 66, 154, 231)", err.getvalue() + out.getvalue())

    def test_json_carries_the_earnings_and_the_floor(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = main(["--json", "arity-ladder", "--max-period", "500"])
        self.assertEqual(rc, 0)
        payload = json.loads(out.getvalue())
        self.assertEqual(payload["earnings"]["3->4"], [462])
        self.assertEqual(payload["first_free"]["parts"], 3)
        self.assertEqual(payload["universe"], "divisors")

    def test_a_short_horizon_is_reported_as_not_measured(self):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(["arity-ladder", "--max-period", "8",
                       "--exhaustive-until", "2"])
        self.assertEqual(rc, 4)
        self.assertIn("NO_MEDIDO", err.getvalue())
