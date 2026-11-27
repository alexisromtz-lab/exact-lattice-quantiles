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
from fractions import Fraction
from unittest import mock

from itertools import combinations

from exact_lattice_quantiles import (arity_ladder_census, block_advantage,
                                     block_defect, block_fee, block_fee_closed,
                                     block_free, carrier_cofactor_window,
                                     carrier_cofactors, carrier_image,
                                     carrier_image_by_law,
                                     carrier_image_census,
                                     carrier_law_census,
                                     carrier_universe_by_law, carriers,
                                     cofactor_block, crossing_scalars,
                                     euler_phi,
                                     fee_frontier_by_blocks,
                                     fee_frontier_census, free_margin,
                                     free_margin_census, free_margin_cross,
                                     grids, in_carrier_universe,
                                     scaling_census, scaling_margin, universe)
from exact_lattice_quantiles.cli import main
from exact_lattice_quantiles.pricing import _best_margin


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


class TheFreeMargin(unittest.TestCase):
    """The distance to a free block, in the integer the sign throws away.

    ``block_free`` is a predicate and a predicate cannot be read as a horizon:
    the census that asks "where does the floor stop being one" answers in fees,
    and a fee is a quotient.  The margin is the numerator, so the two questions
    -- is this block free, and how far is this period from having a free block
    -- part company exactly where a horizon needs them not to.
    """

    @classmethod
    def setUpClass(cls):
        cls.carriers = free_margin_census(max_period=1200)
        cls.divisors = free_margin_census(max_period=1200, reduced=False)

    def test_the_margin_decides_freeness_and_keeps_the_distance(self):
        self.assertEqual(free_margin((6, 10, 15)), 0)      # delta 1, one odd part
        self.assertFalse(block_free((6, 10, 15)))
        self.assertEqual(block_fee((6, 10, 15)), 1)
        self.assertEqual(free_margin((12, 20, 30)), 2)     # delta 2, no odd part
        self.assertTrue(block_free((12, 20, 30)))
        self.assertEqual(block_fee((12, 20, 30)), 0)

    def test_a_fee_is_the_margin_divided_away(self):
        """Fee 1 is charged for margin 0 and margin -2 alike: the sign's distance
        is invisible to a fee, which is why the horizon read on fees stopped
        short of 1080 with the witness in reach."""
        for period in range(2, 61, 2):
            for pool in (universe(period), grids(period)):
                for parts in self._antichains(pool):
                    margin = free_margin(parts)
                    self.assertEqual(margin % 2, 0)
                    self.assertEqual(block_free(parts), margin > 0)
                    self.assertEqual(block_fee(parts) == 0, margin > 0)

    def test_the_search_is_the_answer_on_every_small_universe(self):
        """The pruned walk is checked against enumeration, not trusted."""
        from exact_lattice_quantiles.pricing import _best_margin
        for period in range(2, 121, 2):
            for pool in (universe(period), grids(period)):
                truth = [free_margin(parts) for parts in self._antichains(pool)]
                margin, block = _best_margin(pool)
                if not truth:
                    self.assertIsNone(margin)
                    continue
                self.assertEqual(margin, max(truth), f"period {period}")
                self.assertEqual(free_margin(block), margin)

    def test_scaling_buys_the_crossing_and_never_a_fourth_part(self):
        """delta and the lcm scale, the count of odd parts does not, so
        margin(kA) = k * delta - (k even ? 0 : k(A)) and nine times the Giuga
        triple is free while the triple itself pays one."""
        giuga = (6, 10, 15)
        for k in (1, 2, 3, 4):
            scaled = tuple(k * part for part in giuga)
            self.assertEqual(free_margin(scaled), k if k % 2 == 0 else k - 1)
        self.assertEqual(len(self.carriers["first_free"]["block"]), 3)

    def test_only_four_scalings_enter_the_carriers_below_4000(self):
        """The README's lambda table: which multiples of the Giuga triple the
        carriers admit at all, and that each of them arrives free."""
        pools = [set(universe(period)) for period in range(2, 4001, 2)]
        admitted = [k for k in range(2, 23)
                    if any(all(k * part in pool for part in (6, 10, 15))
                           for pool in pools)]
        self.assertEqual(admitted, [9, 14, 18, 22])
        for k in admitted:
            self.assertGreater(free_margin((6 * k, 10 * k, 15 * k)), 0)

    def test_the_carriers_touch_the_line_and_cross_at_1080(self):
        first = self.carriers["first_free"]
        self.assertEqual((first["period"], first["block"], first["margin"]),
                         (1080, (54, 90, 135), 8))
        self.assertEqual(sorted({r["block"] for r in self.carriers["on_line"]}),
                         [(6, 10, 15)])
        self.assertEqual([r["period"] for r in self.carriers["on_line"]],
                         [p for p in range(60, 1201, 60) if p != 1080])

    def test_a_horizon_read_on_fees_is_eighty_periods_short(self):
        """The number that makes this an instrument and not a restatement: at
        1000 there is no free carrier block, at 1080 there is one."""
        self.assertIsNone(free_margin_census(max_period=1000)["first_free"])
        self.assertEqual(
            free_margin_census(max_period=1080)["first_free"]["period"], 1080)

    def test_the_divisors_zero_on_the_giuga_numbers(self):
        self.assertEqual([r["period"] for r in self.divisors["on_line"]],
                         [30, 858])
        self.assertEqual(self.divisors["first_free"]["period"], 60)

    def test_the_command_crosses_and_says_so(self):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(["free-margin"])
        self.assertEqual((rc, err.getvalue()), (0, ""))
        self.assertIn("1080", out.getvalue())
        self.assertIn("(54, 90, 135)", out.getvalue())

    def test_a_horizon_with_no_block_at_all_is_not_measured(self):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(["free-margin", "--max-period", "2"])
        self.assertEqual(rc, 4)
        self.assertIn("NO_MEDIDO", err.getvalue())

    @staticmethod
    def _antichains(pool):
        pool = tuple(sorted(set(pool)))
        for size in range(2, len(pool) + 1):
            for parts in combinations(pool, size):
                if not any(max(a, b) % min(a, b) == 0
                           for a, b in combinations(parts, 2)):
                    yield parts


class TheSecondReader(unittest.TestCase):
    """A number from one instrument is a claim about that instrument.

    ``1080`` comes out of a pruned depth-first walk, and the bound that makes
    the walk affordable is exactly where a bug would sit: a bound too small
    silently drops the crossing and ``free_margin_census`` reports a floor that
    the universe does not have.  ``free_margin_cross`` asks the same question
    again by grouping blocks by their lcm -- a decomposition with no objective
    bound at all -- and compares period by period.
    """

    def test_both_readers_agree_period_by_period_in_both_universes(self):
        for reduced in (True, False):
            row = free_margin_cross(max_period=300, reduced=reduced)
            self.assertEqual(row["disagreements"], ())
            self.assertEqual(row["broken_witnesses"], ())
            self.assertGreater(row["periods_measured"], 100)

    def test_the_second_reader_alone_reaches_the_crossing(self):
        row = free_margin_cross(max_period=1200)
        first = row["first_free"]
        self.assertEqual((first["period"], first["block"], first["margin"]),
                         (1080, (54, 90, 135), 8))
        self.assertEqual(free_margin(first["block"]), first["margin"])

    def test_the_divisor_line_to_2400_is_the_giuga_trio(self):
        """The prediction stated before the sweep: if the line's zeros are
        Giuga's equation, the third one is there and nothing else is."""
        self.assertEqual(
            free_margin_cross(max_period=2400, reduced=False)["on_line"],
            (30, 858, 1722))
        for primes in ((2, 3, 5), (2, 3, 11, 13), (2, 3, 7, 41)):
            self.assertEqual(free_margin(cofactor_block(primes)), 0)

    def test_the_crossings_to_3000_are_scales_and_nothing_wider(self):
        """The margin of ``k`` times the Giuga triple is ``k`` when ``k`` is
        even and ``k - 1`` when it is odd, so freeness is inherited by scaling
        and no wider block is needed to reach ``18`` at ``2160``."""
        crossed = free_margin_cross(max_period=3000)["crossed"]
        self.assertEqual([r["period"] for r in crossed], [1080, 2160, 2520])
        for (k, row) in zip((9, 18, 14), crossed):
            self.assertEqual(row["block"], (6 * k, 10 * k, 15 * k))
            self.assertEqual(row["margin"], k if k % 2 == 0 else k - 1)

    def test_a_reader_restricted_to_maximal_grids_loses_the_line(self):
        """The plausible simplification: only maximal grids are worth blocking.
        On ``U(60)`` that reads ``-20`` where the truth is ``0``, so the
        restriction -- not the arithmetic -- is what hides the crossing."""
        pool = universe(60)
        maxima = tuple(a for a in pool
                       if not any(b != a and b % a == 0 for b in pool))
        self.assertEqual(maxima, (6, 15, 20))
        self.assertEqual(free_margin(maxima), -20)
        self.assertEqual(free_margin_census(max_period=60)["on_line"][0]["block"],
                         (6, 10, 15))

    def test_the_command_reads_twice_and_says_so(self):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(["cross-read", "--max-period", "120"])
        self.assertEqual((rc, err.getvalue()), (0, ""))
        self.assertIn("by lcm", out.getvalue())
        self.assertIn("disagreements", out.getvalue())

    def test_a_horizon_with_no_block_at_all_is_not_measured(self):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(["cross-read", "--max-period", "2"])
        self.assertEqual(rc, 4)
        self.assertIn("NO_MEDIDO", err.getvalue())
        self.assertNotIn("the two readers agree", out.getvalue())


class TheCarrierLaw(unittest.TestCase):
    """``U(T)`` is a filter over divisors, not a census over grids.

    ``universe(T)`` builds every even grid of the period and unions their
    medians; that is a *forwards* definition, and forwards the answer costs
    ``T / 2`` seats per grid.  The law says the same set can be read
    *backwards*: a divisor ``a`` is in ``U(T)`` iff some even grid ``q * a``
    dividing ``T`` carries it, which is a question about the divisor and its
    cofactor alone.  The two routes must name the same set -- otherwise the law
    is a different universe wearing the same name.
    """

    def test_both_directions_name_the_same_universe(self):
        row = carrier_law_census(600)
        self.assertEqual(row["splits"], ())
        self.assertEqual(row["periods_measured"], 300)
        self.assertIn("same universe", row["verdict"])

    def test_the_divisor_is_in_the_universe_iff_a_grid_carries_it(self):
        for period in (60, 120, 180, 360, 1080, 2520):
            self.assertEqual(carrier_universe_by_law(period),
                             tuple(sorted(universe(period))))
        self.assertTrue(in_carrier_universe(6, 60))
        self.assertFalse(in_carrier_universe(30, 60))
        self.assertFalse(in_carrier_universe(7, 60))

    def test_a_divisor_no_grid_carries_is_in_no_universe_at_all(self):
        """The periods where the *divisor* universe touches the free line are
        ``30``, ``858`` and ``1722``; the first two are carried by no grid, so
        no ``U(T)`` ever holds them.  The carrier universe touches at ``60``,
        which is likewise uncarried -- the line is where a divisor is not."""
        image = carrier_image(900)
        self.assertIn(30, image["absent"])
        self.assertIn(858, image["absent"])
        self.assertIn(60, image["absent"])
        self.assertNotIn(6, image["absent"])
        self.assertNotIn(15, image["absent"])
        for period in (60, 120, 180, 300, 600):
            self.assertNotIn(30, universe(period))
            self.assertNotIn(60, universe(period))

    def test_the_cofactor_route_and_the_membership_test_agree(self):
        for a in (6, 10, 12, 15, 45, 90, 105, 504):
            self.assertEqual(carrier_cofactors(a, limit=16),
                             tuple(q for q in range(1, 17)
                                   if (q * a) % 2 == 0
                                   and a in carriers(q * a)))


class TheCarrierWindow(unittest.TestCase):
    """The cofactors that can carry ``a``, read off ``a`` and proved.

    ``carrier_image`` needed a grid limit because it reads the grids; the
    window makes that limit unnecessary, because whatever carries ``a`` has a
    cofactor inside a range that depends on ``a`` alone.  What is under test is
    that the range is a *theorem*: nothing outside it ever carries, and nothing
    inside it that carries is missed.
    """

    def test_the_window_is_the_sweep(self):
        for a in range(2, 241):
            self.assertEqual(
                carrier_cofactors(a),
                tuple(q for q in range(1, 61)
                      if (q * a) % 2 == 0 and a in carriers(q * a)),
                f"a = {a}")

    def test_the_bound_is_attained_so_a_alone_cannot_lower_it(self):
        """``a = 2`` needs cofactor ``2``, which is the ceiling ``2 phi(2) / 2``:
        the upper bound is reached, so no argument reading only ``a`` sharpens
        it.  Lower it and ``2`` leaves its own image."""
        floor, ceiling = carrier_cofactor_window(2)
        self.assertEqual((floor, ceiling), (1, 2))
        self.assertIn(ceiling, carrier_cofactors(2))
        self.assertIn(2, carriers(2 * ceiling))

    def test_a_power_of_two_carries_itself_at_cofactor_one(self):
        """The tie: when four divides the grid and the mass *meets* the middle,
        the divisor that straddles it carries too.  Reading the floor as a
        strict successor drops ``q = 1`` here, which is how this window was
        first written and wrong."""
        for a in (4, 8, 16, 32, 64, 128, 256, 512):
            self.assertEqual(carrier_cofactor_window(a)[0], 1)
            self.assertIn(1, carrier_cofactors(a))
            self.assertIn(a, carriers(a))

    def test_the_complement_is_not_sized_by_a_over_phi(self):
        """``a / phi(a)`` cannot name the uncarried divisors: ``210`` has ratio
        ``4.375`` and ``630`` carries it at cofactor ``3``, while ``30`` sits at
        ``3.75`` and no grid carries it at all.  The complement is a
        measurement, and this is the measurement refusing a closed form."""
        self.assertEqual(Fraction(210, euler_phi(210)), Fraction(35, 8))
        self.assertEqual(Fraction(30, euler_phi(30)), Fraction(15, 4))
        self.assertIn(210, carriers(3 * 210))
        self.assertIn(3, carrier_cofactors(210))
        self.assertNotIn(30, carrier_image_by_law(300)["least_cofactor"])
        self.assertTrue(all(a % 2 == 0
                            for a in carrier_image_by_law(300)["absent"]))

    def test_the_image_by_law_needs_no_limit_and_matches_the_grids(self):
        census = carrier_image_census(300)
        self.assertEqual(census["missed_by_law"], ())
        self.assertEqual(census["missed_by_grids"], ())
        self.assertEqual(census["by_law"]["absent"],
                         carrier_image(300)["absent"])
        self.assertIn("no limit", census["verdict"])

    def test_a_short_window_is_reported_as_the_theorem_failing(self):
        """The falsifier, mechanised: pretend the window misses the last
        cofactor and the census must say *the window is short*, not *the limit
        is short* -- one of those two is the law being wrong."""
        real = carrier_image_by_law

        def crippled(max_divisor: int = 1000) -> dict:
            full = real(max_divisor)
            least = dict(full["least_cofactor"])
            least.pop(6)
            return {**full, "least_cofactor": least,
                    "absent": full["absent"] + (6,),
                    "carriers_measured": full["carriers_measured"] - 1}

        with mock.patch("exact_lattice_quantiles.pricing.carrier_image_by_law",
                        crippled):
            census = carrier_image_census(60)
        self.assertEqual(census["missed_by_law"], (6,))
        self.assertEqual(census["missed_by_grids"], ())
        self.assertIn("theorem failing", census["verdict"])

    def test_a_short_grid_limit_is_reported_as_the_measurement_short(self):
        """And the other direction: a law-correct set against grids read only
        to the largest divisor must read as the *limit* being short, because
        ``180`` is carried at ``720`` and the grids never saw it."""
        census = carrier_image_census(300, grid_limit=300)
        self.assertEqual(census["missed_by_law"], ())
        self.assertIn(180, census["missed_by_grids"])
        self.assertIn("grid limit is short", census["verdict"])


class TheScalarLaw(unittest.TestCase):
    """Which multiples of ``(6, 10, 15)`` clear the free line, and when.

    The margin is not the question: ``margin(kA) = k - [k odd]``, positive for
    every ``k >= 2``.  What decides a crossing is whether the three scaled
    parts can all land in the same ``U(T)`` -- a membership question, answered
    by the carrier grids, and answered without sweeping a single period.
    """

    def test_the_scaled_margin_is_the_scale_minus_its_parity(self):
        base = (6, 10, 15)
        for k in range(1, 41):
            block = tuple(k * a for a in base)
            self.assertEqual(free_margin(block), k - (k % 2))
            self.assertEqual(scaling_margin(base, k), free_margin(block))

    def test_the_law_forecasts_every_crossing_to_4001(self):
        rows = [row for row in crossing_scalars(4001) if row["margin"] > 0]
        self.assertEqual([(row["k"], row["period"], row["margin"])
                          for row in rows],
                         [(9, 1080, 8), (14, 2520, 14), (18, 2160, 18),
                          (22, 3960, 22), (27, 3240, 26)])
        for row in rows:
            self.assertEqual(row["margin"], free_margin(row["block"]))
            self.assertEqual(row["margin"], row["margin_by_law"])

    def test_the_law_names_a_crossing_no_sweep_to_3000_could_see(self):
        """``k = 27`` crosses at ``3240``, past the horizon every census before
        v0.23.0 used.  The law reads it off the grids of the scaled parts; the
        universe of ``3240`` then confirms it with both readers."""
        rows = {row["k"]: row for row in crossing_scalars(4001)}
        self.assertNotIn(27, {row["k"] for row in crossing_scalars(3000)})
        self.assertEqual(rows[27]["period"], 3240)
        self.assertEqual(rows[27]["block"], (162, 270, 405))
        self.assertEqual(rows[27]["margin"], 26)
        self.assertEqual(free_margin(rows[27]["block"]), 26)
        self.assertEqual(_best_margin(universe(3240))[0], 26)

    def test_the_census_of_the_universe_finds_no_sixth_crossing(self):
        row = scaling_census(1200)
        self.assertEqual(row["reader_splits"], ())
        self.assertEqual([r["period"] for r in row["forecast"]],
                         [r["period"] for r in row["measured"]])
        self.assertEqual([r["period"] for r in row["measured"]], [1080])

    def test_a_scale_that_cannot_land_is_not_a_crossing(self):
        """``k = 5`` has margin ``4`` and is still no crossing: ``6 * 5 = 30``
        is carried by no even grid, so no period holds the scaled block and
        the margin is a number about a block that cannot be built."""
        self.assertEqual(scaling_margin((6, 10, 15), 5), 4)
        self.assertIn(30, carrier_image(100)["absent"])
        self.assertNotIn(5, {row["k"] for row in crossing_scalars(4001)})

    def test_the_command_reads_the_universe_both_ways_and_says_so(self):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(["carrier-law", "--max-period", "120",
                       "--max-divisor", "100"])
        self.assertEqual((rc, err.getvalue()), (0, ""))
        self.assertIn("periods read twice", out.getvalue())
        self.assertIn("no grid carries", out.getvalue())


class TheScalarLawCommand(unittest.TestCase):
    def test_the_command_names_the_scale_and_the_period(self):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(["scalar-law", "--max-period", "1200"])
        self.assertEqual((rc, err.getvalue()), (0, ""))
        self.assertIn("1080 (k=9, margin 8)", out.getvalue())
        self.assertIn("same 1 crossing periods", out.getvalue())

    def test_a_json_run_carries_the_forecast_and_the_census(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = main(["scalar-law", "--max-period", "240", "--json"])
        row = json.loads(out.getvalue())
        self.assertEqual(rc, 0)
        self.assertEqual([r["period"] for r in row["forecast"]], [])
        self.assertEqual([r["period"] for r in row["measured"]], [])
