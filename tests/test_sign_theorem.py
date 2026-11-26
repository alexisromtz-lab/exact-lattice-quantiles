"""The sign reading of the lemma, proved where it can be and banded where it cannot.

The law needs one direction of the lemma: *if no block saves, no antichain
saves*.  Over every period below ``2000`` that was a measurement.  It is now a
theorem on the region ``sum(1/p) <= 1``, and a measurement only on the band
where the block is empty-handed while the reciprocal sum is above one.

What makes the proof work is that a saving is not a sum but a **rounded** sum::

    saving(A) = (sum(a for a in A) - lcm(A) - odd_slack) / 2

with ``odd_slack`` at most one.  Beating zero therefore needs ``sum(a)`` to
exceed the period itself, and the lemma bounds ``sum(a)`` by the weight of the
drop-one block -- ``sum(period // p)`` -- which is at most the period exactly
when the reciprocal sum is at most one.  Integer arithmetic then takes the
last half point: something strictly below ``1/2`` and integral is at most zero.

The residue is where the reciprocal sum is above one and the block still saves
nothing.  That is not a gap in the search but an equation: it needs
``period * (sum(1/p) - 1) <= 1``, which for an even period is
``sum(1/p) = 1 + 1/period`` with the period twice an odd number.  Below
``2000`` exactly three periods satisfy it -- ``30``, ``858``, ``1722`` -- and
their equation is the Giuga equation, because for a squarefree period
``1/period`` is the product of the reciprocals.

These tests pin both halves, and pin that the instrument reports its own
refutation: the census returns a violation rather than a verdict when it finds
one, and the CLI turns that into ``RC 1``.
"""

from __future__ import annotations

import contextlib
import io
import unittest
from fractions import Fraction
from unittest import mock

from exact_lattice_quantiles.cli import main
from exact_lattice_quantiles.pricing import (antichain_saving_ceiling,
                                             block_saving,
                                             reciprocal_prime_sum, sign_regime,
                                             sign_theorem_census, two_adic,
                                             worst_antichain_weight)

BAND = (30, 858, 1722)


class TheReciprocalSum(unittest.TestCase):
    def test_it_is_exact(self):
        self.assertEqual(reciprocal_prime_sum(30), Fraction(31, 30))
        self.assertEqual(reciprocal_prime_sum(12), Fraction(5, 6))
        self.assertEqual(reciprocal_prime_sum(2310), Fraction(2927, 2310))

    def test_858_is_one_plus_its_own_reciprocal(self):
        self.assertEqual(reciprocal_prime_sum(858), 1 + Fraction(1, 858))

    def test_a_period_below_two_is_refused(self):
        with self.assertRaises(ValueError):
            reciprocal_prime_sum(1)

    def test_the_two_adic_valuation_is_what_the_band_turns_on(self):
        self.assertEqual(two_adic(858), 1)
        self.assertEqual(two_adic(60), 2)
        self.assertEqual(two_adic(45), 0)


class TheTheorem(unittest.TestCase):
    def test_no_antichain_saves_where_the_reciprocal_sum_is_at_most_one(self):
        census = sign_theorem_census(500)
        self.assertEqual(census["sign_violations"], [])
        self.assertGreater(census["proven_zero"], 0)
        self.assertEqual(census["proven_zero"],
                         sum(1 for n in range(2, 501)
                             if reciprocal_prime_sum(n) <= 1))

    def test_the_lemma_is_sharp_where_it_applies(self):
        """The block attains the bound, so nothing in the proof is slack."""
        for period in (12, 30, 36, 42, 66, 70):
            weight, witness = worst_antichain_weight(period)
            self.assertEqual(weight, sum(period // p
                                         for p in _primes(period)),
                             f"period {period}")
        self.assertEqual(worst_antichain_weight(36), (30, (18, 12)))

    def test_the_period_itself_is_not_an_antichain_it_is_a_grid(self):
        """``{period}`` weighs the period and saves nothing, so it is excluded.

        Admitting it would break the bound on every period with a reciprocal
        sum below one and say nothing about any cover.
        """
        weight, _ = worst_antichain_weight(858)
        self.assertEqual(weight, 859)
        self.assertEqual(sum(858 // p for p in _primes(858)), 859)
        self.assertEqual(block_saving(858), 0)

    def test_a_period_below_two_is_refused(self):
        with self.assertRaises(ValueError):
            worst_antichain_weight(1)


class TheBand(unittest.TestCase):
    def test_it_is_exactly_three_periods_below_2000(self):
        census = sign_theorem_census(2000)
        self.assertEqual([row["period"] for row in census["band_rows"]],
                         list(BAND))
        for row in census["band_rows"]:
            self.assertEqual(row["regime"], "band")
            self.assertEqual(row["block"], 0)
            self.assertEqual(row["excess"], "1")
            self.assertEqual(row["two_adic"], 1)
            self.assertTrue(row["band_equation"])
            self.assertTrue(row["characterised"])

    def test_nothing_saves_inside_it_but_the_proof_does_not_reach(self):
        for period in BAND:
            self.assertEqual(antichain_saving_ceiling(period), (0, ()))
            self.assertGreater(reciprocal_prime_sum(period), 1)

    def test_the_band_is_the_equation_and_not_a_window(self):
        """A band period needs ``sum(1/p) = 1 + 1/period``, twice an odd number.

        Both directions, over every even composite period below ``1000``: so an
        empty-handed block with a reciprocal sum above one cannot happen
        anywhere else, and the three members are a solution set, not the tail
        of a search.
        """
        for period in range(6, 1001, 2):
            if len(_primes(period)) < 2:
                continue
            equation = reciprocal_prime_sum(period) == 1 + Fraction(1, period)
            self.assertEqual(sign_regime(period)["regime"] == "band", equation,
                             f"period {period}")

    def test_they_are_the_giuga_equation(self):
        for period in BAND:
            product = Fraction(1)
            for prime in _primes(period):
                product *= Fraction(1, prime)
            self.assertEqual(reciprocal_prime_sum(period) - product, 1,
                             f"period {period}")

    def test_an_odd_period_is_a_window_and_not_an_open_question(self):
        # An odd band member is an odd Giuga number, open since 1950, so the
        # characterisation cannot be closed -- but the regime is no longer
        # unexplained either: it is the window 1 <= excess <= omega - 1.
        odd = sign_regime(15)
        self.assertEqual(odd["regime"], "proven-zero")
        self.assertTrue(odd["characterised"])
        self.assertEqual(odd["characterisation"], "giuga-window")
        self.assertLessEqual(odd["excess_integer"], 0)
        self.assertEqual(odd["slack"], 1)          # omega = 2 for 15
        self.assertEqual(sign_regime(105)["slack"], 2)   # omega = 3
        self.assertEqual(sign_regime(3)["slack"], 0)     # omega = 1

    def test_the_block_saves_and_there_is_nothing_left_to_prove(self):
        self.assertEqual(sign_regime(60)["regime"], "block-saves")
        self.assertEqual(sign_regime(180)["regime"], "block-saves")
        self.assertEqual(sign_regime(180)["block"], 3)


class TheCensusReportsItsOwnRefutation(unittest.TestCase):
    def test_a_refuted_lemma_comes_back_as_a_violation(self):
        with mock.patch("exact_lattice_quantiles.pricing.worst_antichain_weight",
                        return_value=(10**6, ())):
            census = sign_theorem_census(12)
        self.assertNotEqual(census["lemma_violations"], [])
        self.assertEqual(census["verdict"], "the lemma is false")

    def test_a_refuted_theorem_comes_back_as_a_violation(self):
        with mock.patch("exact_lattice_quantiles.pricing.antichain_saving_ceiling",
                        return_value=(5, (2, 3))):
            census = sign_theorem_census(12)
        self.assertNotEqual(census["sign_violations"], [])
        self.assertEqual(census["sign_verdict"],
                         "an antichain saves where the theorem forbids it")

    def test_the_cli_turns_a_violation_into_rc_one(self):
        with mock.patch("exact_lattice_quantiles.pricing.worst_antichain_weight",
                        return_value=(10**6, ())):
            with contextlib.redirect_stdout(io.StringIO()):
                rc = main(["sign-theorem", "--max-period", "12"])
        self.assertEqual(rc, 1)

    def test_no_band_member_can_have_two_adic_other_than_one(self):
        # v2(L) >= 2 puts no odd member in the block, so the block saves
        # L * (e1 - 1) / 2 >= 1 there and the band cannot reach it.  Pinned
        # over the whole census rather than argued once.
        census = sign_theorem_census(900)
        self.assertTrue(census["band_rows"])
        for row in census["band_rows"]:
            self.assertEqual(row["two_adic"], 1, row)
            self.assertEqual(row["excess"], "1", row)

    def test_a_period_below_two_is_a_bad_request(self):
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main(["sign-theorem", "--max-period", "1"]), 3)


def _primes(n: int) -> tuple[int, ...]:
    from exact_lattice_quantiles.arithmetic import primes_of
    return primes_of(n)


if __name__ == "__main__":
    unittest.main()
