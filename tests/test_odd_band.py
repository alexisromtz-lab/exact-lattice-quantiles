"""The odd half of the band: a window through Giuga's equation, walked empty.

An odd period has ``omega`` odd entries in its block, so the slack the block may
give away is ``omega - 1`` and not ``1``.  Three claims are pinned here: the
unified law that turns the three regimes into one integer predicate, the
minimal-feasibility theorem that lets the search prune, and the search itself --
with a control that shows the walk is not blind.
"""

from __future__ import annotations

import contextlib
import io
import unittest
from fractions import Fraction
from itertools import combinations
from math import gcd

from exact_lattice_quantiles.arithmetic import primes_of
from exact_lattice_quantiles.cli import RC_ABOVE_CEILING, RC_BAD_REQUEST
from exact_lattice_quantiles.cli import RC_NO_MEDIDO, RC_OK, main
from exact_lattice_quantiles.pricing import (band_slack, band_window,
                                            block_saving, odd_band_census,
                                            odd_members, prime_excess,
                                            reciprocal_prime_sum, sign_regime)

FLOOR = (3, 5, 7, 11, 13, 17, 19, 23, 29)


def _product(primes) -> int:
    out = 1
    for prime in primes:
        out *= prime
    return out


def _excess(primes) -> int:
    """``sum(L / p) - L`` straight from the primes, without factoring a product."""
    product = _product(primes)
    return sum(product // p for p in primes) - product


def _odd_primes(limit: int) -> tuple[int, ...]:
    return tuple(p for p in range(3, limit + 1, 2)
                 if all(p % d for d in range(3, int(p ** 0.5) + 1, 2)))


class TheUnifiedLaw(unittest.TestCase):
    def test_the_block_saving_is_the_excess_minus_the_slack_halved(self):
        checked = 0
        for n in range(6, 1201):
            if len(primes_of(n)) < 2:
                continue
            self.assertEqual(block_saving(n),
                             (prime_excess(n) - band_slack(n)) // 2, f"n = {n}")
            checked += 1
        self.assertGreater(checked, 900)

    def test_the_slack_is_what_the_parity_of_the_period_decides(self):
        self.assertEqual(band_slack(858), 1)   # exactly one factor 2
        self.assertEqual(band_slack(60), 0)    # four divides it
        self.assertEqual(band_slack(45), 1)    # odd: omega - 1 with omega = 2
        self.assertEqual(band_slack(3 * 5 * 7), 2)
        for n in (45, 105, 3 * 5 * 7 * 11, 5 * 7 * 11 * 13):
            self.assertEqual(band_slack(n), len(primes_of(n)) - 1)
        for n in (9, 25, 45, 105, 3 * 5 * 7 * 11):
            self.assertEqual(band_slack(n), odd_members(n) - 1)

    def test_the_library_excess_is_the_formula_it_claims(self):
        for chosen in ((3, 5), (3, 5, 7), (3, 5, 7, 11), FLOOR):
            self.assertEqual(prime_excess(_product(chosen)), _excess(chosen))

    def test_the_excess_recurses_by_one_prime_at_a_time(self):
        # M(P U {q}) = q * M(P) + L(P), the identity the theorem rests on.
        primes = _odd_primes(60)
        for size in range(1, 5):
            for chosen in combinations(primes, size):
                product = _product(chosen)
                excess = _excess(chosen)
                for prime in (p for p in primes if p > chosen[-1]):
                    self.assertEqual(_excess((*chosen, prime)),
                                     prime * excess + product,
                                     f"{chosen} + {prime}")

    def test_the_regime_is_the_predicate_and_not_an_extra_case(self):
        for period in range(2, 601):
            state = sign_regime(period)
            self.assertEqual(state["in_band"], state["regime"] == "band",
                             f"period {period}")
            self.assertTrue(state["characterised"], f"period {period}")
            self.assertEqual(state["in_band"],
                             1 <= prime_excess(period) <= band_slack(period))


class TheWindow(unittest.TestCase):
    def test_the_band_above_a_prefix_is_one_interval_of_primes(self):
        for chosen in ((3, 5), (3, 7), (3, 5, 7), (3, 5, 7, 11), (5, 7, 11)):
            product = _product(chosen)
            excess = _excess(chosen)
            threshold = len(chosen)
            low, high = band_window(product, excess, threshold)
            for prime in _odd_primes(200):
                member = prime * excess + product
                self.assertEqual(1 <= member <= threshold, low <= prime <= high,
                                 f"{chosen} + {prime}")

    def test_the_window_is_never_wider_than_the_slack(self):
        for chosen in ((3, 5), (3, 5, 7), (3, 5, 7, 11), (3, 5, 7, 11, 13)):
            product = _product(chosen)
            excess = _excess(chosen)
            for threshold in (1, 2, 5):
                low, high = band_window(product, excess, threshold)
                self.assertLessEqual(high - low + 1, threshold)

    def test_the_window_above_three_and_five_closes_on_two(self):
        # 2 * 3 * 5 = 30 is the odd shape's corner: excess exactly 1.
        self.assertEqual(band_window(15, prime_excess(15), 2), (2, 2))
        self.assertEqual(prime_excess(30), 1)

    def test_a_non_negative_excess_is_not_a_prefix(self):
        with self.assertRaises(ValueError):
            band_window(30, 1, 2)
        with self.assertRaises(ValueError):
            band_window(15, -7, 0)


class TheFeasibilityFloor(unittest.TestCase):
    def test_eight_odd_primes_cannot_reach_one_and_nine_can(self):
        self.assertLess(reciprocal_prime_sum(_product(FLOOR[:8])), 1)
        self.assertGreater(reciprocal_prime_sum(_product(FLOOR)), 1)

    def test_the_first_witness_is_far_above_any_window(self):
        excess = prime_excess(_product(FLOOR))
        self.assertEqual(excess, 108169298)
        self.assertEqual(gcd(excess, _product(FLOOR)), 1)
        self.assertGreater(excess, len(FLOOR) - 1)

    def test_a_horizon_that_cannot_cross_one_is_not_measured(self):
        census = odd_band_census(8, 200)
        self.assertEqual(census["status"], "below-floor")
        self.assertEqual(census["candidates"], [])
        self.assertIn("NO_MEDIDO", census["verdict"])

    def test_a_cap_below_the_ninth_prime_has_no_floor_either(self):
        census = odd_band_census(12, 23)  # 23 is the eighth odd prime
        self.assertEqual(census["status"], "below-floor")
        self.assertIsNone(census["floor"]["omega"])


class TheWalkIsNotBlind(unittest.TestCase):
    """The control: lift the threshold and the same walk must find members.

    Without it, an empty answer would be consistent with a search that prunes
    everything, which is the failure mode this package exists to avoid.
    """

    def test_exhaustive_enumeration_agrees_that_the_band_is_empty(self):
        for limit, omega in ((60, 5), (90, 7)):
            primes = _odd_primes(limit)
            for size in range(2, omega + 1):
                for chosen in combinations(primes, size):
                    excess = _excess(chosen)
                    self.assertFalse(1 <= excess <= size - 1, f"{chosen}")

    def test_lifting_the_threshold_finds_the_floor_witness(self):
        census = odd_band_census(9, 30, slack=10**9)
        self.assertEqual(len(census["candidates"]), 1)
        member = census["candidates"][0]
        self.assertEqual(tuple(member["primes"]), FLOOR)
        self.assertEqual(member["excess"], 108169298)
        self.assertGreater(member["excess"], member["slack"])

    def test_no_pruned_prefix_could_have_been_a_member(self):
        # The pruning rule, falsified directly: a set with an excess above one
        # either leaves the window itself or has a proper subset that does.
        primes = _odd_primes(60)
        for size in (9, 10):
            for chosen in combinations(primes, size):
                excess = _excess(chosen)
                if excess < 1:
                    continue
                if excess <= size - 1:
                    self.assertTrue(
                        any(_excess(sub) >= 1
                            for sub in combinations(chosen, size - 1)),
                        f"{chosen} is a member with a feasible proper subset")


class TheCensusOfTheOddBand(unittest.TestCase):
    def test_nothing_below_the_horizon_is_in_the_band(self):
        census = odd_band_census(10, 100)
        self.assertEqual(census["status"], "measured")
        self.assertEqual(census["candidates"], [])
        self.assertGreater(census["prefixes"], 1000)
        self.assertIn("no odd period in the band", census["verdict"])

    def test_the_boundary_is_named_and_not_hidden(self):
        census = odd_band_census(9, 30)
        self.assertIn("Giuga", census["boundary"])
        self.assertEqual(census["floor"]["omega"], 9)

    def test_a_bad_horizon_is_refused(self):
        with self.assertRaises(ValueError):
            odd_band_census(1, 100)
        with self.assertRaises(ValueError):
            odd_band_census(10, 2)


class TheCommandLine(unittest.TestCase):
    def _run(self, *argv):
        err = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(err):
            rc = main(list(argv))
        return rc, err.getvalue()

    def test_an_empty_horizon_exits_zero(self):
        self.assertEqual(self._run("odd-band", "--max-omega", "10",
                                   "--max-prime", "60")[0], RC_OK)

    def test_a_horizon_that_cannot_cross_one_exits_no_medido(self):
        self.assertEqual(self._run("odd-band", "--max-omega", "8")[0],
                         RC_NO_MEDIDO)

    def test_a_member_exits_one_because_it_refutes_the_package(self):
        self.assertEqual(self._run("odd-band", "--max-omega", "9",
                                   "--max-prime", "30", "--slack", "1000000000")[0],
                         RC_ABOVE_CEILING)

    def test_a_bad_horizon_exits_three(self):
        self.assertEqual(self._run("odd-band", "--max-omega", "1")[0],
                         RC_BAD_REQUEST)

    def test_json_carries_the_verdict(self):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            self.assertEqual(main(["--json", "odd-band", "--max-omega", "9",
                                   "--max-prime", "30"]), RC_OK)
        self.assertIn('"candidates": []', buffer.getvalue())
        self.assertIn('"floor"', buffer.getvalue())


class TheSignRegimeGainedACharacterisation(unittest.TestCase):
    def test_an_odd_period_is_now_a_window_instead_of_a_none(self):
        odd = sign_regime(45)
        self.assertEqual(odd["characterisation"], "giuga-window")
        self.assertTrue(odd["characterised"])
        self.assertEqual(odd["slack"], 1)
        self.assertEqual(odd["excess_integer"], prime_excess(45))

    def test_the_band_members_still_carry_the_giuga_equation(self):
        for period in (30, 858, 1722):
            state = sign_regime(period)
            self.assertEqual(state["characterisation"], "giuga-equation")
            self.assertEqual(state["excess_integer"], 1)
            self.assertEqual(state["slack"], 1)

    def test_four_dividing_the_period_leaves_no_member_possible(self):
        self.assertEqual(sign_regime(36)["characterisation"],
                         "no-member-possible")
        for period in range(4, 601, 4):
            self.assertNotEqual(sign_regime(period)["regime"], "band",
                                f"period {period}")


if __name__ == "__main__":
    unittest.main()
