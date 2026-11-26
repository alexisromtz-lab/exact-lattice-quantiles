"""The lemma read as an extremum, with the warrant of each region attached.

The package published two readings of the lemma the closed form rests on.  The
**ranking** reading -- the drop-one block saves most -- is false from 176 primes
on (:mod:`test_strong_lemma_frontier`).  The **sign** reading -- no antichain
saves where no block does -- is proved while the reciprocal prime sum stays at
most one.  This file pins a third sentence that neither of them states: the
block's weight ``sum(period / p)`` is not only a bound every antichain respects,
it is the **maximum**, and it is attained by the block itself.

Two lines make that a theorem in the proved region: the block is a candidate, so
the maximum is at least its weight; the lemma caps every antichain at that same
weight; therefore equality.  In the band the upper bound is unproved, so there
the equality is a measurement -- and the measurement is tight, because a saving
needs ``weight >= period + slack + 2`` and the maximum sits two points below it.

The last class keeps the honest edge on the claim: the equality is *not* a
theorem about all periods, and the same primorial that refutes the ranking
reading refutes it, in exact integers, inside this file.
"""

from __future__ import annotations

import contextlib
import io
import json
import unittest
from itertools import combinations
from math import lcm, prod

from exact_lattice_quantiles.arithmetic import primes_of
from exact_lattice_quantiles.cli import main
from exact_lattice_quantiles.pricing import (block_saving, extremal_census,
                                             heaviest_antichain, sign_regime,
                                             strong_lemma_frontier)

HORIZON = 2000
BAND = (30, 858, 1722)


def _lawful(limit: int = HORIZON) -> list[int]:
    return [n for n in range(2, limit + 1) if len(primes_of(n)) >= 2]


class TheBlockIsTheMaximum(unittest.TestCase):
    def test_the_bound_is_attained_on_every_period_never_just_approached(self):
        # The drop-one block is always a candidate, so a maximum below this
        # number would mean the search missed a candidate it cannot miss.
        for n in _lawful(600):
            row = heaviest_antichain(n)
            self.assertEqual(row["weight"], row["bound"], n)
            self.assertEqual(row["gap"], 0, n)

    def test_the_witness_is_the_drop_one_block_and_not_another_tie(self):
        for n in _lawful(600):
            row = heaviest_antichain(n)
            self.assertTrue(row["block_attains"], n)
            self.assertEqual(tuple(row["witness"]), tuple(row["drop_one"]), n)

    def test_the_witness_really_is_an_antichain_that_spans_the_period(self):
        for n in _lawful(200):
            witness = heaviest_antichain(n)["witness"]
            self.assertEqual(lcm(*witness), n, n)
            for a, b in combinations(witness, 2):
                self.assertNotEqual(a % b, 0, (n, a, b))
                self.assertNotEqual(b % a, 0, (n, a, b))
            self.assertEqual(sum(witness), heaviest_antichain(n)["weight"], n)

    def test_the_proved_region_carries_a_theorem_and_not_a_census(self):
        for n in _lawful(HORIZON):
            row = heaviest_antichain(n)
            if row["regime"] == "proven-zero":
                self.assertEqual(row["warrant"], "theorem", n)
            else:
                # The lemma's hypothesis is what makes it a proof; without it
                # nothing here may be called proved.
                self.assertEqual(row["warrant"], "measured", n)

    def test_the_band_is_measured_and_its_miss_is_the_smallest_one_possible(self):
        for n in BAND:
            row = heaviest_antichain(n)
            self.assertEqual(row["regime"], "band", n)
            self.assertEqual(row["warrant"], "measured", n)
            self.assertEqual(row["gap"], 0, n)
            self.assertEqual(row["miss"], 2, n)
            self.assertTrue(row["tightest_miss"], n)
            self.assertEqual(row["weight"], n + row["slack"], n)

    def test_the_miss_is_measured_against_the_first_weight_that_saves(self):
        for n in BAND:
            row = heaviest_antichain(n)
            first_saving = n + row["slack"] + 2
            self.assertEqual(first_saving - row["weight"], row["miss"], n)
            self.assertEqual(row["miss"], 2, n)

    def test_the_miss_is_the_distance_from_the_top_of_the_window(self):
        # miss = slack - m + 2, so tightness is exactly m == slack.  The band is
        # 1 <= m <= slack, and only its top sits two points away.
        for n in _lawful(600):
            row = heaviest_antichain(n)
            self.assertEqual(row["bound"], n + row["excess"], n)
            self.assertEqual(row["miss"], row["slack"] - row["excess"] + 2, n)
            self.assertEqual(row["tightest_miss"],
                             row["excess"] == row["slack"], n)

    def test_one_point_of_weight_is_never_the_shortfall(self):
        # weight and period + slack share a parity, so a miss is always even:
        # "one weight point short" is not a state an antichain can be in.
        for n in _lawful(600):
            row = heaviest_antichain(n)
            self.assertEqual(row["weight"] % 2, (n + row["slack"]) % 2, n)
            self.assertEqual(row["miss"] % 2, 0, n)

    def test_a_prime_period_has_no_antichain_to_weigh(self):
        for n in (1, 2, 4, 7, 49):
            self.assertLess(len(primes_of(n)), 2, n)
            with self.assertRaises(ValueError):
                heaviest_antichain(n)


class TheCensusByRegion(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.census = extremal_census(HORIZON)

    def test_no_antichain_outweighs_the_block_on_this_population(self):
        self.assertEqual(self.census["violations"], [])
        self.assertTrue(self.census["block_attains_everywhere"])

    def test_the_regions_partition_the_population(self):
        self.assertEqual(sum(row["rows"] for row in self.census["regions"].values()),
                         self.census["periods"])
        self.assertEqual(self.census["periods"], len(_lawful()))

    def test_every_region_reports_a_zero_gap_and_its_own_warrant(self):
        regions = self.census["regions"]
        self.assertEqual(sorted(regions),
                         ["band", "block-saves", "proven-zero"])
        for name, row in regions.items():
            self.assertEqual(row["gap_zero"], row["rows"], name)
            self.assertEqual(row["block_attains"], row["rows"], name)
        self.assertEqual(regions["proven-zero"]["warrant"], "theorem")
        self.assertEqual(regions["block-saves"]["warrant"], "measured")
        self.assertEqual(regions["band"]["warrant"], "measured")

    def test_the_band_it_reports_is_the_one_the_regimes_name(self):
        self.assertEqual([row["period"] for row in self.census["band_rows"]],
                         [n for n in _lawful() if sign_regime(n)["regime"] == "band"])

    def test_the_tightest_period_is_the_band_member_that_misses_by_one_unit(self):
        tightest = self.census["tightest"]
        self.assertEqual(tightest["period"], min(BAND))
        self.assertEqual(tightest["miss"], 2)
        self.assertEqual(self.census["regions"]["band"]["tightest_miss"], 3)

    def test_the_band_members_are_tight_because_they_sit_at_the_top_of_the_window(self):
        # The tightness condition and the odd band's search threshold are the
        # same predicate m == slack, so the two instruments agree by construction
        # rather than by both having run on the same periods.
        for row in self.census["band_rows"]:
            self.assertEqual(row["excess"], row["slack"], row["period"])
            self.assertTrue(row["tightest_miss"], row["period"])

    def test_a_horizon_below_the_first_member_is_refused(self):
        with self.assertRaises(ValueError):
            extremal_census(5)


class TheEqualityIsNotATheoremAboutAllPeriods(unittest.TestCase):
    """The primorial that kills the ranking reading kills this claim too.

    Both the block ``{n / p}`` and the pair family ``{n / (p q)}`` are
    antichains with least common multiple ``n``, so their weights are compared
    directly: ``sum(n / (p q))`` against ``sum(n / p)``.  Dividing by ``n``
    leaves the two elementary symmetric functions of the reciprocals, and the
    second overtakes the first at 176 primes.  The refutation is exact integer
    arithmetic, and it lands **inside** the region where the block already
    saves -- not in the band, where the sign reading needs the bound.
    """

    def test_at_176_primes_the_pairs_weigh_more_than_the_block(self):
        primes = _first_primes(176)
        self.assertEqual(max(primes), 1049)
        n = prod(primes)
        block = sum(n // p for p in primes)
        pairs = sum(n // (a * b) for a, b in combinations(primes, 2))
        self.assertGreater(pairs, block)
        # The same direction under the integer rule the package prices by, so
        # the failure is not an artefact of dropping the floors.
        charged = sum((n // (a * b)) // 2 for a, b in combinations(primes, 2))
        self.assertGreater(charged, sum((n // p) // 2 for p in primes))

    def test_at_175_primes_the_block_is_still_the_heaviest(self):
        primes = _first_primes(175)
        n = prod(primes)
        block = sum(n // p for p in primes)
        pairs = sum(n // (a * b) for a, b in combinations(primes, 2))
        self.assertLess(pairs, block)

    def test_so_the_census_publishes_its_frontier_and_not_a_universal(self):
        census = extremal_census(300)
        self.assertIn("176 primes", census["frontier"])
        self.assertIn("not an upper bound in general", census["frontier"])
        self.assertIn("sum(1/p) <= 1", census["theorem"])

    def test_the_region_that_needs_the_bound_is_the_one_measured_true(self):
        # The bound is only *used* where the block saves nothing: elsewhere the
        # block's own witness already decides the sign reading, so the failure
        # above lands outside the statement the closed form depends on.
        for n in BAND:
            self.assertEqual(block_saving(n), 0, n)
        self.assertGreater(strong_lemma_frontier(200)["crossing"]["block_saving"], 0)
        bounded = [n for n in _lawful() if block_saving(n) == 0]
        self.assertTrue(all(sign_regime(n)["regime"] != "block-saves"
                            for n in bounded))


class TheCommandLine(unittest.TestCase):
    def run_cli(self, *argv: str) -> tuple[int, str]:
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(list(argv))
        return rc, out.getvalue() + err.getvalue()

    def test_the_census_answers_zero_and_names_its_regions(self):
        rc, out = self.run_cli("extremal-gap", "--limit", "300")
        self.assertEqual(rc, 0)
        for name in ("proven-zero", "block-saves", "band"):
            self.assertIn(name, out)
        self.assertIn("30", out)

    def test_a_horizon_that_never_reaches_the_band_is_not_a_green(self):
        rc, out = self.run_cli("extremal-gap", "--limit", "20")
        self.assertEqual(rc, 4)
        self.assertIn("NO_MEDIDO", out)

    def test_the_json_carries_the_warrant_of_each_region(self):
        rc, out = self.run_cli("--json", "extremal-gap", "--limit", "300")
        self.assertEqual(rc, 0)
        data = json.loads(out)
        self.assertEqual(data["regions"]["proven-zero"]["warrant"], "theorem")
        self.assertEqual(data["regions"]["band"]["warrant"], "measured")
        self.assertEqual(data["violations"], [])


def _first_primes(count: int) -> list[int]:
    out: list[int] = []
    candidate = 2
    while len(out) < count:
        if all(candidate % p for p in out):
            out.append(candidate)
        candidate += 1
    return out


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
