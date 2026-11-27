"""Cofactor two is the one cofactor a theorem can decide without a grid.

The median test has two sides.  Asking whether ``a`` is the median of the grid
``2a`` asks that ``a`` clear the middle from below *and* from above, and only at
``q = 2`` does the upper side read ``a + X >= a`` for a mass ``X`` that is never
negative: the condition collapses to a single inequality.  That is why cofactor
two has a closed form and cofactor three does not, and it is the whole content
of :func:`carrier_two_law`.

Four things are pinned here.

* The law agrees with the definition.  ``carriers(2 * a)`` is a search; the law
  is a comparison on the divisor's own factorization; they name the same
  ``a`` for every ``a`` to the horizon.  The definition is the referee.
* The closed form is the same statement twice: for even ``a = 2^s m`` the law
  is exactly ``2 * U(m) >= 2 * m - phi(m)``, a half-mass comparison on the odd
  part.  A second reader that runs the same algebra is not a second reader, so
  the point of the test is that both routes are *the same route*, and the
  referee is elsewhere.
* Equality is a hit, not a near miss: ``18`` sits exactly on the tie for the
  grid ``36`` and the grid carries it.  The tie is not a corner case either --
  it holds for twenty divisors to ``600``, starting at ``6`` -- so a strict
  reading loses all twenty and that is the mutant this file keeps honest.
* The family ``6p`` never passes: ``U(3p) = phi(3p) = 2(p - 1)`` while the base
  demands ``3p - (p - 1)``, and the first is smaller for every prime ``p >= 5``.
  So every ``6p`` is a two-gap, and a two-gap is not yet the gap.

The law is a theorem about membership and **not** a speedup: measured, it
replaces ``2482`` of ``231 109`` grid constructions over the same sweep, and
building the factorization to evaluate it costs more than the grid it skips.
The claim is epistemic.  It is not repeated in code as a shortcut, and no test
here asserts a timing.
"""

from __future__ import annotations

import io
import json
import unittest
from contextlib import redirect_stdout

from exact_lattice_quantiles.arithmetic import (euler_phi, lower_half_mass,
                                               odd_part, upper_half_mass)
from exact_lattice_quantiles.cli import RC_OK, main
from exact_lattice_quantiles.lattice import carriers
from exact_lattice_quantiles.pricing import carrier_two_census, carrier_two_law

HORIZON = 600


class TestTheLawAgreesWithTheDefinition(unittest.TestCase):
    """``carriers(2 * a)`` is the referee; the law is what is on trial."""

    def test_the_law_names_the_same_divisors_as_the_grid(self):
        wrong = [a for a in range(2, HORIZON + 1)
                 if carrier_two_law(a) != (a in carriers(2 * a))]
        self.assertEqual(wrong, [])

    def test_the_even_closed_form_is_the_same_statement(self):
        for a in range(2, HORIZON + 1, 2):
            m = odd_part(a)
            self.assertEqual(carrier_two_law(a),
                             2 * upper_half_mass(m) >= 2 * m - euler_phi(m),
                             f"a = {a}, odd part {m}")

    def test_the_two_halves_of_a_totient_sum_are_the_number_itself(self):
        for n in range(1, HORIZON + 1):
            self.assertEqual(upper_half_mass(n) + lower_half_mass(n), n)


class TestEqualityIsAHit(unittest.TestCase):
    """The tie is the measurement, not an edge case to round away."""

    def test_eighteen_sits_on_the_tie_and_the_grid_carries_it(self):
        m = odd_part(18)
        self.assertEqual(m, 9)
        self.assertEqual(2 * upper_half_mass(m), 2 * m - euler_phi(m))
        self.assertTrue(18 in carriers(36))
        self.assertTrue(carrier_two_law(18))

    def test_the_strict_reading_loses_every_tie_and_there_are_twenty(self):
        strict = [a for a in range(2, HORIZON + 1, 2)
                  if carrier_two_law(a) != (
                      2 * upper_half_mass(odd_part(a))
                      > 2 * odd_part(a) - euler_phi(odd_part(a)))]
        self.assertIn(18, strict)
        self.assertNotIn(30, strict)
        for a in strict:
            m = odd_part(a)
            self.assertEqual(2 * upper_half_mass(m), 2 * m - euler_phi(m))
        self.assertEqual(len(strict), 20)


class TestTheFamilySixPIsNeverCarriedAtTwo(unittest.TestCase):
    """A whole family, decided by hand and then by the law."""

    def test_no_six_p_passes_cofactor_two(self):
        from exact_lattice_quantiles.arithmetic import primes_of

        primes = [p for p in primes_of(200) if 5 <= p <= 200]
        self.assertTrue(primes)
        passing = [6 * p for p in primes if carrier_two_law(6 * p)]
        self.assertEqual(passing, [])

    def test_and_the_odd_part_misses_the_bar_it_would_need(self):
        from exact_lattice_quantiles.arithmetic import primes_of

        for p in [q for q in primes_of(200) if 5 <= q <= 200]:
            self.assertEqual(odd_part(6 * p), 3 * p)
            self.assertEqual(upper_half_mass(3 * p), euler_phi(3 * p))


class TestTheCensusCountsTheDisagreement(unittest.TestCase):
    """What the law buys and what it does not: the gap is inside the two-gap."""

    def test_the_census_finds_no_disagreement(self):
        census = carrier_two_census(HORIZON)
        self.assertEqual(census["disagreements"], ())
        self.assertEqual(census["verdict"].count("no disagreement"), 1)

    def test_the_gap_is_a_subset_and_strictly_smaller(self):
        from exact_lattice_quantiles.pricing import carrier_window_census

        gap = set(carrier_window_census(HORIZON)["gap"])
        two_gap = set(carrier_two_census(HORIZON)["two_gap"])
        self.assertLess(gap, two_gap)
        self.assertIn(105, two_gap - gap)
        self.assertNotIn(105, gap)

    def test_every_divisor_coprime_to_six_is_settled_by_cofactor_two(self):
        census = carrier_two_census(HORIZON)
        self.assertEqual(census["coprime_to_six_decided"],
                         census["coprime_to_six"])


class TestTheCLIPublishesIt(unittest.TestCase):
    """The census reaches the printer, and the printer does not recompute it."""

    def run_cli(self, *argv):
        out, err = io.StringIO(), io.StringIO()
        from contextlib import redirect_stderr

        with redirect_stdout(out), redirect_stderr(err):
            rc = main(list(argv))
        return rc, out.getvalue(), err.getvalue()

    def test_cofactor_law_runs_and_reports_the_two_gap(self):
        rc, out, _ = self.run_cli("carrier-law", "--max-period", "60",
                                  "--max-divisor", "120", "--no-cross-check")
        self.assertEqual(rc, RC_OK)
        self.assertIn("cofactor two", out)
        self.assertIn("half-mass comparison", out)

    def test_the_json_payload_carries_the_census(self):
        rc, out, _ = self.run_cli("carrier-law", "--max-period", "60",
                                  "--max-divisor", "120", "--no-cross-check",
                                  "--json")
        self.assertEqual(rc, RC_OK)
        payload = json.loads(out)
        self.assertEqual(payload["cofactor_two"]["disagreements"], [])
        self.assertEqual(payload["cofactor_two"]["max_divisor"], 120)


if __name__ == "__main__":
    unittest.main()
