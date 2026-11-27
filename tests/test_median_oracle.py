"""A third route to ``U(T)``, built to share nothing with the other two.

The package has read its universe twice since v0.21.0 and said, correctly, that
the two readings are not independent: both end at :func:`carriers`, so a bias
in ``mass`` or in the divisor list would survive the double reading and present
itself as agreement.  The third route is the one that makes the claim testable:
it sorts the class values of a grid and takes the middle, and it builds its
grids by division, so the only step it shares with the others is the definition
of a carrier.

Two consequences are pinned here, and the second one is the interesting one.

* The median needs no totient, so the frequency of each divisor among the
  classes has to come out as ``mass(d) / 2``.  That identity is what both other
  routes *assume*; here it is derived, on every divisor of every even grid.
* The tie branch in :func:`carriers` is load-bearing.  A reading that takes a
  single middle value is the natural wrong reading of this route, and it parts
  from the mass at grid ``4`` -- not at some far horizon, but at the second
  grid of the period ``4``.  The mutant is exercised through the census so the
  falsifier is the instrument and not a comment.
"""

from __future__ import annotations

import contextlib
import io
import unittest
from collections import Counter
from unittest import mock

from exact_lattice_quantiles.arithmetic import divisors, mass
from exact_lattice_quantiles.cli import RC_ABOVE_CEILING, main
from exact_lattice_quantiles.lattice import (carriers, class_values, grids,
                                             grids_by_trial, median_carriers,
                                             median_oracle_census, universe,
                                             universe_by_median)


def _single_middle(grid: int) -> tuple[int, ...]:
    """The median route read naively: one middle value, tie branch dropped."""
    values = sorted(class_values(grid))
    return (values[len(values) // 2],)


CENSUS = median_oracle_census(1200)


class TheCarrierIsTheMedianOfTheClasses(unittest.TestCase):
    def test_a_class_meets_the_divisor_the_grid_leaves_of_it(self):
        self.assertEqual(class_values(12), (12, 6, 4, 3, 12, 2))
        self.assertTrue(all(12 % d == 0 for d in class_values(12)))

    def test_the_frequency_of_a_divisor_is_half_its_mass(self):
        """The identity both other routes assume, derived here instead."""
        failures = []
        for grid in range(2, 361, 2):
            counts = Counter(class_values(grid))
            for d in divisors(grid)[1:]:
                if counts.get(d, 0) != mass(d) // 2:
                    failures.append((grid, d, counts.get(d, 0), mass(d) // 2))
        self.assertEqual(failures, [])

    def test_four_dividing_the_grid_puts_the_median_between_two_carriers(self):
        for grid in (4, 8, 12, 36):
            self.assertEqual(median_carriers(grid), carriers(grid))
            self.assertEqual(len(median_carriers(grid)), 2, grid)

    def test_and_a_grid_without_that_middle_has_one(self):
        for grid in (2, 6, 180, 360, 1080):
            self.assertEqual(median_carriers(grid), carriers(grid))
            self.assertEqual(len(median_carriers(grid)), 1, grid)

    def test_an_odd_grid_has_no_middle_and_is_refused_by_name(self):
        with self.assertRaises(ValueError) as caught:
            median_carriers(15)
        self.assertIn("even grid", str(caught.exception))
        self.assertEqual(carriers(15), (15,))

    def test_the_grid_table_built_by_division_is_the_one_built_from_primes(self):
        for period in (2, 60, 180, 720, 1200):
            self.assertEqual(grids_by_trial(period), grids(period))


class TheThirdRouteAgreesAndTheCensusSaysSo(unittest.TestCase):
    def test_no_carrier_splits_and_no_universe_splits_to_1200(self):
        self.assertEqual(CENSUS["grid_splits"], ())
        self.assertEqual(CENSUS["carrier_splits"], ())
        self.assertEqual(CENSUS["universe_splits"], ())
        self.assertEqual(CENSUS["periods_measured"], 600)
        self.assertEqual(CENSUS["grids_measured"], 600)

    def test_the_census_reaches_every_period_and_not_a_sample_of_them(self):
        self.assertEqual(universe_by_median(180), universe(180))
        self.assertEqual(universe_by_median(1080), universe(1080))
        self.assertEqual(CENSUS["largest_universe"], max(
            len(universe(period)) for period in range(2, 1201, 2)))

    def test_a_reading_without_the_tie_branch_splits_at_the_first_grid(self):
        """The falsifier: drop the tie and the census must go red immediately."""
        with mock.patch("exact_lattice_quantiles.lattice.median_carriers",
                        _single_middle):
            split = median_oracle_census(12)
        self.assertEqual([row[0] for row in split["carrier_splits"]][:1], [4])
        self.assertEqual(split["carrier_splits"][0][1:], ((4,), (2, 4)))

    def test_a_union_absorbs_a_wrong_carrier_where_the_grid_does_not(self):
        """``U(T)`` is a union, so it hides what the grid shows.

        Dropping the tie loses carrier ``2`` of grid ``4`` and ``4`` of grid
        ``8``, and both are still carried by other grids of their periods: no
        universe below ``36`` changes, while the grids part at ``4``.  The
        universe is the claim and the grid is the instrument.
        """
        with mock.patch("exact_lattice_quantiles.lattice.median_carriers",
                        _single_middle):
            split = median_oracle_census(1200)
        self.assertEqual(split["universe_splits"][0][0], 36)
        self.assertNotIn(4, {row[0] for row in split["universe_splits"]})
        self.assertEqual(split["carrier_splits"][0][0], 4)

    def test_and_the_command_exits_one_on_that_mutant(self):
        out, err = io.StringIO(), io.StringIO()
        with mock.patch("exact_lattice_quantiles.lattice.median_carriers",
                        _single_middle):
            with contextlib.redirect_stdout(out), \
                    contextlib.redirect_stderr(err):
                rc = main(["median-oracle", "--max-period", "120"])
        self.assertEqual(rc, RC_ABOVE_CEILING)
        self.assertIn("SPLIT: grid 4", err.getvalue())

    def test_the_clean_run_names_all_three_readings_and_exits_zero(self):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = main(["median-oracle", "--max-period", "240"])
        self.assertEqual((rc, err.getvalue()), (0, ""))
        printed = out.getvalue()
        for label in ("grid splits", "carrier splits", "universe splits"):
            self.assertIn(label, printed)
        self.assertIn("a third time", printed)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
