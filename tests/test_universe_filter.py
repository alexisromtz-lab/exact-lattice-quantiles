"""The sign lemma transfers to any universe -- and the census must show it does.

The lemma is proved from ``lcm(A)`` alone, so the universe never enters the
argument: if ``sum(1/p) <= 1``, no antichain saves over *any* set of grids.  A
theorem that transfers is worth exactly the row that exercises it, because the
failure mode here is a proof that is cited where it does not reach.  So this
file pins the two columns separately:

* below one, both universes are clean -- the proof reaching the census;
* above one the proof is silent, and the universe is what decides, which is why
  ``60`` breaks over every divisor and not over the carriers.
"""

import io
import json
import unittest
from contextlib import redirect_stdout

from exact_lattice_quantiles.cli import main
from exact_lattice_quantiles.lattice import universe
from exact_lattice_quantiles.pricing import (carrier_block_of, divisors,
                                             frontier_by_blocks, prime_excess,
                                             reciprocal_filter_census)

CENSUS = reciprocal_filter_census(400)


def _row(period):
    for row in CENSUS["rows"]:
        if row["period"] == period:
            return row
    raise AssertionError(f"period {period} is not in the census")


class TheProofReachesTheCensus(unittest.TestCase):
    def test_no_period_below_one_breaks_in_either_universe(self):
        breaking = [row for row in CENSUS["rows"]
                    if not row["above_one"]
                    and (row["reduced"] > 0 or row["divisors"] > 0)]
        self.assertEqual(breaking, [])
        self.assertEqual(CENSUS["proved_without_break"], CENSUS["proved"])

    def test_the_count_of_proved_periods_is_not_an_assumption(self):
        # 187 is the number of even periods up to 400 whose reciprocal sum is
        # at most one; if the census stopped measuring the column, this fails.
        recomputed = [row for row in CENSUS["rows"] if prime_excess(row["period"]) <= 0]
        self.assertEqual(CENSUS["proved"], len(recomputed))
        self.assertEqual(CENSUS["proved"], 187)
        self.assertEqual(CENSUS["periods"], 200)


class TheFilterIsWhatDecidesAboveOne(unittest.TestCase):
    def test_sixty_breaks_over_every_divisor_and_not_over_the_carriers(self):
        row = _row(60)
        self.assertTrue(row["above_one"])
        self.assertEqual(row["reduced"], 0)
        self.assertEqual(row["reduced_witness"], ())
        self.assertEqual(row["divisors"], 1)
        self.assertEqual(row["divisors_witness"], (12, 20, 30))

    def test_thirty_breaks_in_neither_universe_because_the_band_absorbs_it(self):
        row = _row(30)
        self.assertTrue(row["above_one"])
        self.assertEqual((row["reduced"], row["divisors"]), (0, 0))

    def test_every_period_above_one_is_clean_over_the_carriers(self):
        band = [row for row in CENSUS["rows"] if row["above_one"]]
        self.assertEqual(CENSUS["band"], 13)
        self.assertEqual([row for row in band if row["reduced"] > 0], [])
        self.assertEqual(CENSUS["band_clean_reduced"], CENSUS["band"])
        self.assertEqual(CENSUS["band_breaking_divisors"], 12)


class TheColumnsAreBlockSavings(unittest.TestCase):
    def test_the_divisor_column_is_the_carrier_block_and_not_the_form(self):
        # The closed form's 9 984 discrepancies over D(180) are a different
        # object; the census counts savings of blocks, so the two must not be
        # read as the same column.
        kept, whole = carrier_block_of(60, set(divisors(60)))
        self.assertEqual(sum(d // 2 for d in kept) - whole // 2, _row(60)["divisors"])
        # The carriers are a strict subset here, which is the whole point.
        self.assertLess(len(universe(60)), len(divisors(60)))

    def test_the_census_agrees_with_the_frontier_it_is_a_row_of(self):
        reduced = frontier_by_blocks(400, reduced=True)
        divisors_only = frontier_by_blocks(400, reduced=False)
        self.assertIsNone(reduced["period"])
        self.assertEqual(divisors_only["period"], 60)
        first_breaking = [row["period"] for row in CENSUS["rows"]
                          if row["divisors"] > 0]
        self.assertEqual(first_breaking[0], divisors_only["period"])


class TheCommandLine(unittest.TestCase):
    def _run(self, *argv):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            code = main(list(argv))
        return code, buffer.getvalue()

    def test_the_census_runs_and_exits_zero(self):
        code, out = self._run("universes")
        self.assertEqual(code, 0)
        self.assertIn("sum 1/p at most one", out)
        self.assertIn("the filter decides", out)

    def test_json_carries_the_two_columns_and_the_verdict(self):
        code, out = self._run("--json", "universes", "--max-period", "60")
        self.assertEqual(code, 0)
        census = json.loads(out)
        self.assertEqual(census["band_breaking_divisors"], 1)
        self.assertIn("verdict", census)
        row = [r for r in census["rows"] if r["period"] == 60][0]
        self.assertEqual((row["reduced"], row["divisors"]), (0, 1))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
