"""The CLI is the contract in a shell: same numbers, exit codes that say why."""

from __future__ import annotations

import contextlib
import io
import json
import unittest

from exact_lattice_quantiles.cli import main


def run(*argv: str) -> tuple[int, str]:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = main(list(argv))
    return rc, buf.getvalue()


class TestExitCodes(unittest.TestCase):
    def test_ceiling(self):
        rc, out = run("ceiling", "180")
        self.assertEqual(rc, 0)
        self.assertIn("1/3", out)

    def test_json_ceiling_is_json(self):
        rc, out = run("--json", "ceiling", "180")
        self.assertEqual(rc, 0)
        data = json.loads(out)
        self.assertEqual(data["margin"], "1/3")
        self.assertEqual(data["critical_divisor"], 6)
        self.assertEqual(data["by_product"], "1/3")
        self.assertEqual(data["by_scan"], "1/3")

    def test_above_the_ceiling_exits_one(self):
        rc, out = run("--json", "solve", "180", "--max-cost", "1000",
                      "--probes", "8", "--min-margin", "2/5")
        self.assertEqual(rc, 1)
        self.assertEqual(json.loads(out)["verdict"], "INFEASIBLE_ABOVE_INSTANCE_CEILING")

    def test_over_budget_exits_two(self):
        rc, out = run("--json", "solve", "180", "--max-cost", "46",
                      "--probes", "8", "--min-margin", "1/6")
        self.assertEqual(rc, 2)
        payload = json.loads(out)
        self.assertEqual(payload["verdict"], "OVER_BUDGET")
        self.assertEqual(payload["cheapest"], 47)

    def test_the_frontier_prints_the_crossing(self):
        rc, out = run("strong-lemma", "--max-primes", "200")
        self.assertEqual(rc, 0)
        self.assertIn("crossing at 176 primes", out)
        self.assertIn("1049", out)
        # The savings have hundreds of digits: the headline must not dump them.
        self.assertLess(max(len(token) for token in out.split()), 40)

    def test_a_short_search_is_not_a_verdict(self):
        rc, out = run("strong-lemma", "--max-primes", "100")
        self.assertEqual(rc, 1)
        self.assertIn("NOT_MEASURED", out)

    def test_json_frontier_is_json_with_exact_integers(self):
        rc, out = run("--json", "strong-lemma", "--max-primes", "200")
        self.assertEqual(rc, 0)
        data = json.loads(out)
        self.assertEqual(data["crossing_primes"], 176)
        self.assertEqual(data["crossing"]["digits"], 440)
        self.assertEqual(data["crossing"]["excess"],
                         data["crossing"]["pair_saving"]
                         - data["crossing"]["block_saving"])

    def test_no_command_is_a_bad_request(self):
        rc, _ = run()
        self.assertEqual(rc, 3)

    def test_a_bad_rational_is_a_bad_request(self):
        rc, _ = run("solve", "180", "--max-cost", "64", "--probes", "4",
                    "--min-margin", "1/0")
        self.assertEqual(rc, 3)


class TestSolveOutput(unittest.TestCase):
    def test_solve_prints_the_design_and_the_verdict(self):
        rc, out = run("solve", "180", "--max-cost", "64", "--probes", "4",
                      "--min-margin", "1/3")
        self.assertEqual(rc, 0)
        self.assertIn("1/3", out)
        self.assertIn("OK", out)

    def test_designs_census_is_reproducible(self):
        rc, out = run("--json", "designs", "180")
        self.assertEqual(rc, 0)
        payload = json.loads(out)
        self.assertEqual(payload["designs"], 131072)
        self.assertEqual(payload["covering"], 96768)
        self.assertEqual(payload["cheapest"], 47)

    def test_drift_prints_the_knife_edge(self):
        # The witness margin and its vanishing tolerance in the same block: the
        # pair of numbers is the point of the command.
        rc, out = run("drift", "180", "--max-cost", "64", "--probes", "4")
        self.assertEqual(rc, 0)
        self.assertIn("witness margin   1/15", out)
        self.assertIn("coverage drift   1/15", out)
        self.assertIn("1/15            0", out)

    def test_json_drift_is_json(self):
        rc, out = run("--json", "drift", "180", "--max-cost", "64", "--probes", "4")
        self.assertEqual(rc, 0)
        data = json.loads(out)
        self.assertEqual(data["witness_margin"], "1/15")
        self.assertEqual(data["coverage_drift"], "1/15")
        self.assertEqual(data["kept_margin_to_drift"][0], ["1/15", "0"])
        self.assertEqual(data["kept_margin_to_drift"][-1], ["0", "1/15"])

    def test_census_and_ladder_run(self):
        rc, out = run("census", "2", "400")
        self.assertEqual(rc, 0)
        self.assertIn("167", out)
        rc, out = run("ladder", "180")
        self.assertEqual(rc, 0)
        self.assertIn("1/6", out)


class TestSelfTest(unittest.TestCase):
    def test_self_test_passes(self):
        rc, out = run("--self-test")
        self.assertEqual(rc, 0, out)
        self.assertNotIn("FALLA", out)


if __name__ == "__main__":
    unittest.main()
