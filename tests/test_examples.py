"""The examples are claims: if one stops printing its numbers, the README lies."""

from __future__ import annotations

import os
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE = ROOT / "examples" / "clinic_period.py"


def run_example() -> subprocess.CompletedProcess:
    env = dict(os.environ, PYTHONPATH="src")
    return subprocess.run([sys.executable, str(EXAMPLE)], cwd=ROOT, env=env,
                          capture_output=True, text=True)


class ClinicPeriodTest(unittest.TestCase):
    def test_it_runs_clean(self):
        result = run_example()
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_it_prints_the_numbers_the_readme_quotes(self):
        out = run_example().stdout
        self.assertIn("period               180 minutes", out)
        self.assertIn("instance ceiling     1/3 = 60 minute(s)", out)
        self.assertIn("cheapest design      47 points: cadences [12, 18, 20, 45]", out)
        self.assertIn("    64 points -> margin 1/3 = 60 minute(s)", out)
        self.assertIn("still keeps        17/180 = 17 minute(s)", out)

    def test_the_service_level_stays_inside_the_ceiling(self):
        out = run_example().stdout
        self.assertIn("ceiling check      keep + drift = 1/20 <= 1/3: True", out)


if __name__ == "__main__":
    unittest.main()
