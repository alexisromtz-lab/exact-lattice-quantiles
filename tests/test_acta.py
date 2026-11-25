"""The acta boundary as a gate: every divergence names the field, not the file."""

from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
from unittest import mock

from exact_lattice_quantiles import acta
from exact_lattice_quantiles.acta import (ManifestError, render,
                                          shipped_manifest, verify_file,
                                          verify_manifest)
from exact_lattice_quantiles.cli import main


def run(*argv: str) -> tuple[int, str]:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = main(list(argv))
    return rc, buf.getvalue()


def load() -> dict:
    return json.loads(shipped_manifest().read_text(encoding="utf-8"))


def with_anchor(mutate) -> dict:
    data = load()
    mutate(data["anchors"])
    return data


class TestTheShippedManifest(unittest.TestCase):
    def test_every_anchor_reproduces(self):
        checks = verify_file(shipped_manifest())
        divergent = [str(c) for c in checks if not c.ok]
        self.assertEqual(divergent, [])
        self.assertGreater(len(checks), 0)

    def test_the_three_kinds_are_all_exercised(self):
        kinds = {c.kind for c in verify_file(shipped_manifest())}
        self.assertEqual(kinds, {"ceiling", "design", "frontier"})

    def test_the_acta_table_is_the_one_anchored(self):
        checks = verify_file(shipped_manifest())
        seen = {(c.expected, c.measured) for c in checks if c.field == "margin"}
        for row in ((Fraction(2, 7), Fraction(2, 7)),
                    (Fraction(10, 33), Fraction(10, 33)),
                    (Fraction(4, 15), Fraction(4, 15))):
            self.assertIn(row, seen)


class TestADivergenceNamesTheField(unittest.TestCase):
    def test_a_wrong_margin_is_caught_at_the_ceiling_and_nowhere_else(self):
        data = with_anchor(lambda a: a[5].update({"margin": "5/16"}))
        divergent = [c for c in verify_manifest(data) if not c.ok]
        self.assertEqual(len(divergent), 1)
        self.assertEqual((divergent[0].kind, divergent[0].field),
                         ("ceiling", "margin"))
        self.assertEqual(divergent[0].expected, Fraction(5, 16))
        self.assertEqual(divergent[0].measured, Fraction(1, 3))
        self.assertIn("ceiling.margin", str(divergent[0]))

    def test_a_wrong_cost_names_cost_not_margin(self):
        data = with_anchor(lambda a: a[5].update({"cost_cover": 46}))
        divergent = [c for c in verify_manifest(data) if not c.ok]
        self.assertEqual([c.field for c in divergent], ["cost_cover"])
        self.assertEqual(divergent[0].measured, 47)

    def test_a_design_above_the_ceiling_diverges_instead_of_crashing(self):
        data = with_anchor(lambda a: a[18].update({"margin": "2/5"}))
        divergent = [c for c in verify_manifest(data) if not c.ok]
        self.assertEqual([c.field for c in divergent], ["feasible"])
        self.assertIn("infeasible", str(divergent[0].measured))

    def test_a_moved_frontier_is_caught(self):
        data = with_anchor(lambda a: a[20].update({"period": 1082}))
        divergent = [c for c in verify_manifest(data) if not c.ok]
        self.assertEqual([(c.field, c.measured) for c in divergent],
                         [("period", 1080)])


class TestAManifestThatCannotBeJudged(unittest.TestCase):
    def test_an_empty_manifest_is_not_a_green_light(self):
        with self.assertRaises(ManifestError):
            verify_manifest({"anchors": []})

    def test_a_manifest_without_anchors_is_refused(self):
        with self.assertRaises(ManifestError):
            verify_manifest({"about": "nothing"})

    def test_an_unknown_kind_is_refused_not_skipped(self):
        with self.assertRaises(ManifestError):
            verify_manifest({"anchors": [{"kind": "ceiling-ish", "period": 180}]})

    def test_a_missing_required_field_is_refused(self):
        with self.assertRaises(ManifestError):
            verify_manifest({"anchors": [{"kind": "design", "period": 180,
                                          "max_cost": 64, "probes": 4}]})

    def test_a_float_cannot_stand_for_a_margin(self):
        with self.assertRaises(ManifestError):
            verify_manifest({"anchors": [{"kind": "ceiling", "period": 180,
                                          "margin": 0.333}]})

    def test_a_missing_file_is_a_malformed_request(self):
        with self.assertRaises(ManifestError):
            verify_file(Path(tempfile.gettempdir()) / "no-such-manifest.json")

    def test_a_file_that_is_not_json_is_a_malformed_request(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "broken.json"
            path.write_text("{not json", encoding="utf-8")
            with self.assertRaises(ManifestError):
                verify_file(path)


class TestTheExactnessOfTheReport(unittest.TestCase):
    def test_fractions_render_as_fractions(self):
        self.assertEqual(render(Fraction(1, 3)), "1/3")
        self.assertEqual(render(Fraction(2, 33)), "2/33")
        self.assertEqual(render((12, 18, 20, 45)), "(12, 18, 20, 45)")


class TestTheGateHasTeethInTheDangerousDirection(unittest.TestCase):
    def test_a_drift_in_the_package_turns_the_shipped_manifest_red(self):
        real = acta.instance_ceiling

        def drift(period: int):
            result = real(period)
            return replace(result, margin=Fraction(1, 3)) if period == 210 else result

        with mock.patch.object(acta, "instance_ceiling", side_effect=drift):
            divergent = [c for c in verify_file(shipped_manifest()) if not c.ok]
        self.assertEqual([(c.kind, c.field) for c in divergent],
                         [("ceiling", "margin")])
        self.assertEqual(divergent[0].expected, Fraction(12, 35))
        self.assertEqual(divergent[0].measured, Fraction(1, 3))


class TestTheCommand(unittest.TestCase):
    def test_the_shipped_manifest_is_green(self):
        rc, out = run("verify-acta", str(shipped_manifest()))
        self.assertEqual(rc, 0)
        self.assertIn("0 divergent", out)

    def test_a_divergent_manifest_exits_one(self):
        data = with_anchor(lambda a: a[5].update({"margin": "5/16"}))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "acta.json"
            path.write_text(json.dumps(data), encoding="utf-8")
            rc, out = run("verify-acta", str(path))
        self.assertEqual(rc, 1)
        self.assertIn("FAIL", out)
        self.assertIn("ceiling.margin", out)

    def test_json_reports_one_result_per_field(self):
        rc, out = run("--json", "verify-acta", str(shipped_manifest()))
        self.assertEqual(rc, 0)
        data = json.loads(out)
        self.assertEqual(data["divergent"], 0)
        self.assertEqual(data["checks"], len(data["results"]))
        self.assertEqual(data["results"][0]["package"], "8/15")

    def test_a_malformed_manifest_exits_three(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.json"
            path.write_text(json.dumps({"anchors": []}), encoding="utf-8")
            rc, out = run("verify-acta", str(path))
        self.assertEqual(rc, 3)
        self.assertIn("BAD_REQUEST", out)


if __name__ == "__main__":
    unittest.main()
