"""The edge nobody measured: manifest against the acta it quotes.

Every test here writes its own acta, so none of them depends on the corpus that
produced the manifest: the mechanism is what is under test, not one file.
"""

from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from exact_lattice_quantiles import acta
from exact_lattice_quantiles.acta import (SourceError, inspect_source,
                                          parse_tables, shipped_manifest)
from exact_lattice_quantiles.ceiling import instance_ceiling
from exact_lattice_quantiles.cli import main


def run(*argv: str) -> tuple[int, str]:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = main(list(argv))
    return rc, buf.getvalue()


ACTA = """\
from fractions import Fraction as F

TECHOS = {
    180: F(1, 3),
    252: F(2, 7),
}
"""

MANIFEST = {
    "about": "synthetic: one table, two periods, both anchored",
    "source": {
        "path": "acta_sintetica.py",
        "blob": "",
        "lines": 0,
        "tables": [{"name": "TECHOS", "line": 3, "columns": ["margin"]}],
    },
    "anchors": [
        {"kind": "ceiling", "period": 180, "margin": "1/3"},
        {"kind": "ceiling", "period": 252, "margin": "2/7"},
    ],
}


class Planted(unittest.TestCase):
    """A manifest and an acta on disk, both writable by the test alone."""

    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.root = Path(self._dir.name)
        self.addCleanup(self._dir.cleanup)

    def plant(self, text: str = ACTA, **manifest) -> tuple[Path, Path]:
        source = self.root / "acta_sintetica.py"
        source.write_text(text, encoding="utf-8")
        data = json.loads(json.dumps(MANIFEST))
        data["source"]["blob"] = acta.read_blob(source)
        data["source"]["lines"] = text.count("\n") + (0 if text.endswith("\n") else 1)
        for key, value in manifest.items():
            data[key] = value
        path = self.root / "manifest.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        return path, source


class TestTheEdgeIsMeasured(Planted):
    def test_an_acta_that_agrees_reproduces_every_cell(self):
        manifest, source = self.plant()
        report = inspect_source(source, json.loads(manifest.read_text()))
        self.assertFalse(report.moved)
        self.assertEqual(report.divergent, ())
        self.assertEqual([c.field for c in report.checks],
                         ["TECHOS[180].margin", "TECHOS[252].margin"])

    def test_the_blob_and_the_lines_are_read_from_the_file(self):
        manifest, source = self.plant()
        report = inspect_source(source, json.loads(manifest.read_text()))
        self.assertEqual(report.blob, acta.read_blob(source))
        self.assertEqual(report.lines, ACTA.count("\n"))

    def test_a_changed_cell_is_red_and_names_the_cell(self):
        manifest, source = self.plant(ACTA.replace("F(1, 3)", "F(8, 15)"))
        report = inspect_source(source, json.loads(manifest.read_text()))
        self.assertEqual([c.field for c in report.divergent], ["TECHOS[180].margin"])
        self.assertIn("8/15", str(report.divergent[0]))

    def test_a_moved_blob_alone_does_not_turn_the_numbers_red(self):
        manifest, source = self.plant()
        data = json.loads(manifest.read_text())
        data["source"]["blob"] = "0" * 40
        report = inspect_source(source, data)
        self.assertTrue(report.moved)
        self.assertEqual(report.divergent, ())

    def test_an_unanchored_cell_is_named_not_passed(self):
        manifest, source = self.plant(ACTA.replace("F(2, 7),", "F(2, 7),\n    396: F(1, 4),"))
        report = inspect_source(source, json.loads(manifest.read_text()))
        self.assertEqual(report.unanchored, ("TECHOS[396].margin",))
        self.assertFalse(report.moved)
        self.assertEqual(report.divergent, ())


class TestWhatCannotBeMeasuredIsNotGreen(Planted):
    def test_a_float_in_a_table_is_refused(self):
        manifest, source = self.plant(ACTA.replace("F(1, 3)", "0.3333333"))
        with self.assertRaises(SourceError) as caught:
            inspect_source(source, json.loads(manifest.read_text()))
        self.assertIn("not exact", str(caught.exception))

    def test_a_table_computed_at_run_time_is_refused(self):
        text = ACTA + "\nOTRA = {180: TECHOS[180]}\n"
        with self.assertRaises(SourceError) as caught:
            parse_tables(text, ["OTRA"])
        self.assertIn("never imports the acta", str(caught.exception))

    def test_a_table_the_manifest_cites_and_the_acta_lacks_is_refused(self):
        manifest, source = self.plant()
        data = json.loads(manifest.read_text())
        data["source"]["tables"] = [{"name": "BAJAN", "line": 1, "columns": ["margin"]}]
        with self.assertRaises(SourceError) as caught:
            inspect_source(source, data)
        self.assertIn("BAJAN", str(caught.exception))

    def test_a_manifest_without_source_cannot_measure_one(self):
        manifest, source = self.plant()
        data = json.loads(manifest.read_text())
        del data["source"]
        with self.assertRaises(SourceError) as caught:
            inspect_source(source, data)
        self.assertIn("no 'source'", str(caught.exception))

    def test_a_column_this_package_cannot_compare_is_refused(self):
        manifest, source = self.plant()
        data = json.loads(manifest.read_text())
        data["source"]["tables"] = [{"name": "TECHOS", "line": 3, "columns": ["roughly"]}]
        data["anchors"][0]["roughly"] = 1
        with self.assertRaises(SourceError) as caught:
            inspect_source(source, data)
        self.assertIn("cannot compare", str(caught.exception))

    def test_a_column_no_anchor_declares_is_a_named_blind_spot(self):
        manifest, source = self.plant()
        data = json.loads(manifest.read_text())
        data["source"]["tables"] = [{"name": "TECHOS", "line": 3, "columns": ["roughly"]}]
        report = inspect_source(source, data)
        self.assertEqual(report.checks, ())
        self.assertEqual(report.unanchored,
                         ("TECHOS[180].roughly", "TECHOS[252].roughly"))

    def test_a_missing_file_is_no_medido(self):
        manifest, _ = self.plant()
        with self.assertRaises(SourceError):
            inspect_source(self.root / "no_existe.py", json.loads(manifest.read_text()))


class TestTheCommand(Planted):
    def test_without_source_the_verdict_is_unchanged(self):
        rc, out = run("verify-acta", str(shipped_manifest()))
        self.assertEqual(rc, 0)
        self.assertIn("84 checks over 24 anchors", out)

    def test_a_divergent_acta_exits_one(self):
        manifest, _ = self.plant(ACTA.replace("F(1, 3)", "F(8, 15)"))
        rc, out = run("verify-acta", str(manifest), "--source",
                      str(self.root / "acta_sintetica.py"))
        self.assertEqual(rc, 1)
        self.assertIn("TECHOS[180].margin", out)

    def test_an_unmeasurable_acta_exits_four_and_says_no_medido(self):
        manifest, _ = self.plant(ACTA.replace("F(1, 3)", "0.3333333"))
        rc, out = run("verify-acta", str(manifest), "--source",
                      str(self.root / "acta_sintetica.py"))
        self.assertEqual(rc, 4)
        self.assertIn("NO_MEDIDO", out)

    def test_a_moved_blob_is_reported_as_moved_with_the_numbers_holding(self):
        manifest, source = self.plant()
        data = json.loads(manifest.read_text())
        data["source"]["blob"] = "0" * 40
        manifest.write_text(json.dumps(data), encoding="utf-8")
        rc, out = run("verify-acta", str(manifest), "--source", str(source))
        self.assertEqual(rc, 0)
        self.assertIn("MOVIDO", out)
        self.assertIn("stale", out)

    def test_a_period_keyed_table_from_a_later_sweep_is_red_when_it_drifts(self):
        acta_text = "from fractions import Fraction as F\n\nTECHOS = {720: F(4, 15)}\n"
        manifest, source = self.plant(
            acta_text,
            anchors=[{"kind": "ceiling", "period": 720, "margin": "4/15"}],
        )
        rc, out = run("verify-acta", str(manifest), "--source", str(source))
        self.assertEqual(rc, 0)
        manifest, source = self.plant(
            acta_text.replace("F(4, 15)", "F(4, 13)"),
            anchors=[{"kind": "ceiling", "period": 720, "margin": "4/15"}],
        )
        rc, out = run("verify-acta", str(manifest), "--source", str(source))
        self.assertEqual(rc, 1)
        self.assertIn("TECHOS[720].margin", out)

    def test_json_carries_the_source_block_and_the_subject_per_check(self):
        manifest, source = self.plant()
        rc, out = run("--json", "verify-acta", str(manifest), "--source", str(source))
        self.assertEqual(rc, 0)
        payload = json.loads(out)
        self.assertFalse(payload["source"]["moved"])
        self.assertEqual(payload["source"]["tables"], ["TECHOS"])
        edge = [r for r in payload["results"] if r["kind"] == "source"]
        self.assertEqual({r["subject"] for r in edge}, {"acta"})
        self.assertEqual(edge[0]["acta"], "1/3")


OUT_OF_SAMPLE = (720, 840, 1260, 2520)


class TestTheOutOfSampleClaimIsOnBothEdges(unittest.TestCase):
    """The acta claims the factorised ceiling holds above 400. Both edges.

    Anchoring those rows means the package recomputes the ceiling *and* the
    acta's written cell is re-read; a claim resting on one edge alone is a
    citation nobody checks.
    """

    def data(self) -> dict:
        return json.loads(Path(shipped_manifest()).read_text())

    def test_the_four_periods_above_four_hundred_are_anchored(self):
        anchored = {a["period"] for a in self.data()["anchors"]
                    if a["kind"] == "ceiling"}
        self.assertEqual(anchored & set(OUT_OF_SAMPLE), set(OUT_OF_SAMPLE))

    def test_every_out_of_sample_anchor_is_the_ceiling_the_package_computes(self):
        for anchor in self.data()["anchors"]:
            if anchor["kind"] == "ceiling" and anchor["period"] in OUT_OF_SAMPLE:
                with self.subTest(period=anchor["period"]):
                    self.assertEqual(Fraction(anchor["margin"]),
                                     instance_ceiling(anchor["period"]).margin)


if __name__ == "__main__":
    unittest.main()
