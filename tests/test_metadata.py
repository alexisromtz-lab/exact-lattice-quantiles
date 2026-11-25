"""One version, three files that state it, and nothing that compared them.

The version lives in ``pyproject.toml``, ``CITATION.cff`` and ``__version__``,
and nothing read all three: a bump applied to one of them is a commit that
passes every test, ships the old number to an install, and tells a DOI the wrong
release.  No test in this suite would have noticed.  This is the reader.
"""

from __future__ import annotations

import re
import unittest
from pathlib import Path

import exact_lattice_quantiles

ROOT = Path(__file__).resolve().parents[1]


def stated_version(path: str, pattern: str) -> str:
    found = re.search(pattern, (ROOT / path).read_text(encoding="utf-8"), re.MULTILINE)
    if found is None:
        raise AssertionError(f"{path} states no version matching {pattern!r}")
    return found.group(1)


class TestOneVersion(unittest.TestCase):
    def test_the_package_the_metadata_and_the_citation_agree(self):
        self.assertEqual(
            {
                "src/exact_lattice_quantiles/__init__.py": exact_lattice_quantiles.__version__,
                "pyproject.toml": stated_version("pyproject.toml", r'^version = "([^"]+)"'),
                "CITATION.cff": stated_version("CITATION.cff", r"^version: (\S+)$"),
            },
            dict.fromkeys(
                ("src/exact_lattice_quantiles/__init__.py", "pyproject.toml", "CITATION.cff"),
                exact_lattice_quantiles.__version__),
        )

    def test_the_citation_names_a_release_date_that_parses(self):
        released = stated_version("CITATION.cff", r'^date-released: "([^"]+)"')
        year, month, day = (int(part) for part in released.split("-"))
        self.assertTrue(1 <= month <= 12 and 1 <= day <= 31 and year >= 2026, released)
