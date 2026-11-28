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


def stated(path: str, pattern: str) -> str:
    found = re.search(pattern, (ROOT / path).read_text(encoding="utf-8"), re.MULTILINE)
    if found is None:
        raise AssertionError(f"{path} states no field matching {pattern!r}")
    return found.group(1)


class TestOneVersion(unittest.TestCase):
    def test_the_package_the_metadata_and_the_citation_agree(self):
        self.assertEqual(
            {
                "src/exact_lattice_quantiles/__init__.py": exact_lattice_quantiles.__version__,
                "pyproject.toml": stated("pyproject.toml", r'^version = "([^"]+)"'),
                "CITATION.cff": stated("CITATION.cff", r"^version: (\S+)$"),
            },
            dict.fromkeys(
                ("src/exact_lattice_quantiles/__init__.py", "pyproject.toml", "CITATION.cff"),
                exact_lattice_quantiles.__version__),
        )

    def test_the_citation_names_a_release_date_that_parses(self):
        released = stated("CITATION.cff", r'^date-released: "([^"]+)"')
        year, month, day = (int(part) for part in released.split("-"))
        self.assertTrue(1 <= month <= 12 and 1 <= day <= 31 and year >= 2026, released)


class TestAnnouncementSurface(unittest.TestCase):
    """The index is a fourth place the name is stated, and nothing read it.

    The package reached PyPI while the only install line in the README was a
    ``git+https://`` URL: the artifact existed, the announcement surface did not
    name it, and every test stayed green.  An artifact nobody is told how to
    install is published and mute at the same time.
    """

    def test_the_readme_installs_the_package_by_the_name_the_manifest_declares(self):
        name = stated("pyproject.toml", r'^name = "([^"]+)"')
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        missing = [c for c in (f"pip install {name}", f"uvx {name}") if c not in readme]
        self.assertEqual(missing, [], f"README announces none of {missing}")
