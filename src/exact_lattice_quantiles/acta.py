"""Recompute an acta's anchors: the boundary as a gate, not as memory.

This package is quoted from places that keep their own record of what it says --
an *acta*: a table of instance ceilings, the price of reaching one, the first
period a closed form misprices.  Nothing recomputes those tables, so the record
and the library can drift apart while both stay green: the acta keeps its
numbers, the library keeps evolving, and the only witness to the boundary is
somebody remembering it.

A **manifest** writes the boundary down as JSON: a list of anchors, each one a
claim about one quantity at one period.  :func:`verify_manifest` recomputes every
anchor in exact arithmetic and yields one :class:`Check` **per field**, so a
divergence names the quantity instead of the file.

Three rules keep the reader from being decorative:

* exact values only.  A JSON float is refused: it cannot express ``1/3``, so the
  comparison would be a rounding experiment rather than a measurement;
* an empty manifest is an error, not a pass -- verifying nothing is not a
  green light;
* an unknown kind or a missing required field is an error, not a skip -- a check
  that silently disappears is worse than one that fails.

Anchor kinds.  Each states what it needs and what it may add::

    ceiling   period, margin            critical_divisor, grids, universe_size,
                                        cost_cover, cost_ceiling
    design    period, max_cost, probes, margin
                                        cost, design, probes_needed,
                                        robust_margin, verdict
    frontier  period                    max_period, reduced, witness, bound,
                                        cover, saving, breaks

``cost_cover`` is the cheapest design that covers the universe at margin ``0``;
``cost_ceiling`` is the cheapest that reaches the ceiling itself.  Both are
``sum(M // 2)`` over the design, the cost axis of this package.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

from .ceiling import InfeasibleAboveInstanceCeiling, instance_ceiling
from .lattice import cell_widths, grids, universe
from .coverage import minimum_cost_design
from .pricing import exactness_frontier
from .solve import InfeasibleWithinBudget, solve_lattice_schedule

KINDS = ("ceiling", "design", "frontier")


class ManifestError(Exception):
    """The manifest cannot be read as a manifest: the request is malformed."""


@dataclass(frozen=True)
class Check:
    """One field of one anchor, recomputed.

    ``expected`` is what the manifest claims, ``measured`` what the package
    computes now.  Both are exact (``Fraction``, ``int``, ``str`` or a tuple of
    those), never floats.
    """

    anchor: int
    kind: str
    field: str
    expected: object
    measured: object

    @property
    def ok(self) -> bool:
        return self.expected == self.measured

    def __str__(self) -> str:
        mark = "ok  " if self.ok else "FAIL"
        return (f"{mark}  #{self.anchor} {self.kind}.{self.field}: "
                f"manifest {render(self.expected)}  package {render(self.measured)}")


def render(value: object) -> str:
    """An exact value as text: fractions stay fractions, never decimals."""
    if isinstance(value, tuple):
        return "(" + ", ".join(render(v) for v in value) + ")"
    return str(value)


def _require(anchor: dict, key: str, where: str) -> object:
    if key not in anchor:
        raise ManifestError(f"{where}: missing required field {key!r}")
    return anchor[key]


def _int(value: object, where: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ManifestError(f"{where}: integer expected, got {value!r}")
    return value


def _period(value: object, where: str) -> int:
    period = _int(value, where)
    if period < 2:
        raise ManifestError(f"{where}: period must be at least 2, got {period}")
    return period


def _exact(value: object, where: str) -> Fraction:
    if isinstance(value, bool) or isinstance(value, float):
        raise ManifestError(
            f"{where}: exact value expected, got the float {value!r}; "
            "a fraction is written as a string, e.g. \"1/3\"")
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str):
        try:
            return Fraction(value)
        except (ValueError, ZeroDivisionError):
            raise ManifestError(f"{where}: {value!r} is not a rational") from None
    raise ManifestError(f"{where}: exact value expected, got {type(value).__name__}")


def _ints(value: object, where: str) -> tuple[int, ...]:
    if not isinstance(value, list):
        raise ManifestError(f"{where}: list of integers expected, got {value!r}")
    return tuple(sorted(_int(v, f"{where}[{i}]") for i, v in enumerate(value)))


def _check_ceiling(index: int, anchor: dict) -> list[Check]:
    where = f"anchor {index} (ceiling)"
    period = _period(_require(anchor, "period", where), f"{where}: period")
    result = instance_ceiling(period)
    margin = _exact(_require(anchor, "margin", where), f"{where}: margin")
    checks = [Check(index, "ceiling", "margin", margin, result.margin)]
    if "critical_divisor" in anchor:
        checks.append(Check(index, "ceiling", "critical_divisor",
                            _int(anchor["critical_divisor"], f"{where}: critical_divisor"),
                            result.critical_divisor))
    if "grids" in anchor:
        checks.append(Check(index, "ceiling", "grids",
                            _int(anchor["grids"], f"{where}: grids"),
                            len(grids(period))))
    if "universe_size" in anchor:
        checks.append(Check(index, "ceiling", "universe_size",
                            _int(anchor["universe_size"], f"{where}: universe_size"),
                            len(result.universe)))
    if "cost_cover" in anchor or "cost_ceiling" in anchor:
        period_grids, divisors = grids(period), universe(period)
        width = cell_widths(period_grids, divisors)
        for field, cut in (("cost_cover", Fraction(0)), ("cost_ceiling", result.margin)):
            if field in anchor:
                cost, _ = minimum_cost_design(period_grids, divisors, width, cut)
                checks.append(Check(index, "ceiling", field,
                                    _int(anchor[field], f"{where}: {field}"), cost))
    return checks


def _check_design(index: int, anchor: dict) -> list[Check]:
    where = f"anchor {index} (design)"
    period = _period(_require(anchor, "period", where), f"{where}: period")
    max_cost = _int(_require(anchor, "max_cost", where), f"{where}: max_cost")
    probes = _int(_require(anchor, "probes", where), f"{where}: probes")
    margin = _exact(_require(anchor, "margin", where), f"{where}: margin")
    try:
        schedule = solve_lattice_schedule(period, max_cost, probes, margin)
    except (InfeasibleAboveInstanceCeiling, InfeasibleWithinBudget) as exc:
        return [Check(index, "design", "feasible", True, f"infeasible: {exc}")]
    checks = [Check(index, "design", "witness_margin", margin, schedule.witness_margin)]
    if "cost" in anchor:
        checks.append(Check(index, "design", "cost",
                            _int(anchor["cost"], f"{where}: cost"), schedule.cost))
    if "probes_needed" in anchor:
        checks.append(Check(index, "design", "probes_needed",
                            _int(anchor["probes_needed"], f"{where}: probes_needed"),
                            schedule.probes_needed))
    if "design" in anchor:
        checks.append(Check(index, "design", "design",
                            _ints(anchor["design"], f"{where}: design"),
                            tuple(sorted(schedule.design))))
    if "robust_margin" in anchor:
        checks.append(Check(index, "design", "robust_margin",
                            _exact(anchor["robust_margin"], f"{where}: robust_margin"),
                            schedule.robust_margin))
    if "verdict" in anchor:
        if not isinstance(anchor["verdict"], str):
            raise ManifestError(f"{where}: verdict must be a string, got {anchor['verdict']!r}")
        checks.append(Check(index, "design", "verdict", anchor["verdict"], schedule.verdict))
    return checks


def _check_frontier(index: int, anchor: dict) -> list[Check]:
    where = f"anchor {index} (frontier)"
    period = _period(_require(anchor, "period", where), f"{where}: period")
    max_period = _int(anchor.get("max_period", 8000), f"{where}: max_period")
    reduced = anchor.get("reduced", True)
    if not isinstance(reduced, bool):
        raise ManifestError(f"{where}: reduced must be true or false, got {reduced!r}")
    found = exactness_frontier(max_period, reduced=reduced)
    checks = [Check(index, "frontier", "period", period, found["period"])]
    if "witness" in anchor:
        checks.append(Check(index, "frontier", "witness",
                            _ints(anchor["witness"], f"{where}: witness"),
                            tuple(sorted(found["witness"]))))
    for field in ("bound", "cover", "saving", "breaks", "scanned"):
        if field in anchor:
            checks.append(Check(index, "frontier", field,
                                _int(anchor[field], f"{where}: {field}"), found[field]))
    return checks


_CHECKERS = {
    "ceiling": _check_ceiling,
    "design": _check_design,
    "frontier": _check_frontier,
}


def verify_manifest(data: object) -> list[Check]:
    """Recompute every anchor of a parsed manifest.

    Raises :class:`ManifestError` when the manifest cannot be judged at all --
    not a JSON object, no ``anchors`` list, an empty list, an unknown kind or a
    required field missing.  A manifest that breaks no rule of form but whose
    numbers do not reproduce returns checks with ``ok`` false; that is a
    divergence (RC 1), not a malformed request (RC 3).
    """
    if not isinstance(data, dict):
        raise ManifestError(f"manifest must be a JSON object, got {type(data).__name__}")
    anchors = data.get("anchors")
    if not isinstance(anchors, list):
        raise ManifestError("manifest has no 'anchors' list")
    if not anchors:
        raise ManifestError("manifest lists no anchors: an empty manifest verifies nothing")
    checks: list[Check] = []
    for index, anchor in enumerate(anchors):
        if not isinstance(anchor, dict):
            raise ManifestError(f"anchor {index}: must be a JSON object")
        kind = anchor.get("kind")
        if kind not in _CHECKERS:
            raise ManifestError(
                f"anchor {index}: unknown kind {kind!r}; known kinds are "
                + ", ".join(KINDS))
        checks.extend(_CHECKERS[kind](index, anchor))
    return checks


def load_manifest(path: str | Path) -> object:
    """Read a manifest from disk.  Unreadable or non-JSON is a malformed request."""
    try:
        text = Path(path).read_text(encoding="utf-8")
    except OSError as exc:
        raise ManifestError(f"cannot read {path}: {exc.strerror or exc}") from None
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise ManifestError(f"{path} is not JSON: {exc}") from None


def verify_file(path: str | Path) -> list[Check]:
    """Read a manifest and recompute every anchor in it."""
    return verify_manifest(load_manifest(path))


def shipped_manifest() -> Path:
    """The manifest this package ships: the acta it is quoted by, field by field."""
    return Path(__file__).with_name("manifests") / "acta_margen_secciones_20261125.json"
