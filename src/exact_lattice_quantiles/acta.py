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
    ladder    period, design, probes, margin
    frontier  period                    max_period, reduced, witness, bound,
                                        cover, saving, breaks

``ladder`` takes the design as an **input**, unlike ``design`` where ``design``
is an echoed output: a row of a margin ladder is quoted per ``(design, k)``, and
the design is the half of that pair the manifest has to state.

``cost_cover`` is the cheapest design that covers the universe at margin ``0``;
``cost_ceiling`` is the cheapest that reaches the ceiling itself.  Both are
``sum(M // 2)`` over the design, the cost axis of this package.
"""

from __future__ import annotations

import ast
import json
import subprocess
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

from .ceiling import InfeasibleAboveInstanceCeiling, instance_ceiling
from .lattice import cell_widths, grids, seats, universe
from .coverage import best_probes, covered_divisors, minimum_cost_design, witness_margin
from .pricing import exactness_frontier
from .solve import InfeasibleWithinBudget, solve_lattice_schedule

KINDS = ("ceiling", "design", "ladder", "frontier")


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
    subject: str = "package"

    @property
    def ok(self) -> bool:
        return self.expected == self.measured

    def __str__(self) -> str:
        mark = "ok  " if self.ok else "FAIL"
        return (f"{mark}  #{self.anchor} {self.kind}.{self.field}: "
                f"manifest {render(self.expected)}  {self.subject} {render(self.measured)}")


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


def _check_ladder(index: int, anchor: dict) -> list[Check]:
    """One row of a margin ladder: the margin ``k`` shared sections buy on a design.

    The design is stated by the anchor rather than looked up, because a ladder
    row is a claim about a *fixed* design; a row that let the design float would
    be the cheapest-design question the ``design`` kind already answers.
    """
    where = f"anchor {index} (ladder)"
    period = _period(_require(anchor, "period", where), f"{where}: period")
    design = _ints(_require(anchor, "design", where), f"{where}: design")
    probes = _int(_require(anchor, "probes", where), f"{where}: probes")
    margin = _exact(_require(anchor, "margin", where), f"{where}: margin")
    U = universe(period)
    carried = {d for M in design for d, _, _ in seats(M)}
    missing = [d for d in U if d not in carried]
    if missing:
        raise ManifestError(
            f"{where}: the design carries no block for {missing} of U({period}), "
            "so its ladder would report the margin of a divisor nobody holds")
    reachable, witness = best_probes(design, U, probes)
    held: Fraction | None = None
    if witness and set(covered_divisors(design, witness, U)) == set(U):
        held = reachable
    return [Check(index, "ladder", "margin", margin, held)]


_CHECKERS = {
    "ceiling": _check_ceiling,
    "design": _check_design,
    "ladder": _check_ladder,
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


# ── the source: a blob declared and never read is a copied field ──────────────
#
# The manifest has always carried ``source.blob`` and ``source.lines``, and for a
# while nothing read them.  That is the defect this section closes: a pointer
# whose target is never measured is a copy, and two copies diverge.  The reader
# below measures the acta itself, so the chain *acta -> manifest -> package* has
# both of its edges checked by an instrument instead of by memory.

_EXACT_CALLS = ("F", "Fraction")
_INT_FIELDS = frozenset({"grids", "universe_size", "critical_divisor",
                         "cost_cover", "cost_ceiling"})


class SourceError(Exception):
    """The acta cannot be measured at all.

    Distinct from a divergence on purpose: a divergence (RC 1) is a measurement
    that says *no*; this is a measurement that could not be taken, and it must
    never be read as agreement.
    """


@dataclass(frozen=True)
class SourceReport:
    """What the acta says, next to what the manifest claims it says."""

    path: str
    blob: str
    lines: int
    declared_blob: object
    declared_lines: object
    tables: tuple[str, ...]
    checks: tuple[Check, ...]
    unanchored: tuple[str, ...]

    @property
    def moved(self) -> bool:
        """The file is not the object the manifest cites."""
        return self.declared_blob is not None and self.declared_blob != self.blob

    @property
    def divergent(self) -> tuple[Check, ...]:
        return tuple(check for check in self.checks if not check.ok)


def _exact_node(node: ast.AST, where: str) -> object:
    """Evaluate one literal of an acta, exactly.  No names, no calls but ``F``."""
    if isinstance(node, ast.Constant):
        if isinstance(node.value, bool) or isinstance(node.value, float):
            raise SourceError(
                f"{where}: {node.value!r} is not exact; a table quoted as a "
                "measurement holds integers and fractions, never a rounding")
        if isinstance(node.value, (int, str)):
            return node.value
        raise SourceError(f"{where}: {type(node.value).__name__} is not exact arithmetic")
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        value = _exact_node(node.operand, where)
        if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
            raise SourceError(f"{where}: cannot negate {value!r}")
        return -value
    if isinstance(node, ast.Dict):
        out: dict[object, object] = {}
        for key_node, value_node in zip(node.keys, node.values):
            if key_node is None:
                raise SourceError(f"{where}: dict unpacking is not a table")
            out[_exact_node(key_node, where)] = _exact_node(value_node, where)
        return out
    if isinstance(node, (ast.Tuple, ast.List)):
        return tuple(_exact_node(element, where) for element in node.elts)
    if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id in _EXACT_CALLS and not node.keywords):
        args = [_exact_node(arg, where) for arg in node.args]
        if not 1 <= len(args) <= 2:
            raise SourceError(f"{where}: {node.func.id}() takes one or two arguments")
        for arg in args:
            if isinstance(arg, bool) or not isinstance(arg, int):
                raise SourceError(f"{where}: {node.func.id}{tuple(args)!r} is not a rational")
        return Fraction(*args)
    raise SourceError(
        f"{where}: this is not a literal.  The reader evaluates literals only and "
        "never imports the acta, so a table computed at run time has to be written "
        "down before it can be quoted")


def parse_tables(text: str, names: Sequence[str]) -> dict[str, object]:
    """Read ``name = <literal>`` at module level.  The acta is parsed, not imported."""
    try:
        tree = ast.parse(text)
    except SyntaxError as exc:
        raise SourceError(f"the acta is not Python: {exc}") from None
    found: dict[str, object] = {}
    for node in tree.body:
        target = None
        value = None
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target, value = node.targets[0], node.value
        elif isinstance(node, ast.AnnAssign):
            target, value = node.target, node.value
        if isinstance(target, ast.Name) and target.id in names and value is not None:
            found[target.id] = _exact_node(value, f"{target.id} (l. {node.lineno})")
    missing = [name for name in names if name not in found]
    if missing:
        raise SourceError("the acta assigns no " + ", ".join(missing)
                          + ": the manifest cites a table that is not there")
    return found


def _declared_tables(data: object) -> list[dict]:
    source = data.get("source") if isinstance(data, dict) else None
    if not isinstance(source, dict):
        raise SourceError("the manifest declares no 'source' to measure the acta against")
    tables = source.get("tables")
    if not isinstance(tables, list) or not tables:
        raise SourceError("the manifest's 'source' names no 'tables'")
    out: list[dict] = []
    for index, spec in enumerate(tables):
        if not isinstance(spec, dict):
            raise SourceError(f"source.tables[{index}]: must be a JSON object")
        name, columns = spec.get("name"), spec.get("columns")
        if not isinstance(name, str):
            raise SourceError(f"source.tables[{index}]: needs a string 'name'")
        if not isinstance(columns, list) or not columns:
            raise SourceError(f"source.tables[{index}]: needs a non-empty 'columns' list")
        for column in columns:
            if not isinstance(column, str):
                raise SourceError(f"source.tables[{index}]: column names must be strings")
        out.append(spec)
    return out


def _anchor_exact(anchor: dict, column: str, where: str) -> object:
    if column in _INT_FIELDS:
        return _int(anchor[column], f"{where}: {column}")
    if column in ("margin", "robust_margin"):
        return _exact(anchor[column], f"{where}: {column}")
    raise SourceError(f"{where}: the manifest's tables declare a column {column!r} "
                      "this package cannot compare")


def anchor_table_checks(data: object, raw: dict[str, object]) -> tuple[list[Check], tuple[str, ...]]:
    """Compare every anchor field that the acta's tables also declare.

    The manifest was already checked against the package; this is the missing
    edge, manifest against acta.  Cells that no anchor declares come back as
    ``unanchored``: an unmeasured cell is a blind spot with a name, not a pass.
    """
    anchors = data.get("anchors") if isinstance(data, dict) else None
    if not isinstance(anchors, list):
        raise SourceError("manifest has no 'anchors' list")
    checks: list[Check] = []
    covered: set[tuple[str, object, str]] = set()
    for spec in _declared_tables(data):
        name, columns = spec["name"], tuple(spec["columns"])
        rows = raw[name]
        if not isinstance(rows, dict):
            raise SourceError(f"{name}: a table cited by period must map periods to rows")
        for index, anchor in enumerate(anchors):
            if not isinstance(anchor, dict) or anchor.get("kind") != "ceiling":
                continue
            period = anchor.get("period")
            if isinstance(period, bool) or not isinstance(period, int) or period not in rows:
                continue
            row = rows[period]
            values = row if len(columns) > 1 else (row,)
            if not isinstance(values, tuple) or len(values) != len(columns):
                raise SourceError(
                    f"{name}[{period}]: the manifest declares {len(columns)} columns, "
                    f"the acta's row has {len(values) if isinstance(values, tuple) else 1}")
            for column, value in zip(columns, values):
                if column not in anchor:
                    continue
                covered.add((name, period, column))
                expected = _anchor_exact(anchor, column, f"anchor {index} (ceiling)")
                checks.append(Check(index, "source", f"{name}[{period}].{column}",
                                    expected, value, subject="acta"))
    unanchored: list[str] = []
    for spec in _declared_tables(data):
        name, columns = spec["name"], tuple(spec["columns"])
        for period in raw[name]:
            for column in columns:
                if (name, period, column) not in covered:
                    unanchored.append(f"{name}[{period}].{column}")
    return checks, tuple(unanchored)


def read_blob(path: str | Path) -> str:
    """The git object id of the file: the address the manifest cites."""
    try:
        done = subprocess.run(["git", "hash-object", str(path)],
                              capture_output=True, text=True)
    except OSError as exc:
        raise SourceError(f"cannot run git to measure the blob: {exc}") from None
    if done.returncode != 0:
        raise SourceError(f"git hash-object refused {path}: "
                          f"{done.stderr.strip() or done.returncode}")
    blob = done.stdout.strip()
    if len(blob) != 40:
        raise SourceError(f"git hash-object returned {blob!r}, not an object id")
    return blob


def inspect_source(path: str | Path, data: object) -> SourceReport:
    """Measure a live acta against the manifest that quotes it."""
    target = Path(path)
    try:
        text = target.read_text(encoding="utf-8")
    except OSError as exc:
        raise SourceError(f"cannot read {path}: {exc.strerror or exc}") from None
    declared = data.get("source") if isinstance(data, dict) else None
    if not isinstance(declared, dict):
        raise SourceError("the manifest declares no 'source' to measure the acta against")
    specs = _declared_tables(data)
    raw = parse_tables(text, [spec["name"] for spec in specs])
    checks, unanchored = anchor_table_checks(data, raw)
    return SourceReport(
        path=str(target),
        blob=read_blob(target),
        lines=text.count("\n") + (0 if text.endswith("\n") else 1),
        declared_blob=declared.get("blob"),
        declared_lines=declared.get("lines"),
        tables=tuple(spec["name"] for spec in specs),
        checks=tuple(checks),
        unanchored=unanchored,
    )
