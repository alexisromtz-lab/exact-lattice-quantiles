"""One call: from a period and a budget to a schedule, or to a proof of no.

The contract is worth stating precisely, because the two budgets are not
interchangeable and conflating them is the classic mistake.

* ``max_cost`` is a budget on **grids**, measured in points (``M // 2`` each).
* ``k_probes`` is a budget on **shared sections**, and a section is not a point.
* ``min_margin`` is a requirement on the **design**: the margin the grids
  guarantee when each divisor may be read in whichever of them suits it best.

A design meeting ``min_margin`` may still need more than ``k_probes`` sections
to *realise* that margin, and the returned schedule says exactly how many it
needs.  That is not a failure of the design and it is not a silent success:
it is a named outcome, ``PROBES_INSUFFICIENT``.

Above ``nu(T)`` nothing helps.  The exception carries the ceiling, so the
caller learns the number that made the request impossible instead of a
negative search result.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .ceiling import Ceiling, InfeasibleAboveInstanceCeiling, instance_ceiling
from .coverage import (best_probes, design_margin, drift_radius,
                       margin_ladder, minimum_cost_design, minimum_probes,
                       witness_margin)
from .lattice import cell_widths, grids, universe

OK = "OK"
PROBES_INSUFFICIENT = "PROBES_INSUFFICIENT"


class InfeasibleWithinBudget(Exception):
    """The margin exists but not for ``max_cost`` points.

    Distinct from :class:`InfeasibleAboveInstanceCeiling` on purpose: this one
    can be fixed by paying, the other one cannot be fixed at all.
    """

    rc = 2

    def __init__(self, period: int, needed: Fraction, cheapest: int, max_cost: int) -> None:
        self.period = period
        self.needed = needed
        self.cheapest = cheapest
        self.max_cost = max_cost
        super().__init__(
            f"margin {needed} needs at least {cheapest} points for period "
            f"{period}, and the budget is {max_cost}"
        )


@dataclass(frozen=True)
class Schedule:
    """A design, the sections that serve it, and what each of them proves.

    ``design_margin`` is the P51 currency: each divisor read in its best grid.
    ``witness_margin`` is what the ``k`` chosen sections actually hold, with
    each divisor leaning on its own pair.  ``robust_margin`` is the exact
    radius of common drift the sections still *cover* -- a weaker question than
    the margin, because drift may push a section into a block it never captured.
    Quote the witness margin when a schedule is being promised; quote
    ``drift_tolerance`` when the question is how much common clock error the
    sections keep that margin under.
    """

    period: int
    max_cost: int
    k_probes: int
    design: tuple[int, ...]
    cost: int
    design_margin: Fraction
    probes: tuple[Fraction, ...]
    witness_margin: Fraction
    robust_margin: Fraction
    probes_needed: int | None
    verdict: str
    ceiling: Ceiling

    @property
    def meets_requirement(self) -> bool:
        """True when the realised margin, not just the design margin, holds."""
        return self.verdict == OK

    def as_dict(self) -> dict[str, object]:
        """A JSON-safe view: every rational as ``"n/d"``, never a float."""
        return {
            "period": self.period,
            "max_cost": self.max_cost,
            "k_probes": self.k_probes,
            "design": list(self.design),
            "cost": self.cost,
            "design_margin": str(self.design_margin),
            "probes": [str(q) for q in self.probes],
            "witness_margin": str(self.witness_margin),
            "robust_margin": str(self.robust_margin),
            "probes_needed": self.probes_needed,
            "verdict": self.verdict,
            "meets_requirement": self.meets_requirement,
            "ceiling": {
                "margin": str(self.ceiling.margin),
                "critical_divisor": self.ceiling.critical_divisor,
                "dropped": list(self.ceiling.dropped),
                "prime_product": str(self.ceiling.prime_product),
                "radical_bound_hit": self.ceiling.radical_bound_hit,
                "universe": list(self.ceiling.universe),
            },
        }


def solve_lattice_schedule(period: int, max_cost: int, k_probes: int,
                           min_margin: Fraction | int = Fraction(0)) -> Schedule:
    """Cheapest design of ``period`` reaching ``min_margin``, served by ``k_probes``.

    Raises :class:`InfeasibleAboveInstanceCeiling` when ``min_margin`` exceeds
    ``nu(T)`` -- no design, at any price -- and :class:`InfeasibleWithinBudget`
    when the cheapest design reaching it exceeds ``max_cost``.
    """
    if k_probes < 1:
        raise ValueError(f"a schedule needs at least one section, got {k_probes}")
    required = Fraction(min_margin)
    ceiling = instance_ceiling(period)
    if required > ceiling.margin:
        raise InfeasibleAboveInstanceCeiling(period, required, ceiling.margin)

    gr = grids(period)
    U = universe(period)
    widths = cell_widths(gr, U)
    cost, design = minimum_cost_design(gr, U, widths, required)
    if cost is None:  # pragma: no cover - ruled out by the ceiling check above
        raise InfeasibleAboveInstanceCeiling(period, required, ceiling.margin)
    if cost > max_cost:
        raise InfeasibleWithinBudget(period, required, cost, max_cost)

    _, sections = best_probes(design, U, k_probes)
    probes = tuple(sections)
    realised = min(witness_margin(design, probes, U).values()) if probes else Fraction(0)
    held = min(drift_radius(design, probes, U).values()) if probes else Fraction(0)
    needed, _ = minimum_probes(design, U, required) if required > 0 else (0, ())

    if required > 0 and (needed is None or needed > k_probes):
        verdict = PROBES_INSUFFICIENT
    else:
        verdict = OK

    return Schedule(
        period=period,
        max_cost=max_cost,
        k_probes=k_probes,
        design=design,
        cost=cost,
        design_margin=design_margin(design, U, widths),
        probes=probes,
        witness_margin=realised,
        robust_margin=held,
        probes_needed=needed,
        verdict=verdict,
        ceiling=ceiling,
    )


def ladder(period: int) -> dict[int, Fraction]:
    """``cost -> best design margin`` over the whole instance.

    The curve quoted whenever the instance is argued about: the cheapest
    point, where buying stops helping, and the ceiling the ladder converges
    to and cannot pass.
    """
    gr = grids(period)
    U = universe(period)
    return margin_ladder(gr, U, cell_widths(gr, U))
