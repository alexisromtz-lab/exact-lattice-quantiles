"""Coverage duality: the ceiling as a hypergraph problem, not as a search.

Two currencies and one dual:

* **Cost.**  Grid ``M`` costs ``M // 2`` points.  A design is a set of grids
  and its cost is the sum.
* **Witness margin.**  Under a set ``Q`` of shared sections, divisor ``d``
  survives if some pair ``(M, q)`` holds it far from the edges of its block.
  Each ``d`` picks its own pair, so this is a witness condition, not a
  robust one.
* **The dual.**  "margin at least ``v``" is "the grids whose blocks are at
  least ``v`` wide cover the universe ``U(T)``".  A coverage condition.  Its
  minimum cost is a DP over the ``2^|U|`` masks, and when the hypergraph of
  the threshold has no cover, ``v`` is infeasible for *every* design at *every*
  price.  That is a bound; a search is a sample.

The robust reading -- one common drift ``delta`` applied to all sections at
once -- lives in :func:`drift_radius`, which is exact for **coverage**, and in
:func:`drift_tolerance`, which is the same question asked of the **margin**
and answers with a strictly smaller number.

**One currency.**  Every margin here is the P51 currency: the measure of the
block that a section holds, which is ``2 * min(q - L, R - q)`` when a section
``q`` sits in ``(L, R]``.  ``doubled=False`` exists because one figure (``1/33``
for the four-grid design at ``k = 4``) circulated in half units; it is a view
of the same measurement, not a second measurement.
"""

from __future__ import annotations

import itertools
from array import array
from dataclasses import dataclass
from fractions import Fraction

from .lattice import _interval_map, cell_widths, grids, mask_table, universe


def _position(U: tuple[int, ...]) -> dict[int, int]:
    return {d: i for i, d in enumerate(U)}


def cover_masks(gr: tuple[int, ...], U: tuple[int, ...],
                widths: dict[tuple[int, int], Fraction],
                threshold: Fraction) -> dict[int, int]:
    """``M -> bitmask`` of the divisors ``M`` covers at ``threshold``."""
    pos = _position(U)
    out = {}
    for M in gr:
        mask = 0
        for d in _interval_map(M):
            if d in pos and widths[(M, d)] >= threshold:
                mask |= 1 << pos[d]
        out[M] = mask
    return out


def coverage_kernel(gr: tuple[int, ...], U: tuple[int, ...],
                    widths: dict[tuple[int, int], Fraction],
                    threshold: Fraction) -> tuple[tuple[int, ...], dict[int, int]]:
    """The divisors whose fate no other divisor's fate already decides.

    Write ``cols[d]`` for the bitmask over ``gr`` of the grids whose block for
    ``d`` is at least ``threshold`` wide.  A design holds ``d`` exactly when it
    contains a grid of ``cols[d]``, so if ``cols[e]`` is contained in
    ``cols[d]`` then every design that holds ``e`` also holds ``d``: the
    constraint on ``d`` is implied, and dropping it cannot move the optimum.
    Equal columns imply each other, and the smallest divisor of that class is
    the one kept.

    The kernel is what a cheapest-design DP has to distinguish, and it is
    strictly smaller than the universe.  Measured: ``T = 2520`` at ``1/6`` has
    ``|U| = 32`` and ``|K| = 20``; ``T = 1680`` has 27 and 19; ``T = 720`` has
    22 and 14.  It shrinks with the threshold, not with the period: near the
    ceiling almost nothing is implied any more (``T = 2520`` keeps 18 of 32 at
    ``1/126`` and 23 at its ceiling), and only *above* the ceiling does every
    column empty out, which leaves one empty column that contains all the
    others.

    Returns ``(kernel, columns)`` with the columns still keyed by every divisor
    of ``U``, so a caller can read a dropped divisor's coverage.
    """
    columns: dict[int, int] = {}
    for d in U:
        mask = 0
        for i, M in enumerate(gr):
            # A pair that does not exist is not a wide block: ``widths`` is
            # keyed only by the pairs where ``d`` is a block of ``M``, and
            # reading a missing pair as zero would make every grid cover every
            # divisor at threshold ``0``.
            if (M, d) in widths and widths[(M, d)] >= threshold:
                mask |= 1 << i
        columns[d] = mask
    kernel = []
    for d in U:
        mine = columns[d]
        implied = False
        for e in U:
            if e == d:
                continue
            theirs = columns[e]
            if theirs & ~mine:
                continue
            if theirs != mine or e < d:
                implied = True
                break
        if not implied:
            kernel.append(d)
    return tuple(kernel), columns


def minimum_cost_design(gr: tuple[int, ...], U: tuple[int, ...],
                        widths: dict[tuple[int, int], Fraction],
                        threshold: Fraction,
                        period: int | None = None) -> tuple[int | None, tuple[int, ...] | None]:
    """``(cost, design)`` of the cheapest design with margin at least ``threshold``.

    ``(None, None)`` when no design reaches it, which by the duality means the
    threshold is above the ceiling of the instance.

    The state is a set of *covered kernel divisors*, not a set of covered
    divisors: :func:`coverage_kernel` proves that the dropped ones are implied,
    so the DP distinguishes ``2^|K|`` states instead of ``2^|U|`` and needs no
    table over the whole universe.  Each state carries the grids that bought
    its cost, so no layer history is kept and the walk back is a read instead
    of a reconstruction -- a predecessor chase can visit a state the DP never
    reached.

    ``period`` is carried only so that a refusal can name the instance it
    belongs to: the geometry is fully described by ``gr`` and ``U``.
    """
    kernel, columns = coverage_kernel(gr, U, widths, threshold)
    # Two arrays, not two lists: 4 bytes of cost and 8 bytes of grid mask per
    # state, measured at 20.8 B/state against 80.0 for the list version.
    pick_bytes = 8 if len(gr) <= 64 else 8 * ((len(gr) + 63) // 64)
    full = mask_table(len(kernel), "the cheapest-design DP", period=period,
                      width=4 + pick_bytes) - 1
    INF = (1 << 30) - 1
    cost = array("i", [INF]) * (full + 1)
    pick = array("Q", [0]) * (full + 1) if len(gr) <= 64 else [0] * (full + 1)
    cost[0] = 0
    covers = []
    for i, M in enumerate(gr):
        reached = 0
        for j, d in enumerate(kernel):
            if columns[d] >> i & 1:
                reached |= 1 << j
        if reached:
            covers.append((M // 2, reached, 1 << i))
    # A grid whose kernel coverage is contained in a cheaper grid's is never
    # worth buying: swapping it for that grid keeps the coverage and does not
    # raise the cost.  Measured on ``T = 2520``: 47 grids, 15 to 23 of them
    # undominated depending on the threshold.
    covers.sort(key=lambda cover: (cover[0], cover[2]))
    kept: list[tuple[int, int, int]] = []
    for cover in covers:
        if any(cover[1] & ~other[1] == 0 for other in kept):
            continue
        kept.append(cover)
    for price, reached, bit in kept:
        # Descending, so a state this grid improved is never read back as a
        # source by the same grid: the transition only ever lands on a strictly
        # larger mask, which the loop has already passed.
        for mask in range(full, -1, -1):
            c = cost[mask]
            if c >= INF:
                continue
            nxt = mask | reached
            if nxt == mask:
                continue
            c += price
            if c < cost[nxt]:
                cost[nxt] = c
                pick[nxt] = pick[mask] | bit
    total = cost[full]
    if total >= INF:
        return None, None
    chosen = tuple(M for i, M in enumerate(gr) if int(pick[full]) >> i & 1)
    return total, chosen


def margin_ladder(gr: tuple[int, ...], U: tuple[int, ...],
                  widths: dict[tuple[int, int], Fraction],
                  period: int | None = None) -> dict[int, Fraction]:
    """``cost -> best margin`` over every width in the instance.

    Monotone by construction: buying a grid never makes a design worse.  The
    reversal -- a better margin never costs fewer points -- is the one that
    fails, and that failure is the third currency.
    """
    levels = sorted({v for v in widths.values()})
    out: dict[int, Fraction] = {}
    for v in levels:
        cost, _ = minimum_cost_design(gr, U, widths, v, period)
        if cost is None:
            continue
        prev = out.get(cost)
        out[cost] = v if prev is None else max(prev, v)
    return out


def design_margin(G: tuple[int, ...], U: tuple[int, ...],
                  widths: dict[tuple[int, int], Fraction]) -> Fraction:
    """Margin guaranteed by a *set* of grids, each divisor read in its best one.

    ``min`` over the universe of ``max`` over the chosen grids.  The maximum
    and not the minimum: the grid is already paid for, so reading ``d`` in
    whichever grid seats it widest costs nothing extra.

    A design that does not seat some ``d`` at all has margin zero, not an
    error: not covering the universe is the weakest margin there is, and the
    coverage DP reads it the same way.
    """
    worst: Fraction | None = None
    for d in U:
        best = max((widths[(M, d)] for M in G if (M, d) in widths), default=None)
        if best is None:
            return Fraction(0)
        if worst is None or best < worst:
            worst = best
    return Fraction(0) if worst is None else worst


def covered_divisors(G: tuple[int, ...], Q: tuple[Fraction, ...],
                     U: tuple[int, ...]) -> tuple[int, ...]:
    """The divisors some section actually lands inside, in at least one grid.

    :func:`witness_margin` answers ``0`` for two different things -- *held, with
    no slack* and *never held at all* -- and only the first is a measurement.
    This is the reader that tells them apart, so a schedule that covers less
    than its universe cannot report the uncovered part as a margin of zero.
    """
    held = []
    for d in U:
        for M in G:
            block = _interval_map(M).get(d)
            if block is None:
                continue
            L, R = block
            if any(L < q <= R for q in Q):
                held.append(d)
                break
    return tuple(held)


def witness_margin(G: tuple[int, ...], Q: tuple[Fraction, ...],
                   U: tuple[int, ...], doubled: bool = True) -> dict[int, Fraction]:
    """Margin held by each divisor ``d`` under a shared set of sections ``Q``.

    A pair ``(M, q)`` holds ``d`` when ``q`` lands inside ``d``'s block, and
    its slack is the distance to the nearer edge.  Each ``d`` takes its own
    best pair, which is what a witness is.  ``2 * min(q - L, R - q)`` equals
    the block width when ``q`` sits at the middle, so this is the same
    currency as :func:`minimum_cost_design` -- made of sections instead of
    grids.
    """
    scale = 2 if doubled else 1
    out: dict[int, Fraction] = {}
    for d in U:
        best = Fraction(0)
        for M in G:
            block = _interval_map(M).get(d)
            if block is None:
                continue
            L, R = block
            for q in Q:
                if L < q <= R:
                    slack = min(q - L, R - q) * scale
                    if slack > best:
                        best = slack
        out[d] = best
    return out


def drift_radius(G: tuple[int, ...], Q: tuple[Fraction, ...],
                 U: tuple[int, ...], doubled: bool = True) -> dict[int, Fraction]:
    """Radius of a **common** drift ``delta`` the design still absorbs.

    A shift ``q -> q + delta`` moves every section at once, so ``d`` survives
    when some pair ``(M, q)`` ends up holding it: the admissible ``delta`` are
    the spans ``(L - q, R - q)``, and the radius is how far their union reaches
    from zero on both sides.

    Two corrections to the obvious version, both measured.  Every span is a
    candidate even when ``q`` starts outside the block -- a drift can push a
    section *into* ``d``'s block, and keeping only the pairs that already
    capture ``d`` under-reported the radius (``1/12`` against the exact
    ``7/30`` for ``(2, 3, 4, 6, 9, 12, 15, 20, 36, 45, 180)`` at ``k = 4``).
    And the reach is the component grown through overlaps, not a maximum:
    ``min(max(q - L), max(R - q))`` agrees only while every span straddles
    zero, and it stops agreeing as soon as one of them does not.

    Read the comparison with :func:`witness_margin` carefully, because it goes
    the way one expects least: the radius is at least the witness slack.  That
    is not optimism -- capture is a weaker condition than slack, and keeping
    the **margin** under drift is :func:`drift_tolerance`.
    """
    return drift_tolerance(G, Q, U, Fraction(0), doubled)


def _reach_from_zero(spans: list[tuple[Fraction, Fraction]]) -> Fraction:
    """Largest ``r`` with ``(-r, r)`` inside the union of the open ``spans``.

    Only a span that straddles zero survives both drift directions: one that
    ends at zero covers no ``delta > 0``, one that starts at zero covers no
    ``delta < 0``.  From that seed the component grows through overlaps, which
    is why a chain of two spans can reach further than either one alone --
    the reason the reach is computed and not taken as a maximum.
    """
    live = [(a, b) for a, b in spans if a < b]
    seed = [(a, b) for a, b in live if a < 0 < b]
    if not seed:
        return Fraction(0)
    lo = min(a for a, _ in seed)
    hi = max(b for _, b in seed)
    grew = True
    while grew:
        grew = False
        for a, b in live:
            if a < hi < b:
                hi, grew = b, True
            if lo < b and a < lo:
                lo, grew = a, True
    return min(-lo, hi)


def drift_tolerance(G: tuple[int, ...], Q: tuple[Fraction, ...],
                    U: tuple[int, ...], margin: Fraction = Fraction(0),
                    doubled: bool = True) -> dict[int, Fraction]:
    """Common drift every divisor absorbs **while still holding** ``margin``.

    :func:`drift_radius` answers a weaker question: it asks only that ``d``
    stays inside *some* block, so it can promise a drift while the margin
    collapses.  Shrinking every block by half the margin turns the same union
    into the condition that the margin itself survives, and then the drift has
    a price -- measured, not assumed: at the witness margin the price is zero
    in every design checked, because the binding divisor sits exactly on its
    edge.  A margin is a knife edge; the tolerance is what it costs to leave it.

    The price is linear, and then it is not: ``tolerance(v) = radius - v`` for
    every ``v`` below :func:`drift_law`'s critical margin of that divisor, and
    strictly below the line at it.  Extrapolating the linear part is how a
    design over-promises.

    ``margin`` and the result are in the same currency as
    :func:`witness_margin`; ``drift_tolerance(..., 0)`` equals
    :func:`drift_radius`.  ``doubled=False`` halves *both*, so the identity
    ``tolerance(v) = radius - v`` reads the same in either view -- a view is a
    change of units, never a change of measurement.
    """
    scale = 2 if doubled else 1
    half = Fraction(margin) / scale
    out: dict[int, Fraction] = {}
    for d in U:
        out[d] = _reach_from_zero([(a + half, b - half) for a, b in _spans(G, Q, d)]) * scale
    return out


def _largest_shrink(spans: list[tuple[Fraction, Fraction]],
                    need: Fraction) -> Fraction:
    """Largest shrink whose reach still covers ``need``, in span currency.

    :func:`drift_tolerance` maps a shrink to a reach; this inverts it.  The
    reach is linear with slope ``-1`` between two values of :func:`_shrinks`
    and jumps *down* at each of them, so the valid set is a prefix of every
    interval and the answer is read from the top down.  The endpoints are
    tested separately because a jump lands on the knot itself: a margin exactly
    at a knot can still be held when every larger one cannot.
    """
    radius = _reach_from_zero(spans)
    if need <= 0:
        return radius

    def reach_at(shrink: Fraction) -> Fraction:
        return _reach_from_zero([(a + shrink, b - shrink) for a, b in spans])

    levels = sorted({Fraction(0), radius} | {k for k in _shrinks(spans) if k < radius})
    for low, high in reversed(list(zip(levels, levels[1:]))):
        # The reach is non-increasing, so the highest interval whose floor still
        # covers ``need`` holds the supremum.  It is the line's own answer, and
        # it is a supremum when the reach jumps at ``high``: there the shrink is
        # approached, not held.
        if reach_at(low) >= need:
            return min(high, radius - need)
    return Fraction(0)


def margin_under_drift(G: tuple[int, ...], Q: tuple[Fraction, ...],
                       U: tuple[int, ...], delta: Fraction = Fraction(0),
                       doubled: bool = True) -> dict[int, Fraction]:
    """Largest margin each divisor keeps while absorbing a common drift.

    The inverse of :func:`drift_tolerance`, and the reading a service level
    actually needs: not "how much drift does this margin tolerate" but "how
    much margin survives the drift I have".  A common drift moves every section
    at once, so the two are one quantity read from either end.

    It is a **supremum**: at a break the tolerance jumps down, so the last
    margin before the jump is approached and not held.  The guarantee is strict
    -- any margin below the answer tolerates ``delta``, the answer itself may
    not.

    The design-level reading is the ``min`` over divisors, as everywhere else.
    For a design whose tolerance is ``radius - v`` all the way to zero -- every
    design the solver returns in the census of :func:`design_law` -- the answer
    is exactly ``witness_margin - delta`` and the exchange is one-for-one.  A
    design that is *not* cheapest can carry a coverage radius above its witness
    margin, and that extra reach is not spendable as margin: it breaks at a
    break below the radius, where the tolerance falls to zero instead of paying
    the line.
    """
    if delta < 0:
        raise ValueError("delta is a drift magnitude and cannot be negative")
    scale = 2 if doubled else 1
    need = Fraction(delta) / scale
    return {d: _largest_shrink(_spans(G, Q, d), need) * scale for d in U}


def guaranteed_margin(G: tuple[int, ...], Q: tuple[Fraction, ...],
                      U: tuple[int, ...], delta: Fraction = Fraction(0),
                      doubled: bool = True) -> Fraction:
    """Margin this design keeps on every divisor while absorbing ``delta``.

    One number for the whole schedule: the ``min`` of
    :func:`margin_under_drift`, and a supremum like it -- every margin strictly
    below the answer tolerates ``delta``.

    Drift is not a second budget: it is spent from the margin.  A design that
    keeps ``v`` under a drift of ``delta`` is the same design as one kept at
    ``v + delta`` with no drift, so a service level is one margin to search and
    not two budgets to trade.
    """
    return min(margin_under_drift(G, Q, U, delta, doubled).values())


def _spans(G: tuple[int, ...], Q: tuple[Fraction, ...],
           d: int) -> list[tuple[Fraction, Fraction]]:
    """Admissible drifts ``(L - q, R - q)`` for one divisor, at margin zero."""
    spans: list[tuple[Fraction, Fraction]] = []
    for M in G:
        block = _interval_map(M).get(d)
        if block is None:
            continue
        L, R = block
        for q in Q:
            spans.append((L - q, R - q))
    return spans


def _shrinks(spans: list[tuple[Fraction, Fraction]]) -> list[Fraction]:
    """Every ``h`` where the reach can change slope or jump.

    Shrinking each span by ``h`` pushes both of its ends inward, so the reach
    stays linear in ``h`` between two of these values and only meets one here.
    Three families: a span collapsing (``a``, ``b``), a span stopping to
    straddle zero (``-a``), and two ends meeting (``(b1 - a2) / 2``).
    """
    knots = {Fraction(0)}
    for a, b in spans:
        for edge in (a, b, -a):
            if edge > 0:
                knots.add(edge)
    for a1, b1 in spans:
        for a2, b2 in spans:
            for gap in (b1 - a2, b2 - a1):
                if gap > 0:
                    knots.add(gap / 2)
    return sorted(knots)


def _linear_regime(spans: list[tuple[Fraction, Fraction]], radius: Fraction) -> Fraction:
    """Largest shrink below which the reach is exactly ``radius - h``.

    Tests one interior point per interval, which is enough because both sides
    are linear there.  Beyond the returned knot the reach drops *below* the
    line and does not come back: the condition that grows a chain is strict, so
    a link never re-forms once it has broken.
    """
    knots = [k for k in _shrinks(spans) if k <= radius]
    last = Fraction(0)
    for knot in knots:
        middle = (last + knot) / 2
        if _reach_from_zero([(a + middle, b - middle) for a, b in spans]) != radius - middle:
            return last
        last = knot
    tail = (last + radius) / 2
    if _reach_from_zero([(a + tail, b - tail) for a, b in spans]) != radius - tail:
        return last
    return radius


@dataclass(frozen=True)
class DriftLaw:
    """A divisor's drift tolerance read as a function of the margin demanded.

    ``tolerance(v) = radius - v`` for every ``v < critical_margin`` and *below*
    the line at ``critical_margin`` -- so the price of a margin is one-for-one
    right up to a knife edge and a cliff after it.  ``cliff`` is how far the
    tolerance sits under the line at that edge, which is the part a linear
    extrapolation promises and the geometry does not pay.

    All three fields are in the same currency as :func:`drift_tolerance`'s
    result, so ``critical_margin`` can be handed back to it unchanged.
    """

    radius: Fraction
    critical_margin: Fraction
    cliff: Fraction


def drift_law(G: tuple[int, ...], Q: tuple[Fraction, ...],
              U: tuple[int, ...], doubled: bool = True) -> dict[int, DriftLaw]:
    """The critical margin of every divisor, and the cliff that follows it.

    Measured over the 200 even periods at ``k = 4`` and ``k = 6``: strictly
    below the critical margin the line holds in 1680 of 1680 divisor readings,
    and at it the tolerance is under the line in 369 of 369 -- never on it.
    The critical margin is a property of the span geometry, not of the budget:
    it is where a chain of overlapping spans loses a link.
    """
    scale = 2 if doubled else 1
    out: dict[int, DriftLaw] = {}
    for d in U:
        spans = _spans(G, Q, d)
        radius = _reach_from_zero(spans)
        edge = _linear_regime(spans, radius)
        here = _reach_from_zero([(a + edge, b - edge) for a, b in spans])
        out[d] = DriftLaw(radius * scale, edge * scale, (radius - edge - here) * scale)
    return out


@dataclass(frozen=True)
class DesignLaw:
    """A whole design read as a function of the margin demanded.

    ``tolerance(v) == radius - v`` on the **entire** interval
    ``0 <= v <= radius``: the design has no cliff.  The one :func:`drift_law`
    reports per divisor never survives the minimum over divisors, because the
    divisor that breaks first is not the divisor that fixes the radius.

    ``first_divisor_break`` is that earliest per-divisor break.  It is a
    pessimistic reading of the design, not a boundary of it: for 78 of the 240
    designs of the 200 even periods at ``k = 4`` and ``k = 6`` it is strictly
    smaller than the radius, and the design pays nothing at it.

    What the line rests on is one inequality per divisor: a drop of ``cliff``
    from ``radius_d`` reaches the line only while ``cliff <= radius_d - radius``
    at a critical margin inside the range, and that is the whole obstruction.

    Measured over 3076 designs -- 200 even periods, budgets
    ``{16, 24, 32, 47, 64, 96}``, ``k`` in ``{2, 3, 4, 5, 6, 7, 8}`` -- the line
    holds on every interval between consecutive knots, and the inequality holds
    with it, no exception.  It is not slack either: 1234 readings in that census
    are divisors that *do* have a cliff to jump and land exactly on the line,
    so equality is the typical case rather than a corner.  ``T = 12`` with the
    design ``(4, 6)`` is the smallest: divisor ``4`` breaks at ``1/6`` and drops
    ``1/3`` of a radius ``2/3``, which is exactly the design's excess.
    """

    radius: Fraction
    first_divisor_break: Fraction


def design_law(G: tuple[int, ...], Q: tuple[Fraction, ...],
               U: tuple[int, ...], doubled: bool = True) -> DesignLaw:
    """The one-for-one regime of a design, which is the whole of its range.

    The radius is ``min`` over divisors, and so is the tolerance at any margin;
    a minimum of lines of slope ``-1`` is a line of slope ``-1``, and a divisor
    that breaks under its own line only loses the competition when it also
    holds the radius.  Hence the design is exact where the divisors are not.
    """
    laws = drift_law(G, Q, U, doubled)
    if not laws:
        return DesignLaw(Fraction(0), Fraction(0))
    return DesignLaw(
        min(law.radius for law in laws.values()),
        min(law.critical_margin for law in laws.values()),
    )


def _probe_bands(G: tuple[int, ...], U: tuple[int, ...],
                 threshold: Fraction) -> tuple[list[tuple[int, int, Fraction, Fraction]], list[Fraction]]:
    """``(M, d, a, z)``: the sections that give ``d`` margin at least the threshold.

    ``a = L + v / 2`` and ``z = R - v / 2``: a section exactly ``v / 2`` away
    from both edges.  When ``v`` equals the block width the band collapses to a
    single point and that point is still a valid section -- discarding it
    lowers a ceiling of ``1/6`` to ``16/99``.

    The candidates are the right edges ``z``, because moving a section right
    never loses coverage.
    """
    bands, candidates = [], set()
    for d in U:
        for M in G:
            block = _interval_map(M).get(d)
            if block is None:
                continue
            L, R = block
            a, z = L + threshold / 2, R - threshold / 2
            if a <= z:
                bands.append((M, d, a, z))
                candidates.add(z)
    return bands, sorted(candidates)


def _holds(q: Fraction, L: Fraction, R: Fraction, threshold: Fraction) -> bool:
    """``q`` lands inside ``(L, R]`` holding slack at least ``threshold``.

    The left edge is open, and at ``threshold = 0`` that has to be enforced as a
    strict inequality: a section sitting exactly on ``L`` captures nothing, and
    an inclusive test there reports three sections where a brute force over
    every edge shows four are needed.  Zero slack and no coverage are not the
    same thing.  For ``threshold > 0`` the band edge ``L + threshold / 2`` is
    already inside the block, so the two tests agree.
    """
    return (q > L and q <= R
            and 2 * (q - L) >= threshold and 2 * (R - q) >= threshold)


def minimum_probes(G: tuple[int, ...], U: tuple[int, ...],
                   threshold: Fraction,
                   period: int | None = None) -> tuple[int | None, tuple[Fraction, ...] | None]:
    """Fewest shared sections that hold every divisor at ``threshold``.

    Each divisor asks to land inside one of its bands: choosing sections that
    cover all of ``U`` is set cover, solved exactly by a DP over the ``2^|U|``
    masks.  ``(None, None)`` means some divisor has no band wide enough, and
    that is infeasible for every ``k`` -- a property of the design, not a
    budget.

    ``period`` is carried only so that a refusal can name the instance: the
    geometry is fully described by ``G`` and ``U``.
    """
    bands, candidates = _probe_bands(G, U, threshold)
    if len({d for _, d, _, _ in bands}) < len(U):
        return None, None
    # Two tables of ``2^|U|`` -- the cost of a cover and the step it came from.
    full = mask_table(len(U), "the set-cover DP", period=period, width=2) - 1
    pos = _position(U)
    blocks = {(M, d): _interval_map(M)[d] for M, d, _, _ in bands}
    # Both edges of both objects: the block edges where a section can gain a
    # divisor, and the band edges where it can gain slack.  Keeping only the
    # band right edges was enough while the closed left edge let an unusable
    # point masquerade as a candidate; once that was fixed the restricted set
    # gave five sections for a universe a brute force covers with four.
    every_edge = {e for (_, _), (L, R) in blocks.items() for e in (L, R)}
    every_edge |= {a for _, _, a, _ in bands} | set(candidates)
    candidates = sorted(every_edge)
    cover = []
    for q in candidates:
        mask = 0
        for M, d, _, _ in bands:
            L, R = blocks[(M, d)]
            if _holds(q, L, R, threshold):
                mask |= 1 << pos[d]
        cover.append(mask)
    INF = len(U) + 1
    dp = [INF] * (full + 1)
    parent: list[tuple[int, int] | None] = [None] * (full + 1)
    dp[0] = 0
    for mask in range(full + 1):
        c = dp[mask]
        if c >= INF:
            continue
        free = (mask ^ full) & -(mask ^ full)
        for i, m in enumerate(cover):
            if not m & free:
                continue
            nxt = mask | m
            if c + 1 < dp[nxt]:
                dp[nxt], parent[nxt] = c + 1, (mask, i)
    if dp[full] >= INF:
        return None, None
    Q, mask = [], full
    while mask:
        step = parent[mask]
        assert step is not None
        prev, i = step
        Q.append(candidates[i])
        mask = prev
    return dp[full], tuple(sorted(Q))


def best_probes(G: tuple[int, ...], U: tuple[int, ...], k: int,
                period: int | None = None) -> tuple[Fraction, tuple[Fraction, ...]]:
    """Best margin reachable with at most ``k`` shared sections, exact.

    A vertex argument bounds the search: once an assignment is fixed, the
    maximum of a minimum of linear forms sits at a vertex, and at a vertex two
    constraints are equal, so the value is the difference of two block edges.
    ``minimum_probes`` is monotone on that ordered list, so it is searched by
    bisection instead of scanned.
    """
    edges = sorted({e for M in G for L, R in _interval_map(M).values() for e in (L, R)})
    levels = sorted({b - a for a in edges for b in edges if b > a})
    base, Q0 = minimum_probes(G, U, Fraction(0), period=period)
    best, witness = (Fraction(0), Q0) if base is not None and base <= k else (None, ())
    lo, hi = 0, len(levels) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        n, Q = minimum_probes(G, U, levels[mid], period=period)
        if n is not None and n <= k:
            best, witness = levels[mid], Q
            lo = mid + 1
        else:
            hi = mid - 1
    if best is None:
        best = Fraction(0)
    return best, witness


def probe_ladder(G: tuple[int, ...], U: tuple[int, ...],
                 ks: tuple[int, ...],
                 period: int | None = None) -> dict[int, tuple[Fraction, tuple[Fraction, ...]]]:
    """``k -> (best margin, witness sections)`` for each ``k`` in ``ks``."""
    return {k: best_probes(G, U, k, period=period) for k in ks}


def coverage_census(period: int) -> dict[str, int]:
    """Count over all ``2^|grids|`` designs: how many cover, how many are cheap.

    A flat census, not a search: it exists so that a number quoted in prose
    can be reproduced by a command instead of a memory.  Its inner loop reads
    every grid for every design, so the reserve is ``|grids| * 2^|grids|``
    cells -- the widest per-mask cost in the package, and the reason a census
    of ``720`` is refused while ``360`` runs.
    """
    gr = grids(period)
    U = universe(period)
    widths = cell_widths(gr, U)
    cover = cover_masks(gr, U, widths, Fraction(0))
    masks = [cover[M] for M in gr]
    costs = [M // 2 for M in gr]
    designs = mask_table(len(gr), "the design census", period=period,
                         width=len(gr))
    full = (1 << len(U)) - 1
    covering = 0
    cheapest = None
    for bits in range(designs):
        reached, cost = 0, 0
        for i in range(len(gr)):
            if bits >> i & 1:
                reached |= masks[i]
                cost += costs[i]
        if reached != full:
            continue
        covering += 1
        if cheapest is None or cost < cheapest:
            cheapest = cost
    return {"designs": designs, "covering": covering, "cheapest": cheapest or 0}


def all_designs(gr: tuple[int, ...]) -> "itertools.chain":
    """Every design of the grids, by increasing size.  For controls only."""
    return itertools.chain.from_iterable(
        itertools.combinations(gr, r) for r in range(1, len(gr) + 1))
