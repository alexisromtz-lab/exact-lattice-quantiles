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
from fractions import Fraction

from .lattice import _interval_map, cell_widths, grids, universe


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


def minimum_cost_design(gr: tuple[int, ...], U: tuple[int, ...],
                        widths: dict[tuple[int, int], Fraction],
                        threshold: Fraction) -> tuple[int | None, tuple[int, ...] | None]:
    """``(cost, design)`` of the cheapest design with margin at least ``threshold``.

    ``(None, None)`` when no design reaches it, which by the duality means the
    threshold is above the ceiling of the instance.  The DP keeps the parent
    it travelled through; walking the path back by clearing bits takes a
    predecessor the DP never visited.
    """
    n = len(U)
    full = (1 << n) - 1
    cover = cover_masks(gr, U, widths, threshold)
    INF = 10**9
    dp = [[INF] * (1 << n) for _ in range(len(gr) + 1)]
    parent: list[dict[int, int]] = [dict() for _ in range(len(gr) + 1)]
    dp[0][0] = 0
    for i, M in enumerate(gr):
        prev, nxt, par = dp[i], dp[i + 1], parent[i + 1]
        for mask in range(1 << n):
            c = prev[mask]
            if c >= INF:
                continue
            if c < nxt[mask]:
                nxt[mask] = c
                par[mask] = mask
            m2 = mask | cover[M]
            if c + M // 2 < nxt[m2]:
                nxt[m2] = c + M // 2
                par[m2] = mask
    total = dp[len(gr)][full]
    if total >= INF:
        return None, None
    chosen, mask = [], full
    for i in range(len(gr), 0, -1):
        step = parent[i].get(mask)
        if step is None:
            break
        if step != mask:
            chosen.append(gr[i - 1])
        mask = step
    return total, tuple(sorted(chosen))


def margin_ladder(gr: tuple[int, ...], U: tuple[int, ...],
                  widths: dict[tuple[int, int], Fraction]) -> dict[int, Fraction]:
    """``cost -> best margin`` over every width in the instance.

    Monotone by construction: buying a grid never makes a design worse.  The
    reversal -- a better margin never costs fewer points -- is the one that
    fails, and that failure is the third currency.
    """
    levels = sorted({v for v in widths.values()})
    out: dict[int, Fraction] = {}
    for v in levels:
        cost, _ = minimum_cost_design(gr, U, widths, v)
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
    capture ``d`` under-reported the radius (``1/6`` against the exact
    ``37/210`` for ``(2, 3, 4, 6, 9, 12, 15, 20, 36, 45, 180)`` at ``k = 4``).
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

    ``margin`` and the result are in the same currency as
    :func:`witness_margin`; ``drift_tolerance(..., 0)`` equals
    :func:`drift_radius`.
    """
    scale = 2 if doubled else 1
    half = Fraction(margin) / 2
    out: dict[int, Fraction] = {}
    for d in U:
        spans: list[tuple[Fraction, Fraction]] = []
        for M in G:
            block = _interval_map(M).get(d)
            if block is None:
                continue
            L, R = block
            for q in Q:
                spans.append((L - q + half, R - q - half))
        out[d] = _reach_from_zero(spans) * scale
    return out


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


def minimum_probes(G: tuple[int, ...], U: tuple[int, ...],
                   threshold: Fraction) -> tuple[int | None, tuple[Fraction, ...] | None]:
    """Fewest shared sections that hold every divisor at ``threshold``.

    Each divisor asks to land inside one of its bands: choosing sections that
    cover all of ``U`` is set cover, solved exactly by a DP over the ``2^|U|``
    masks.  ``(None, None)`` means some divisor has no band wide enough, and
    that is infeasible for every ``k`` -- a property of the design, not a
    budget.
    """
    bands, candidates = _probe_bands(G, U, threshold)
    if len({d for _, d, _, _ in bands}) < len(U):
        return None, None
    pos = _position(U)
    cover = []
    for q in candidates:
        mask = 0
        for _, d, a, z in bands:
            if a <= q <= z:
                mask |= 1 << pos[d]
        cover.append(mask)
    full = (1 << len(U)) - 1
    INF = len(U) + 1
    dp = [INF] * (1 << len(U))
    parent: list[tuple[int, int] | None] = [None] * (1 << len(U))
    dp[0] = 0
    for mask in range(1 << len(U)):
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


def best_probes(G: tuple[int, ...], U: tuple[int, ...], k: int) -> tuple[Fraction, tuple[Fraction, ...]]:
    """Best margin reachable with at most ``k`` shared sections, exact.

    A vertex argument bounds the search: once an assignment is fixed, the
    maximum of a minimum of linear forms sits at a vertex, and at a vertex two
    constraints are equal, so the value is the difference of two block edges.
    ``minimum_probes`` is monotone on that ordered list, so it is searched by
    bisection instead of scanned.
    """
    edges = sorted({e for M in G for L, R in _interval_map(M).values() for e in (L, R)})
    levels = sorted({b - a for a in edges for b in edges if b > a})
    base, Q0 = minimum_probes(G, U, Fraction(0))
    best, witness = (Fraction(0), Q0) if base is not None and base <= k else (None, ())
    lo, hi = 0, len(levels) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        n, Q = minimum_probes(G, U, levels[mid])
        if n is not None and n <= k:
            best, witness = levels[mid], Q
            lo = mid + 1
        else:
            hi = mid - 1
    if best is None:
        best = Fraction(0)
    return best, witness


def probe_ladder(G: tuple[int, ...], U: tuple[int, ...],
                 ks: tuple[int, ...]) -> dict[int, tuple[Fraction, tuple[Fraction, ...]]]:
    """``k -> (best margin, witness sections)`` for each ``k`` in ``ks``."""
    return {k: best_probes(G, U, k) for k in ks}


def coverage_census(period: int) -> dict[str, int]:
    """Count over all ``2^|grids|`` designs: how many cover, how many are cheap.

    A flat census, not a search: it exists so that a number quoted in prose
    can be reproduced by a command instead of a memory.
    """
    gr = grids(period)
    U = universe(period)
    widths = cell_widths(gr, U)
    cover = cover_masks(gr, U, widths, Fraction(0))
    masks = [cover[M] for M in gr]
    costs = [M // 2 for M in gr]
    full = (1 << len(U)) - 1
    covering = 0
    cheapest = None
    for bits in range(1 << len(gr)):
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
    return {"designs": 1 << len(gr), "covering": covering, "cheapest": cheapest or 0}


def all_designs(gr: tuple[int, ...]) -> "itertools.chain":
    """Every design of the grids, by increasing size.  For controls only."""
    return itertools.chain.from_iterable(
        itertools.combinations(gr, r) for r in range(1, len(gr) + 1))
