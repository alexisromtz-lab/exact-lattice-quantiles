"""What a request costs: an exact minimum, and the closed form that is only local.

The tempting formula reads: cover a request by buying one grid per *maximal*
divisor, so the price is ``sum(d // 2 for d in maximal(needed))``.  It is a
legal design -- grid ``d`` always covers ``d`` -- so it is an **upper bound**,
and over the reduced universe ``U(T)`` of ``T = 180`` it is the exact minimum
for all ``4 095`` non-empty requests (``47`` points, zero discrepancy).

It is not the price.  Over the seventeen divisors of ``180`` above one it
overcharges ``9 984`` of the ``131 071`` non-empty requests, because a block of
grids can be bought whole for less than the sum of its parts: ``{12, 20, 30}``
costs ``31`` bought separately and ``30`` bought as the single grid ``60``,
their least common multiple.  Two incomparable divisors never save
(``lcm(a, b) >= a + b``), so the formula is exact on every request whose every
frontier has no triple, and wrong everywhere else.

So the export is two functions and one census: :func:`price_exact` is the
minimum (a subset DP, no formula), :func:`price_one_per_maximal` is the closed
form, and :func:`price_form_census` states over how many requests of an
instance the two agree -- separately for both universes, because the gap
between them is the finding.  A closed form whose agreement is not censused is
a guess with a proof on one side.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import gcd

from .lattice import cell_widths, grids, universe

_INF = 10**9


def maximal_divisors(needed: tuple[int, ...]) -> tuple[int, ...]:
    """The elements of ``needed`` that no other element is a multiple of."""
    return tuple(d for d in needed
                 if not any(e != d and e % d == 0 for e in needed))


def price_one_per_maximal(needed: tuple[int, ...]) -> int:
    """The closed form: buy one grid per maximal divisor.

    Always a feasible price for ``needed``, never below the minimum.
    """
    return sum(d // 2 for d in maximal_divisors(needed))


def _maximal_table(U: tuple[int, ...]) -> list[int]:
    """``mask -> bitmask`` of its maximal positions, built one bit at a time.

    Adding ``d`` to a set drops the maximals it divides and keeps ``d`` itself
    only when no surviving maximal is a multiple of it.
    """
    n = len(U)
    table = [0] * (1 << n)
    for mask in range(1, 1 << n):
        low = mask & -mask
        i = low.bit_length() - 1
        kept = table[mask ^ low]
        rest = kept
        while rest:
            bit = rest & -rest
            j = bit.bit_length() - 1
            if U[i] % U[j] == 0:
                kept &= ~bit
            rest ^= bit
        dominated = False
        rest = kept
        while rest:
            bit = rest & -rest
            j = bit.bit_length() - 1
            if U[j] % U[i] == 0:
                dominated = True
                break
            rest ^= bit
        table[mask] = kept if dominated else kept | low
    return table


@lru_cache(maxsize=None)
def _cheapest(period: int, reduced: bool
              ) -> tuple[tuple[int, ...], tuple[int, ...], tuple[tuple[int, int], ...]]:
    """``(U, best, covers)`` where ``best[mask]`` is the cheapest cover of a superset.

    ``covers`` is ``((grid, mask), ...)`` at threshold zero: a grid covers
    exactly the divisors it is a multiple of.
    """
    gr = grids(period)
    U = universe(period) if reduced else gr
    widths = cell_widths(gr, U)
    pos = {d: i for i, d in enumerate(U)}
    covers = []
    for M in gr:
        mask = 0
        for d in U:
            if widths.get((M, d), Fraction(0)) > Fraction(0):
                mask |= 1 << pos[d]
        covers.append((M, mask))

    size = 1 << len(U)
    f = [_INF] * size
    f[0] = 0
    for M, mask in covers:
        nxt = f[:]
        for s in range(size):
            v = f[s]
            if v < _INF and v + M // 2 < nxt[s | mask]:
                nxt[s | mask] = v + M // 2
        f = nxt
    # A request is met by any cover of a superset of it.
    for i in range(len(U)):
        bit = 1 << i
        for s in range(size):
            if not s & bit and f[s | bit] < f[s]:
                f[s] = f[s | bit]
    return U, tuple(f), tuple(covers)


def price_exact(period: int, needed: tuple[int, ...], reduced: bool = True) -> int:
    """Cheapest design covering every element of ``needed``.

    The minimum over designs, by dynamic programming over the ``2^|U|``
    coverage masks.  No closed form.
    """
    U, best, _ = _cheapest(period, reduced)
    mask = 0
    for d in needed:
        mask |= 1 << U.index(d)
    cost = best[mask]
    if cost >= _INF:  # pragma: no cover - every divisor is covered by its grid
        raise ValueError(f"request {needed} is not coverable in period {period}")
    return cost


@dataclass(frozen=True)
class PriceFormCensus:
    """How often the closed form is the minimum, on one universe of one instance."""

    period: int
    reduced: bool
    universe_size: int
    consumers: int
    discrepancies: int
    frontiers: int
    saving_frontiers: int
    cheapest_cover: int
    one_grid: int

    @property
    def share(self) -> Fraction:
        return Fraction(self.discrepancies, self.consumers)

    @property
    def adequate(self) -> bool:
        """True when the closed form is the price on every request."""
        return self.discrepancies == 0

    def as_dict(self) -> dict[str, object]:
        return {
            "period": self.period,
            "universe": "reduced" if self.reduced else "divisors",
            "universe_size": self.universe_size,
            "consumers": self.consumers,
            "discrepancies": self.discrepancies,
            "share": str(self.share),
            "frontiers": self.frontiers,
            "saving_frontiers": self.saving_frontiers,
            "cheapest_cover": self.cheapest_cover,
            "one_grid": self.one_grid,
            "adequate": self.adequate,
        }


def price_form_sweep(low: int, high: int, reduced: bool = True,
                     max_masks: int = 1 << 22) -> dict[str, object]:
    """The census over every even period of a window, skipping the unaffordable.

    ``skipped`` is published rather than swallowed: a sweep that silently drops
    the instances it cannot afford reports a rate over a population it did not
    declare.  Measured over ``2..400`` on ``U(T)``: 200 even periods, none
    skipped, the closed form exact in all 200, and cheaper than a single grid in
    192 of them.
    """
    periods = adequate = cheaper = 0
    skipped: list[int] = []
    costliest = (0, 0)
    for period in range(low, high + 1):
        if period % 2:
            continue
        n = len(universe(period) if reduced else grids(period))
        if 1 << n > max_masks:
            skipped.append(period)
            continue
        row = price_form_census(period, reduced=reduced)
        periods += 1
        if row.adequate:
            adequate += 1
        if row.cheapest_cover < row.one_grid:
            cheaper += 1
        if row.cheapest_cover > costliest[0]:
            costliest = (row.cheapest_cover, period)
    return {
        "low": low,
        "high": high,
        "universe": "reduced" if reduced else "divisors",
        "periods": periods,
        "adequate": adequate,
        "inadequate": periods - adequate,
        "cheaper_than_one_grid": cheaper,
        "costliest": costliest,
        "skipped": skipped,
    }


def price_form_census(period: int, reduced: bool = True) -> PriceFormCensus:
    """Census the closed form against the exact minimum, over every request.

    The cost is ``O(2^|U| * |grids|)``, and the size of ``U`` is what decides
    whether the question is affordable at all: over ``U(T)`` the masks are
    ``2^12`` for ``180`` and ``2^17`` for ``360``; over every divisor they are
    ``2^17`` for ``180`` and ``2^23`` for ``360``, which is where this
    implementation stops.  A census that cannot be run is not a census: the
    universes are named so that the caller knows which one was measured.

    A *frontier* is the set of maximal divisors of a request.  Requests sharing
    a frontier share the closed-form price, so ``frontiers`` counts how many
    distinct answers the formula can give, and a frontier *saves* when the
    exact price of that very frontier is below its closed form.
    """
    U, best, covers = _cheapest(period, reduced)
    n = len(U)
    table = _maximal_table(U)
    form_of: dict[int, int] = {}
    seen: set[int] = set()
    saving: set[int] = set()
    consumers = discrepancies = 0
    for mask in range(1, 1 << n):
        exact = best[mask]
        maxmask = table[mask]
        form = form_of.get(maxmask)
        if form is None:
            form = sum(U[i] // 2 for i in range(n) if maxmask >> i & 1)
            form_of[maxmask] = form
        consumers += 1
        if exact != form:
            discrepancies += 1
            saving.add(maxmask)
        seen.add(maxmask)
    full = (1 << n) - 1
    one_grid = min((M // 2 for M, mask in covers if mask == full), default=0)
    return PriceFormCensus(
        period=period,
        reduced=reduced,
        universe_size=n,
        consumers=consumers,
        discrepancies=discrepancies,
        frontiers=len(seen),
        saving_frontiers=len(saving),
        cheapest_cover=best[full],
        one_grid=one_grid,
    )


def _lcm(a: int, b: int) -> int:
    return a * b // gcd(a, b)


def best_block_saving(needed: tuple[int, ...]) -> tuple[int, tuple[int, ...], int, int]:
    """Largest saving any block inside ``needed`` can make, by exact search.

    Returns ``(saving, witness, bound, cover)``: ``saving`` is how many points
    the closed form overcharges on ``witness``, ``bound`` is that closed-form
    price and ``cover`` is the price of the single grid ``lcm(witness)``.  A
    saving of zero means no block is cheaper than its parts.

    The witness is an **antichain**.  The naive predicate -- "any subset whose
    least common multiple costs less than its parts" -- is wrong: if ``b``
    divides ``b'`` and both are in the subset, the grid ``b'`` already covers
    ``b``, so ``b // 2`` is charged twice.  Over the seventeen divisors of
    ``180`` that predicate reports a deficit of ``180`` points with the whole
    universe as its witness; restricted to antichains it reports ``3`` with
    ``{36, 60, 90}``, which is the witness the exact DP also finds, and
    ``180 // 2 = 90`` is the price of the one grid that covers it.

    Three is the smallest block worth searching: two divisors never save, since
    ``lcm(a, b) >= a + b`` whenever both are at least two, so a pair always
    costs at least as much merged as it costs apart.
    """
    N = tuple(sorted(set(needed)))
    n = len(N)
    suffix = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        suffix[i] = suffix[i + 1] + N[i] // 2
    best = [0, (), 0, 0]

    def walk(start: int, chosen: list[int], merged: int, total: int) -> None:
        if len(chosen) >= 3 and total - merged // 2 > best[0]:
            best[0], best[1] = total - merged // 2, tuple(chosen)
            best[2], best[3] = total, merged // 2
        # Nothing reachable from here can pay less than the merged grid.
        if start >= n or total + suffix[start] <= merged // 2:
            return
        for j in range(start, n):
            d = N[j]
            if any(d % c == 0 or c % d == 0 for c in chosen):
                continue
            merged_with = _lcm(merged, d)
            if total + d // 2 + suffix[j + 1] <= merged_with // 2:
                continue
            chosen.append(d)
            walk(j + 1, chosen, merged_with, total + d // 2)
            chosen.pop()

    walk(0, [], 1, 0)
    return best[0], best[1], best[2], best[3]


def exactness_frontier(max_period: int = 8000) -> dict:
    """The first even period whose reduced universe the closed form misprices.

    The census this package published says the closed form *is* the minimum
    over ``U(T)`` for all 200 even periods to 400.  True, and bounded: the
    first period where it fails is ``1080``, where ``{54, 90, 135}`` costs
    ``139`` bought separately and ``135`` as the single grid ``270``.  A census
    without its boundary reads as a law, and the boundary is not where a static
    reading predicts it -- ``2520`` is the same block with ``m = 9`` at more
    than twice the period.

    The witness is certified by the cover, not by a DP: ``{lcm(witness)}`` is a
    legal design that covers the block, so its price is an upper bound the DP
    would only confirm.  At ``1080`` the DP is ``2^23`` masks and at ``2520``
    it is ``2^32``; the certificate has to survive where the search cannot run.
    """
    breaks: list[tuple[int, int, tuple[int, ...], int, int]] = []
    for period in range(2, max_period + 1, 2):
        saving, witness, bound, cover = best_block_saving(universe(period))
        if saving > 0:
            breaks.append((period, saving, witness, bound, cover))
    scanned = max_period // 2
    if not breaks:
        return {"period": None, "saving": 0, "witness": (), "bound": 0,
                "cover": 0, "grid": 0, "breaks": 0, "scanned": scanned}
    period, saving, witness, bound, cover = breaks[0]
    return {
        "period": period,
        "saving": saving,
        "witness": witness,
        "bound": bound,
        "cover": cover,
        "grid": _lcm_grid(witness),
        "breaks": len(breaks),
        "scanned": scanned,
    }


def _lcm_grid(witness: tuple[int, ...]) -> int:
    merged = 1
    for d in witness:
        merged = _lcm(merged, d)
    return merged
