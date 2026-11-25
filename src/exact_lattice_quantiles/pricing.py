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


def price_form_census(period: int, reduced: bool = True) -> PriceFormCensus:
    """Census the closed form against the exact minimum, over every request.

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
