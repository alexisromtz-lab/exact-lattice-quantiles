"""Grids, seats, carriers and the universe ``U(T)``.

The objects here are the geometry the rest of the package reasons about.

* A **grid** ``M`` is a divisor of the period above one.  Grid ``M`` cuts the
  phase into ``N = M // 2`` section *classes*.
* A **seat** is the block of classes that a divisor occupies.  Divisor ``d``
  sits in grid ``M`` when ``M`` is a multiple of ``d``; the classes are
  ``j`` with ``1 <= j <= M // 2`` and ``M / gcd(j, M) == d``.
* A **carrier** of ``M`` is the divisor whose seat holds the median class: the
  first one where the accumulated mass crosses the middle.
* The **universe** ``U(T)`` is the set of carriers of the even grids of ``T``.

Grids may be even or odd.  Only even grids contribute to ``U(T)``, because a
section is a position in a phase of ``M / 2`` classes and an odd ``M`` has no
such middle.
"""

from __future__ import annotations

import math
from fractions import Fraction
from functools import lru_cache

from .arithmetic import divisors, mass


@lru_cache(maxsize=None)
def grids(period: int) -> tuple[int, ...]:
    """Every grid of ``period``: its divisors above one, ascending.

    This is the whole universe of grids, not a list of candidates.  The
    minimum the package takes over grids is a minimum over a universe.
    """
    return tuple(M for M in divisors(period) if M >= 2)


@lru_cache(maxsize=None)
def _seat_map(M: int) -> dict[int, tuple[int, int]]:
    """``d -> (lo, hi)``, the block of classes divisor ``d`` occupies in ``M``."""
    if M < 2:
        raise ValueError(f"a grid below 2 has no classes, got {M}")
    count: dict[int, int] = {}
    for j in range(1, M // 2 + 1):
        d = M // math.gcd(j, M)
        count[d] = count.get(d, 0) + 1
    total, seats = 0, {}
    for d in sorted(count, reverse=True):
        seats[d] = (total + 1, total + count[d])
        total += count[d]
    return seats


def seats(M: int) -> tuple[tuple[int, int, int], ...]:
    """``(d, lo, hi)`` for every divisor ``d`` of the grid ``M``, ascending."""
    return tuple((d, lo, hi) for d, (lo, hi) in sorted(_seat_map(M).items()))


@lru_cache(maxsize=None)
def _interval_map(M: int) -> dict[int, tuple[Fraction, Fraction]]:
    """``d -> (L, R]``, the same block on the space of sections.

    The section ``q`` captures ``d`` when ``q`` lands in ``(lo - 1) / N, hi / N]``
    with ``N = M // 2``.  The left edge is open and the right edge is closed;
    opening both moved the margin of the thin divisors by a whole class.
    """
    N = M // 2
    return {d: (Fraction(lo - 1, N), Fraction(hi, N))
            for d, (lo, hi) in _seat_map(M).items()}


def cell_intervals(M: int) -> tuple[tuple[int, Fraction, Fraction], ...]:
    """``(d, L, R)`` for every divisor of ``M``, ascending by divisor."""
    return tuple((d, L, R) for d, (L, R) in sorted(_interval_map(M).items()))


def cell_width(M: int, d: int) -> Fraction:
    """The measure of the sections that capture ``d`` in grid ``M``.

    It is ``2a / M`` and not ``a / M``: the section space is the space of
    classes, so a block of ``a`` classes out of ``N = M // 2`` measures
    ``a / N = 2a / M``.  Dividing by ``M`` halves every margin in the package.
    """
    lo, hi = _seat_map(M)[d]
    return Fraction(2 * (hi - lo + 1), M)


@lru_cache(maxsize=None)
def carriers(M: int) -> tuple[int, ...]:
    """The carrier (or carriers) of grid ``M``, read off the mass.

    The seat of ``d`` in ``M`` holds ``mass(d) / 2`` classes, so the median is
    the first divisor where the accumulated mass crosses the middle:
    ``M // 2`` when four divides ``M``, ``M // 2 + 1`` when it does not.

    The ``+ 1`` in the second branch is redundant and stays for legibility
    against the seat-by-seat definition it replaces: the accumulated mass is
    always even -- ``mass(2) = 2`` and ``phi(d)`` is even for ``d > 2`` --
    while ``M // 2`` is odd when ``M = 2 (mod 4)``, so both thresholds select
    the same divisor.

    A tie returns two divisors.  Four dividing ``M`` and the mass landing
    exactly on a quarter is the case where two blocks straddle the median.
    """
    ds = divisors(M)[1:]
    target = M // 2 if M % 4 == 0 else M // 2 + 1
    accumulated = 0
    carrier, half = ds[0], 0
    for d in ds:
        accumulated += mass(d)
        if accumulated >= target:
            carrier, half = d, accumulated // 2
            break
    if M % 4 == 0 and half == M // 4:
        return (carrier, min(d for d in ds if d > carrier))
    return (carrier,)


@lru_cache(maxsize=None)
def universe_from_factorization(period: int) -> tuple[int, ...]:
    """``U(T)`` read from the factorization, without building any grid.

    ``U(T)`` is the set of carriers of the even grids of ``T``.  Both the grids
    and their carriers follow from the primes of ``T``, so the whole universe
    is available before a single seat is computed.  ``U(T)`` is never empty:
    the grid ``T`` itself is even for every even ``T``.
    """
    found: set[int] = set()
    for M in divisors(period):
        if M >= 2 and M % 2 == 0:
            found.update(carriers(M))
    return tuple(sorted(found))


def universe(period: int) -> tuple[int, ...]:
    """``U(T)``, by the same route as :func:`universe_from_factorization`.

    Kept as a separate name so the two routes can be compared: an independent
    count over the grids is what makes the factorization route falsifiable.
    """
    found: set[int] = set()
    for M in grids(period):
        if M % 2 == 0:
            found.update(carriers(M))
    return tuple(sorted(found))


@lru_cache(maxsize=None)
def cell_widths(period_grids: tuple[int, ...], U: tuple[int, ...]) -> dict[tuple[int, int], Fraction]:
    """``(M, d) -> width`` for every pair with ``d`` in ``U``.

    ``period_grids`` and ``U`` are tuples so the result can be cached; call it
    through :func:`widths` with a period and it will do the rest.
    """
    keep = set(U)
    return {(M, d): cell_width(M, d)
            for M in period_grids for d in _seat_map(M) if d in keep}


def widths(period: int) -> dict[tuple[int, int], Fraction]:
    """``(M, d) -> width`` over the grids and universe of ``period``."""
    return dict(cell_widths(grids(period), universe(period)))
