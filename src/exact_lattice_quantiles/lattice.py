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

The tables that answer these questions are indexed by bitmask, so their size is
``2^n`` in the number of divisors and it is the one limit of the implementation
that is not a property of the instance.  :data:`TABLE_BUDGET` names it and
:func:`mask_table` refuses over it by name, before the allocation, instead of
letting the kernel kill the process with no verdict at all.
"""

from __future__ import annotations

import math
from fractions import Fraction
from functools import lru_cache

from .arithmetic import divisors, mass

TABLE_BUDGET = 1 << 29
"""Cells one call may allocate: ``2^29`` slots, at most a few GB resident.

A *cell* is one slot of one table row, so a DP keeping ``rows`` tables of
``2^n`` masks asks for ``rows * 2^n`` cells.  The number is measured, not
chosen: the largest instance answered today is ``T = 1260`` at ``36 * 2^23``
cells (302 million, 1.46 GB peak RSS), and the next class up is ``T = 1800``
at ``35 * 2^24`` (604 million), which was left running past every wall clock
this machine allows.  ``2^29`` sits between them, so every instance that has
ever returned still returns and every one that thrashed is refused by name.
"""


class MaskTableTooLarge(ValueError):
    """The answer exists; this implementation will not allocate its table.

    Not a property of the instance and not a budget the caller set: the
    schedule is there, and the table that would find it is above
    :data:`TABLE_BUDGET`.  Distinct from :class:`OddPeriodUnsupported` for that
    reason -- an odd period has no answer to give, this one has an answer and
    no room to compute it.

    A ``ValueError`` so that a caller already catching the package's refusals
    keeps working; ``rc = 5`` so that a shell can tell the two apart.
    """

    rc = 5

    def __init__(self, site: str, count: int, width: int, period: int | None) -> None:
        self.site = site
        self.count = count
        self.width = width
        self.period = period
        self.cells = width << count
        where = f"period {period}: " if period is not None else ""
        rows = "row" if width == 1 else "rows"
        super().__init__(
            f"{where}{site} would allocate 2^{count} masks over {width} {rows} = "
            f"{self.cells} cells, above the budget of {TABLE_BUDGET}; the answer "
            f"exists and this implementation will not allocate the table")


def mask_table(count: int, site: str, period: int | None = None,
               width: int = 1) -> int:
    """``1 << count``, the size of a mask table, if it fits :data:`TABLE_BUDGET`.

    ``width`` is how many slots the loop body touches for each mask -- one for a
    table read once, the number of rows for a DP that keeps a table per row --
    so the reserve asked for is the work the call is about to do and not only
    the memory it will hold.  Refuses with :class:`MaskTableTooLarge`, naming
    the arithmetic, before a single list is built.
    """
    if width < 1:
        raise ValueError(f"a table has at least one row, got {width}")
    if width << count > TABLE_BUDGET:
        raise MaskTableTooLarge(site, count, width, period)
    return 1 << count


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
    """``d -> (L, R]``, the same block on the period, in the one currency.

    The section ``q`` captures ``d`` when ``q`` lands in
    ``(2 * (lo - 1) / M, 2 * hi / M]``, so a block of ``a`` classes measures
    ``2a / M`` -- the same number :func:`cell_width` reports, which is the
    currency every margin in the package is quoted in.

    Dividing by ``N = M // 2`` instead is the tempting error: the blocks then
    tile the section circle to exactly ``1`` for every grid, but each odd grid
    is stretched by ``M / (M - 1)``, and a witness on an odd grid comes out
    wider than the design that contains it.  The ``1 / M`` the period currency
    leaves uncovered on an odd grid is the carrier point -- the phase that is
    exactly zero -- which no section occupies; hiding it under a rescaling is
    what made the two currencies disagree.

    The left edge is open and the right edge is closed; opening both moved the
    margin of the thin divisors by a whole class.
    """
    return {d: (Fraction(2 * (lo - 1), M), Fraction(2 * hi, M))
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


def class_values(grid: int) -> tuple[int, ...]:
    """The divisor each class of ``grid`` meets, in class order.

    Class ``j`` meets ``grid // gcd(j, grid)``; that is the whole geometry, and
    it is available without a seat table, without :func:`mass` and without a
    totient.  Sorted, the tuple's middle is the carrier -- see
    :func:`median_carriers`.
    """
    if grid < 2:
        raise ValueError(f"a grid below 2 has no classes, got {grid}")
    return tuple(grid // math.gcd(j, grid) for j in range(1, grid // 2 + 1))


def median_carriers(grid: int) -> tuple[int, ...]:
    """The carriers of ``grid`` as the median of its class values.

    The third route to the same answer, and the one that shares the least: the
    mass route accumulates ``euler_phi`` over the divisors, the factorization
    route walks the primes, and this one sorts the classes and takes the middle
    -- it depends on ``gcd`` and on nothing else in the package.  A bias in
    ``mass``, in ``divisors`` or in the seat ordering would travel with both of
    the other routes and not with this one, so the three together decide
    whether ``U(T)`` is a property of the lattice or of the way it is read.
    The phase count ``grid // gcd(j, grid)`` is the referee the width contract
    already trusts -- ``test_solver_contract.py``, the currency test -- so this
    is that same referee pointed at the carrier instead of at the block.

    The tie stops being a branch.  A grid with ``4 | grid`` has an even number
    of classes, so its median sits *between* two of them, and the carriers are
    the **distinct** values at those two positions: two when they differ and
    one when they do not, which is exactly what :func:`carriers` spells out.

    An odd grid is refused by name: a section is a position in a phase of
    ``grid // 2`` classes and an odd grid has no middle one.  The mass route
    accepts odd grids because the mass never asks for a middle; this one does,
    so it stops where the geometry stops instead of guessing a convention.
    """
    if grid < 2:
        raise ValueError(f"a grid below 2 has no classes, got {grid}")
    if grid % 2:
        raise ValueError(f"the median of a phase needs an even grid, got {grid}")
    values = sorted(class_values(grid))
    middle = len(values) // 2
    if len(values) % 2:
        return (values[middle],)
    return tuple(sorted({values[middle - 1], values[middle]}))


def period_of(durations: tuple[int, ...] | list[int]) -> int:
    """The period a set of service lengths generates: their least common multiple.

    A system whose lengths are ``15, 20, 30, 45, 60`` does not repeat every 15
    or every 60 units -- it repeats every 180, and the candidate cadences of
    that repetition are the divisors of 180 above one.  So the period is the
    single number that turns a schedule somebody already runs into the lattice
    the rest of this package prices, and the ceiling ``nu(T)`` is then a
    property of *their* system rather than of a period picked by example.

    The lengths are a property of the instance, so a non-positive one is a
    request defect and not an infeasibility: ``ValueError``, never a number.
    """
    period = 1
    for d in durations:
        if d < 1:
            raise ValueError(f"a duration must be a positive number of units: {d}")
        period = period * d // math.gcd(period, d)
    return period


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


def grids_by_trial(period: int) -> tuple[int, ...]:
    """Every grid of ``period``, by trying every candidate below it.

    :func:`grids` builds the divisors from the factorization; this one divides.
    Same set, and the only route in the package that would notice a
    factorization that lost a prime or an exponent.
    """
    if period < 2:
        raise ValueError(f"a period is at least two, got {period}")
    return tuple(M for M in range(2, period + 1) if period % M == 0)


def universe_by_median(period: int) -> tuple[int, ...]:
    """``U(period)`` read from the classes, on a grid table built by division.

    The third reading of the same set.  It shares no step with
    :func:`universe` beyond the definition itself: the grids come from trial
    division instead of from the factorization, and the carrier of each grid
    from the median of its class values instead of from accumulated mass.
    """
    if period < 2:
        raise ValueError(f"a period is at least two, got {period}")
    found: set[int] = set()
    for M in grids_by_trial(period):
        if M % 2 == 0:
            found.update(median_carriers(M))
    return tuple(sorted(found))


def median_oracle_census(max_period: int = 600) -> dict:
    """Both sets read a third time, on every even period to ``max_period``.

    The package has shipped two routes to ``U(T)`` for three versions and has
    said, correctly, that they are not independent: both end at
    :func:`carriers`, so a bias in ``mass`` or in ``divisors`` would survive
    the double reading and look like agreement.  This is the census that
    settles it, and it is deliberately built from the cheapest possible
    machinery -- ``gcd``, a sort and a division -- so that there is nothing
    left for the two routes to share.

    Both sets are read, not just the universe, and that is not belt and braces:
    ``U(T)`` is a **union** over the grids, so it absorbs a wrong carrier until
    a period arrives whose other grids do not repair it.  Read the tie branch
    out of :func:`median_carriers` and no universe below ``36`` changes -- the
    grids of ``4``, ``8`` and ``12`` are covered by other grids -- while the
    grid parts at ``4`` itself.  The universe is the claim; the grid is the
    instrument.

    A grid is the same object in every period that holds it, so each is read
    once, and ``grids_measured`` counts distinct grids rather than occurrences.
    """
    if max_period < 2:
        raise ValueError(f"a horizon is at least two, got {max_period}")
    grid_splits: list[int] = []
    carrier_splits: list[tuple[int, tuple[int, ...], tuple[int, ...]]] = []
    universe_splits: list[tuple[int, tuple[int, ...], tuple[int, ...]]] = []
    periods = grids_read = 0
    largest_universe = 0
    for grid in range(2, max_period + 1, 2):
        grids_read += 1
        by_median, by_mass = median_carriers(grid), carriers(grid)
        if by_median != by_mass:
            carrier_splits.append((grid, by_median, by_mass))
    for period in range(2, max_period + 1, 2):
        periods += 1
        if grids_by_trial(period) != grids(period):
            grid_splits.append(period)
        by_median = universe_by_median(period)
        by_mass = universe(period)
        largest_universe = max(largest_universe, len(by_mass))
        if by_median != by_mass:
            universe_splits.append((period, by_median, by_mass))
    if grid_splits:
        verdict = (f"the grids disagree at {grid_splits[0]}: trial division and "
                   f"the factorization do not name the same set, which is the "
                   f"factorization failing and not the horizon")
    elif carrier_splits:
        grid, by_median, by_mass = carrier_splits[0]
        verdict = (f"the carriers disagree at grid {grid}: the median of the "
                   f"class values reads {list(by_median)} and the mass reads "
                   f"{list(by_mass)}, which is a bias in the mass or in the "
                   f"divisor list and not a horizon")
    elif universe_splits:
        period, by_median, by_mass = universe_splits[0]
        verdict = (f"the universes disagree at {period}, though no single grid "
                   f"does, which is a union being read differently and not a "
                   f"carrier being read differently")
    elif periods == 0:
        verdict = "NO_MEDIDO: no even period to read"
    else:
        verdict = (f"the median of the class values names every carrier the "
                   f"mass names, over {grids_read} even grids of {periods} "
                   f"periods to {max_period}, so U(T) is a property of the "
                   f"lattice and not of the route that reads it")
    return {
        "max_period": max_period,
        "periods_measured": periods,
        "grids_measured": grids_read,
        "largest_universe": largest_universe,
        "grid_splits": tuple(grid_splits),
        "carrier_splits": tuple(carrier_splits),
        "universe_splits": tuple(universe_splits),
        "verdict": verdict,
    }


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
