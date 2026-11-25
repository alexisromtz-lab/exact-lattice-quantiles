"""The ceiling of one instance: ``nu(T)``, and why it is that number.

The margin of a design is the worst block it can guarantee across the
universe.  Its maximum over every design is not a constant of the problem: it
belongs to the period, and it is read in closed form.

The chain of identities, each measured against the one before it:

1. The seat of ``d`` in grid ``M`` does not depend on ``M``:

       w(M, d) = mass(d) / M          (``mass(2) = 2``, ``mass(d) = phi(d)``)

2. So the widest block ``d`` can ever get is in ``M = d``, and

       nu(T) = min_{d in U(T)} mass(d) / d

3. And since every ``d`` in ``U(T)`` divides ``T``, ``d`` is ``T`` without a
   set ``F`` of primes, which turns the minimum into a product:

       nu(T) = P(T) * prod_{q in F} q / (q - 1),      P(T) = prod_{p | T} (1 - 1/p)

``1/3`` is not a law.  It is the ceiling of ``T = 180`` and of 29 other even
periods below 400 -- while 167 of the 200 exceed it and 3 fall below.  The
ceiling is a function of the instance, and the exponent of a prime decides it
only through which primes ``U(T)`` can afford to drop.

The subject is the even periods.  ``U(T)`` is built from the even grids, and
every divisor of an odd ``T`` is odd, so ``U(T)`` is empty for every odd ``T``:
there is no minimum to take and no ceiling to state.  That boundary is refused
by name (:class:`OddPeriodUnsupported`) at each entry point rather than left to
``min()`` on an empty iterable, whose message names the arithmetic and not the
period.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .arithmetic import dropped_primes, mass, prime_product
from .lattice import universe_from_factorization


@dataclass(frozen=True)
class Ceiling:
    """Everything the ceiling of one period is made of.

    ``margin`` is ``nu(T)``.  ``critical_divisor`` is a divisor that attains
    it -- the one that ties the ceiling, and the only one worth naming when a
    design is being argued about.  ``dropped`` is the set ``F`` of primes the
    universe cannot afford; when it is empty the ceiling is exactly ``P(T)``
    and ``radical_bound_hit`` says so.
    """

    period: int
    margin: Fraction
    critical_divisor: int
    dropped: tuple[int, ...]
    prime_product: Fraction
    universe: tuple[int, ...]
    radical_bound_hit: bool

    @property
    def percent(self) -> Fraction:
        """The margin as a percentage, exactly (``1/3`` -> ``100/3``)."""
        return self.margin * 100


class OddPeriodUnsupported(ValueError):
    """The period is not in the theory's domain, and the refusal says why.

    A section is a class of the phase ``M // 2``, so a grid is a carrier only
    when it is even: an odd ``M`` has no middle class, and its blocks are the
    ones the phase cannot seat.  Every divisor of an odd ``T`` is odd, so
    ``U(T)`` is empty for every odd ``T`` -- there is no universe to take a
    minimum over, and therefore no ceiling and no schedule.  The subject of
    this package is the even periods, and this is the name of that boundary.

    A ``ValueError`` because that is what every entry point raised before --
    except that it raised CPython's ``min() iterable argument is empty``, which
    names the arithmetic and not the subject.
    """

    def __init__(self, period: int) -> None:
        self.period = period
        super().__init__(
            f"period {period} is odd: U(T) is empty because every grid of an odd "
            f"period is odd and only an even grid carries a divisor, so there is "
            f"no ceiling and no schedule -- the subject is the even periods"
        )


class InfeasibleAboveInstanceCeiling(Exception):
    """The margin asked for does not exist at any budget.

    ``needed`` exceeded ``ceiling``.  A caller that reads this and searches
    anyway is searching for a schedule that cannot exist.
    """

    rc = 1

    def __init__(self, period: int, needed: Fraction, ceiling: Fraction) -> None:
        self.period = period
        self.needed = needed
        self.ceiling = ceiling
        super().__init__(
            f"margin {needed} is above the instance ceiling {ceiling} of "
            f"period {period}: no schedule exists at any cost"
        )


def instance_ceiling(period: int) -> Ceiling:
    """``nu(T)`` with the divisor that attains it, in exact arithmetic.

    No grid is built and no seat is scanned: both ``U(T)`` and ``mass(d) / d``
    come from the factorization of ``period``.
    """
    if period < 2:
        raise ValueError(f"a period needs at least two phases, got {period}")
    U = universe_from_factorization(period)
    if not U:
        raise OddPeriodUnsupported(period)
    best = min(Fraction(mass(d), d) for d in U)
    critical = min(d for d in U if Fraction(mass(d), d) == best)
    dropped = dropped_primes(period, critical)
    P = prime_product(period)
    return Ceiling(
        period=period,
        margin=best,
        critical_divisor=critical,
        dropped=dropped,
        prime_product=P,
        universe=U,
        radical_bound_hit=not dropped,
    )


def ceiling_by_product(period: int) -> Fraction:
    """``nu(T)`` written as a product over the primes ``U(T)`` drops.

    ``P(T)`` times a factor ``q / (q - 1)`` for each dropped prime ``q``.  This
    is the same number as :func:`instance_ceiling`; the package keeps both
    routes so that a disagreement is a measurement and not a silence.

    ``T = 2`` is the one exception and it is not a rounding: the whole identity
    rests on ``mass(d) / d`` being a ratio, and at ``d = 2`` it is not --
    ``mass(2) = 2`` while ``phi(2) = 1``.  There ``nu = 1`` and ``P(T) = 1/2``.
    """
    if period == 2:
        return Fraction(1)
    c = instance_ceiling(period)
    out = c.prime_product
    for q in c.dropped:
        out *= Fraction(q, q - 1)
    return out


def ceiling_by_scan(period: int) -> Fraction:
    """``nu(T)`` by the definition, scanning every grid and every divisor.

    ``max`` over grids of the width of ``d``'s block, ``min`` over the
    universe.  Slower than :func:`instance_ceiling` on purpose: it is the
    measurement that route has to reproduce.

    Every grid counts, not only the even ones.  ``U(T)`` is built from the even
    grids, but an odd grid still seats its own divisors, and ``M = d`` is where
    ``d`` is widest.  Reading only even grids halves the ceiling of every odd
    divisor -- ``4/15`` instead of ``1/3`` for ``T = 180``, measured.
    """
    from .lattice import cell_width, grids, universe

    U = universe(period)
    if not U:
        raise OddPeriodUnsupported(period)
    gr = tuple(grids(period))
    return min(max(cell_width(M, d) for M in gr if d in _divisors_of(M)) for d in U)


def _divisors_of(n: int) -> tuple[int, ...]:
    from .arithmetic import divisors

    return divisors(n)


def ceiling_census(low: int, high: int, step: int = 2) -> dict[str, object]:
    """How many periods sit above, on and below ``1/3``, and what binds them.

    The classification is by comparison with ``Fraction(1, 3)`` and nothing
    else; it exists because ``1/3`` circulated as if it were a law of the
    problem, and a census is the cheapest way to keep it from doing so again.
    """
    third = Fraction(1, 3)
    above = same = below = 0
    lowest: tuple[Fraction, int] | None = None
    for T in range(low, high + 1, step):
        nu = instance_ceiling(T).margin
        if nu > third:
            above += 1
        elif nu == third:
            same += 1
        else:
            below += 1
        if lowest is None or nu < lowest[0]:
            lowest = (nu, T)
    return {
        "periods": (high - low) // step + 1,
        "above": above,
        "equal": same,
        "below": below,
        "lowest": lowest,
    }
