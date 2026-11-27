"""Exact integer arithmetic on the divisor lattice.

Every quantity in this package is an :class:`int` or a
:class:`fractions.Fraction`.  There is no floating point anywhere: a schedule
either has margin at least ``1/6`` or it does not, and the comparison is
decided exactly.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache


@lru_cache(maxsize=None)
def prime_factorization(n: int) -> tuple[tuple[int, int], ...]:
    """``((prime, exponent), ...)`` of ``n`` by trial division, ascending."""
    if n < 1:
        raise ValueError(f"n must be a positive integer, got {n}")
    factors: list[tuple[int, int]] = []
    x, p = n, 2
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            factors.append((p, e))
        p += 1 if p == 2 else 2
    if x > 1:
        factors.append((x, 1))
    return tuple(factors)


@lru_cache(maxsize=None)
def divisors(n: int) -> tuple[int, ...]:
    """The positive divisors of ``n``, ascending, from its factorization."""
    ds = [1]
    for p, a in prime_factorization(n):
        ds = [d * p**i for d in ds for i in range(a + 1)]
    return tuple(sorted(ds))


@lru_cache(maxsize=None)
def primes_of(n: int) -> tuple[int, ...]:
    """The distinct primes dividing ``n``, ascending."""
    return tuple(p for p, _ in prime_factorization(n))


@lru_cache(maxsize=None)
def radical(n: int) -> int:
    """The product of the distinct primes of ``n``, without exponents."""
    r = 1
    for p in primes_of(n):
        r *= p
    return r


@lru_cache(maxsize=None)
def euler_phi(d: int) -> int:
    """Euler's totient of ``d``, by the same division as the factorization."""
    if d < 1:
        raise ValueError(f"d must be positive, got {d}")
    res, x, p = d, d, 2
    while p * p <= x:
        if x % p == 0:
            while x % p == 0:
                x //= p
            res -= res // p
        p += 1
    if x > 1:
        res -= res // x
    return res


_PHI_PREFIX: list[int] = [0, 1]


def totient_sum(n: int) -> int:
    """``euler_phi(1) + ... + euler_phi(n)``, grown one term at a time.

    The carrier bound asks this once per divisor, so it is kept as a running
    prefix: a new ``n`` costs one totient, not ``n`` of them.
    """
    if n < 1:
        raise ValueError(f"n must be a positive integer, got {n}")
    while len(_PHI_PREFIX) <= n:
        _PHI_PREFIX.append(_PHI_PREFIX[-1] + euler_phi(len(_PHI_PREFIX)))
    return _PHI_PREFIX[n]


@lru_cache(maxsize=None)
def mass(d: int) -> int:
    """The number of section classes that capture divisor ``d`` in grid ``d``.

    ``phi(d)`` for every ``d > 2``, but ``2`` is special.  The classes are the
    ``i`` with ``1 <= i <= d // 2`` and ``gcd(i, d) == 1``; for ``d = 2`` that
    set is ``{1}`` -- one class, not half of one.  Writing ``phi(2) / 2 = 1/2``
    here collapses every ceiling in this package, because ``2`` divides every
    even period.  The function returns a count, never a ratio.
    """
    if d < 2:
        raise ValueError(f"a divisor below 2 carries no classes, got {d}")
    return 2 if d == 2 else euler_phi(d)


@lru_cache(maxsize=None)
def odd_part(n: int) -> int:
    """``n`` with every factor of two taken out.

    The cofactor-two law splits by parity, and the even half of it depends on
    ``n`` only through its odd part: doubling or halving a number moves which
    divisors of ``2n`` fall below ``n``, but never the totient ratio that
    decides the comparison.
    """
    if n < 1:
        raise ValueError(f"n must be a positive integer, got {n}")
    while n % 2 == 0:
        n //= 2
    return n


@lru_cache(maxsize=None)
def upper_half_mass(n: int) -> int:
    """The mass of the divisors of ``n`` at or above half of ``n``.

    ``sum(phi(d) for d | n) == n`` always, so this and
    :func:`lower_half_mass` are the two halves of ``n`` and not two readings of
    it: ``upper_half_mass(n) + lower_half_mass(n) == n`` to the unit.  ``n = 1``
    is the degenerate case that keeps both identities honest -- its only
    divisor is at or above half of itself, so the upper half is all of it.
    """
    if n < 1:
        raise ValueError(f"n must be a positive integer, got {n}")
    return sum(euler_phi(d) for d in divisors(n) if 2 * d >= n)


@lru_cache(maxsize=None)
def lower_half_mass(n: int) -> int:
    """The mass of the divisors of ``n`` strictly below half of ``n``."""
    return n - upper_half_mass(n)


@lru_cache(maxsize=None)
def dropped_primes(period: int, d: int) -> tuple[int, ...]:
    """The primes of ``period`` that ``d`` does not have.

    ``d`` divides ``period``, so ``d`` is ``period`` minus a set of primes;
    this is that set, and it is the whole reason the ceiling is a product.
    """
    return tuple(sorted(set(primes_of(period)) - set(primes_of(d))))


@lru_cache(maxsize=None)
def prime_product(period: int) -> Fraction:
    """``P(T) = prod_{p | T} (1 - 1/p)``, the ceiling when nothing is dropped.

    Note what this is *not*: ``mass(radical(T)) / radical(T)``.  For
    ``radical(T) = 2`` that form evaluates to ``1``, which is the one-class
    convention of :func:`mass` at ``d = 2`` and not a ratio of anything.  Every
    power of two has ``P(T) = 1/2`` and a ceiling of ``1/2``; the copied form
    would say ``1`` and be wrong for all of them.
    """
    v = Fraction(1)
    for p in primes_of(period):
        v *= Fraction(p - 1, p)
    return v
