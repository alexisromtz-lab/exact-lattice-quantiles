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

The closed form also has a hypothesis nobody wrote down: the price of a grid is
``grid // 2`` and nothing else.  Charge a fixed cost ``setup`` per grid bought
and the formula stops being the price at ``setup = 1`` -- :func:`setup_fee_frontier`
finds that first fee and the request that breaks it, and the theorem behind it
is one line: folding a block of ``k`` grids into its lcm saves ``k - 1`` fees
and pays the rest of the bill with ``sum(b // 2) - lcm(b) // 2``.

That identity has a consequence worth exporting: the fee at which a request
breaks the formula is decided by an antichain alone, because the formula prices
the *maximals* of a request, and :func:`block_fee` inverts the identity in one
integer division.  So the frontier of a period is a minimum over its antichains
rather than a walk over fees.

The size of the block that attains it is a *ladder with three rungs*, and the
top rung is free: a pair needs ``2`` and attains it at ``(2, 3)``, a triple can
need ``1``, and a quadruple can need ``0``.  The first period that reaches zero
is ``462 = 2 * 3 * 7 * 11``, on ``(42, 66, 154, 231)`` -- the cofactors of its
primes.  Nothing above four buys anything, because a fee is a count and ``0``
is its floor.  :func:`arity_ladder_census` measures that ladder, and
:func:`block_free` explains its last rung: a block is free exactly when its
parts are *dearer* than their lcm by more than its odd parts are numerous,
which for the cofactors of ``L`` is Giuga's inequality ``sum 1/p > 1`` read
strictly past the boundary -- on the boundary itself the block pays exactly
``1``, and :func:`free_margin` turns the distance into the integer a horizon
can be checked against.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt, lcm

from .arithmetic import divisors, prime_factorization, primes_of
from .lattice import cell_widths, grids, mask_table, universe

_INF = 10**9


def maximal_divisors(needed: tuple[int, ...]) -> tuple[int, ...]:
    """The elements of ``needed`` that no other element is a multiple of."""
    return tuple(d for d in needed
                 if not any(e != d and e % d == 0 for e in needed))


def price_one_per_maximal(needed: tuple[int, ...], setup: int = 0) -> int:
    """The closed form: buy one grid per maximal divisor.

    Always a feasible price for ``needed``, never below the minimum -- and it
    buys one grid per maximal, so a fixed cost ``setup`` per grid is paid once
    per maximal.
    """
    return sum(setup + d // 2 for d in maximal_divisors(needed))


def block_advantage(parts: tuple[int, ...], setup: int = 0) -> int:
    """Points a block of grids saves over buying its parts, fee included.

    Positive means the block is cheaper than its parts, so the closed form --
    which never buys a block -- is not the price on that request.  The identity
    is ``setup * (len(parts) - 1) + sum(p // 2) - lcm(parts) // 2``: folding
    ``k`` grids into one lcm refunds ``k - 1`` fees and pays the rest of the
    bill with the parts' weight against the lcm's.
    """
    if not parts:
        raise ValueError("a block needs at least one grid")
    return sum(setup + p // 2 for p in parts) - (setup + lcm(*parts) // 2)


def _maximal_table(U: tuple[int, ...]) -> list[int]:
    """``mask -> bitmask`` of its maximal positions, built one bit at a time.

    Adding ``d`` to a set drops the maximals it divides and keeps ``d`` itself
    only when no surviving maximal is a multiple of it.
    """
    n = len(U)
    table = [0] * mask_table(n, "the maximal-divisor table")
    for mask in range(1, len(table)):
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
def _cheapest(period: int, reduced: bool, setup: int = 0
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

    size = mask_table(len(U), "the cheapest-cover table", period=period)
    f = [_INF] * size
    f[0] = 0
    for M, mask in covers:
        nxt = f[:]
        for s in range(size):
            v = f[s]
            cost = setup + M // 2
            if v < _INF and v + cost < nxt[s | mask]:
                nxt[s | mask] = v + cost
        f = nxt
    # A request is met by any cover of a superset of it.
    for i in range(len(U)):
        bit = 1 << i
        for s in range(size):
            if not s & bit and f[s | bit] < f[s]:
                f[s] = f[s | bit]
    return U, tuple(f), tuple(covers)


def price_exact(period: int, needed: tuple[int, ...], reduced: bool = True,
                setup: int = 0) -> int:
    """Cheapest design covering every element of ``needed``.

    The minimum over designs, by dynamic programming, and the state is the
    **antichain of maximal elements of the request**, not the universe: a grid
    that covers ``b`` covers every divisor of ``b`` for free, so a request costs
    exactly what its maximals cost.  Covering more never costs less, so the
    minimum over covers of a superset is the minimum over covers of the
    maximals, which is why the two agree on every request -- checked against the
    ``2^|U|`` table on all 8 388 607 requests over the divisors of 360.

    The state space is therefore ``2^k`` with ``k`` the number of maximals, and
    that is what makes the price of ``{54, 90, 135}`` in period ``1080``
    answerable: the full table there is ``2^23`` masks, the maximals are three.
    """
    U = universe(period) if reduced else grids(period)
    for d in needed:
        if d not in U:
            raise ValueError(
                f"request {needed} is not coverable in period {period}")
    A = maximal_divisors(tuple(sorted(set(needed))))
    k = len(A)
    if k == 0:
        return 0
    size = mask_table(k, "the frontier price table", period=period)
    best = [_INF] * size
    best[0] = 0
    for M in grids(period):
        mask = 0
        for i, a in enumerate(A):
            if M % a == 0:
                mask |= 1 << i
        if not mask:
            continue
        cost = setup + M // 2
        nxt = list(best)
        for s in range(size):
            v = best[s]
            if v < _INF and v + cost < nxt[s | mask]:
                nxt[s | mask] = v + cost
        best = nxt
    return best[size - 1]


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
    setup: int = 0
    witness: tuple[int, ...] = ()
    witness_gap: int = 0
    witness_exact: int = 0
    witness_form: int = 0
    pair_failures: int = 0

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
            "setup": self.setup,
            "witness": list(self.witness),
            "witness_gap": self.witness_gap,
            "witness_exact": self.witness_exact,
            "witness_form": self.witness_form,
            "pair_failures": self.pair_failures,
        }


def price_form_sweep(low: int, high: int, reduced: bool = True,
                     max_masks: int = 1 << 22,
                     setup: int = 0) -> dict[str, object]:
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
        row = price_form_census(period, reduced=reduced, setup=setup)
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
        "setup": setup,
    }


def price_form_census(period: int, reduced: bool = True,
                      setup: int = 0) -> PriceFormCensus:
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

    ``setup`` is a fixed cost charged once per grid bought.  At ``setup = 0``
    the price of a grid is ``grid // 2`` and nothing else, which is the
    hypothesis the closed form never states: ``witness`` and ``pair_failures``
    publish where that hypothesis is what fails, and :func:`setup_fee_frontier`
    finds the first ``setup`` that breaks it.
    """
    U, best, covers = _cheapest(period, reduced, setup)
    n = len(U)
    table = _maximal_table(U)
    form_of: dict[int, int] = {}
    seen: set[int] = set()
    saving: set[int] = set()
    consumers = discrepancies = 0
    witness: tuple[int, ...] = ()
    witness_gap = witness_exact = witness_form = 0
    pair_failures: set[int] = set()
    for mask in range(1, len(table)):
        exact = best[mask]
        maxmask = table[mask]
        form = form_of.get(maxmask)
        if form is None:
            form = sum(setup + U[i] // 2 for i in range(n) if maxmask >> i & 1)
            form_of[maxmask] = form
        consumers += 1
        if exact != form:
            discrepancies += 1
            saving.add(maxmask)
            if form - exact > witness_gap:
                witness_gap = form - exact
                witness = tuple(U[i] for i in range(n) if maxmask >> i & 1)
                witness_exact, witness_form = exact, form
            if bin(maxmask).count("1") == 2:
                pair_failures.add(maxmask)
        seen.add(maxmask)
    full = len(best) - 1
    one_grid = min((setup + M // 2 for M, mask in covers if mask == full),
                   default=0)
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
        setup=setup,
        witness=witness,
        witness_gap=witness_gap,
        witness_exact=witness_exact,
        witness_form=witness_form,
        pair_failures=len(pair_failures),
    )


def setup_fee_frontier(period: int, reduced: bool = True,
                       max_setup: int = 6) -> dict[str, object]:
    """The least fixed cost per grid at which the closed form stops being the price.

    The closed form buys one grid per maximal divisor and never buys a block, so
    its price for a request with maximals ``A`` is ``sum(setup + a // 2)``.  A
    block ``B`` costs ``setup + lcm(B) // 2``, and the difference is a one-line
    identity (:func:`block_advantage`)::

        advantage(B, setup) = setup * (len(B) - 1) + sum(b // 2) - lcm(B) // 2

    Two consequences decide the frontier.  First, a block of ``k`` grids folds
    ``k - 1`` fees into one, so the fee **favours the big block**: the requests
    that break first are not the small pairs but the triples.  Second, for an
    incomparable pair ``lcm(a, b) >= a + b`` gives ``sum(b // 2) <= lcm // 2``,
    strictly, so ``advantage <= -1`` and **no pair can win below a fee of two**:
    a census that only looked at pairs would report the closed form surviving
    one fee further than it does.  The bound is attained rather than slack --
    ``(2, 3)`` sits at exactly ``0`` once a fee is charged, so a fee of two is
    not merely sufficient for a pair to win, it is the first that can.

    So the frontier is measured over every request of the instance rather than
    argued from pairs: ``rows`` is the census at each fee from zero to
    ``max_setup``, ``frontier`` is the first row that is not adequate, and
    ``pair_frontier`` is the first row where a pair is what fails (``None``
    inside the horizon, which is the theorem above holding on the population).
    """
    if max_setup < 0:
        raise ValueError("max_setup counts units of a fee, so it cannot be negative")
    rows = tuple(price_form_census(period, reduced=reduced, setup=s)
                 for s in range(max_setup + 1))
    broken = next((r for r in rows if not r.adequate), None)
    pairs = next((r for r in rows if r.pair_failures), None)
    if broken is None:
        verdict = (
            f"NO_MEDIDO: no fee up to {max_setup} breaks the closed form on "
            f"period {period}; the frontier is beyond the horizon, not absent")
    else:
        verdict = (
            f"the closed form stops being the price at a fee of {broken.setup} "
            f"per grid: {broken.discrepancies} of {broken.consumers} requests, "
            f"witness {broken.witness} costs {broken.witness_exact} exact "
            f"against {broken.witness_form} by the formula")
    return {
        "period": period,
        "universe": "reduced" if reduced else "divisors",
        "universe_size": rows[0].universe_size,
        "max_setup": max_setup,
        "rows": rows,
        "frontier": None if broken is None else broken.setup,
        "discrepancies": 0 if broken is None else broken.discrepancies,
        "consumers": rows[0].consumers,
        "witness": () if broken is None else broken.witness,
        "witness_gap": 0 if broken is None else broken.witness_gap,
        "witness_exact": 0 if broken is None else broken.witness_exact,
        "witness_form": 0 if broken is None else broken.witness_form,
        "pair_frontier": None if pairs is None else pairs.setup,
        "theorem": (
            "advantage(B, setup) = setup * (len(B) - 1) + sum(b // 2) - lcm(B) // 2; "
            "an incomparable pair sits at -1 or below on a zero fee, so no fee "
            "below 2 lets a pair win"),
        "verdict": verdict,
    }


def block_defect(parts: tuple[int, ...]) -> int:
    """How far a block is from perfect: its sum minus its lcm.

    Positive is the normal case -- the lcm of a set divides the product, so a
    block almost always costs less than its parts and the difference is what the
    block refunds. ``0`` is a *perfect* block (lcm equal to the sum) and negative
    is a block whose lcm exceeds its sum, which can only happen with parts that
    share no common multiple smaller than the product.
    """
    if len(parts) < 2:
        raise ValueError("a block of one grid has no defect to measure")
    total = 0
    for a in parts:
        total += a
    return total - lcm(*parts)


def block_fee_closed(parts: tuple[int, ...]) -> int:
    """The fee of ``block_fee`` as one division on the defect, exactly.

    ``block_advantage`` is affine in the fee with slope ``len(parts) - 1``, so
    the least fee that makes it positive is ``floor(B / (n - 1)) + 1`` where
    ``B = lcm // 2 - sum(p // 2)``.  ``B`` is not free: writing ``delta`` for the
    defect and ``k`` for the number of odd parts,

    - when the lcm is even, ``floor(lcm / 2) = lcm / 2`` and every odd part
      contributes ``(p - 1) / 2``, so ``B = (k - delta) / 2``;
    - when the lcm is odd every part is odd, so ``B = (n - 1 - delta) / 2``.

    Twice ``B`` is therefore an integer, and the fee is that integer divided by
    ``2 * (n - 1)``, plus one.  The naive version of this -- reading the fee off
    the defect alone -- is wrong, and ``(6, 10, 15)`` is the witness: its defect
    is one, so ``-delta / (2 * (n - 1))`` is negative and the block would look
    free, when the least fee that breaks the formula with it is one.  The
    correction is the odd-part count ``k``: a block of even grids gets
    ``p // 2 = p / 2`` exactly, an odd part loses half a unit to the floor and
    pays it back as ``k``.
    """
    if len(parts) < 2:
        raise ValueError("a block of one grid saves nothing and is not a block")
    n = len(parts)
    lcm_all = lcm(*parts)
    delta = block_defect(parts)
    if lcm_all % 2 == 0:
        twice_b = sum(1 for p in parts if p % 2) - delta
    else:
        twice_b = (n - 1) - delta
    if twice_b < 0:
        return 0
    return twice_b // (2 * (n - 1)) + 1


def block_fee(parts: tuple[int, ...]) -> int:
    """The least fixed cost per grid at which the block beats buying its parts.

    ``block_advantage`` is affine in the fee with slope ``len(parts) - 1``, so
    the least fee that turns it positive is one integer division away, and a
    block that is already cheaper with no fee charged reports ``0``.
    """
    if len(parts) < 2:
        raise ValueError("a block of one grid saves nothing and is not a block")
    advantage = block_advantage(parts)
    return 0 if advantage > 0 else (-advantage) // (len(parts) - 1) + 1


def block_free(parts: tuple[int, ...]) -> bool:
    """Whether the block already beats its parts with no fee charged at all.

    Read off the defect this is backwards on purpose: a free block is one whose
    parts are *dearer* than their lcm.  With ``delta`` the defect and ``k`` the
    number of odd parts, ``2 * B = k - delta`` when the lcm is even, so the free
    case ``B < 0`` needs ``delta > k``.  What such a block saves is the fee, not
    the grid -- it pays more for its parts than one grid of their lcm costs and
    still comes out ahead, because it buys several grids and pays the fee once.
    A pair can never do it: ``lcm(a, b) > a + b`` for incomparable ``a``, ``b``,
    so a pair's defect is negative and its fee is never below ``2``.
    """
    if len(parts) < 2:
        raise ValueError("a block of one grid is not a block")
    return block_advantage(parts) > 0


def free_margin(parts: tuple[int, ...]) -> int:
    """How far a block is from free, as an integer that is zero on the line.

    :func:`block_free` is this number's sign and nothing else, and a sign is
    all a verdict needs.  The distance is what a *horizon* needs, and a fee
    cannot supply it: ``block_fee`` divides the margin by twice one less part,
    so on a triple the single fee ``1`` is charged for margin ``0`` and for
    margin ``-2`` alike.  A census that reads fees therefore cannot tell a
    period sitting exactly on the free line -- ready to cross at the next
    multiple -- from one two below it, and that is how a published horizon
    stopped at ``1000`` with the first carrier crossing at ``1080`` in plain
    sight.

    Doubling the advantage is the exact reading: ``2 * advantage = delta - k``
    when the lcm is even and ``delta - (n - 1)`` when it is odd, so a block on
    an odd lcm is free as soon as ``delta >= k`` while one on an even lcm needs
    ``delta > k``.  No block sits on that difference in either universe up to
    ``1200``, so the shorter reading ``delta - k`` measures the same thing
    everywhere this package looks -- but the definition takes the exact form,
    because a horizon argument is an inequality and not a census.
    """
    if len(parts) < 2:
        raise ValueError("a block of one grid has no margin to be free by")
    return 2 * block_advantage(parts)


def cofactor_block(primes: tuple[int, ...]) -> tuple[int, ...]:
    """The block of the cofactors of ``L = prod(primes)``: one part per prime.

    Every part is ``L // p``, so the block is the set of maximals of the divisor
    lattice of the squarefree ``L`` and its lcm is ``L`` itself.  Its defect is
    ``L * (sum 1/p) - L`` exactly, which makes the whole Giuga family a
    statement about a reciprocal sum: ``delta = 1`` is Giuga's equation
    ``sum 1/p - 1/L = 1`` (``30``, ``858``, ``1722`` are the known ones),
    ``delta = -1`` is the primary pseudoperfect ``sum 1/p + 1/L = 1``, and the
    block is free once ``delta`` climbs past the number of odd parts -- which
    for a set containing ``2`` exactly once is ``delta >= 2``.
    """
    if len(set(primes)) != len(primes) or len(primes) < 2:
        raise ValueError("a cofactor block needs two distinct primes at least")
    L = 1
    for p in primes:
        if p < 2:
            raise ValueError("a cofactor block is indexed by primes, not by units")
        L *= p
    return tuple(sorted(L // p for p in primes))


def _antichains(U: tuple[int, ...], arity: int):
    """Every antichain of ``U`` with two to ``arity`` parts, in universe order.

    Only an antichain can make the closed form wrong: it prices the maximals of
    a request, so a part that divides another is never priced separately.
    """
    for k in range(2, min(arity, len(U)) + 1):
        for A in combinations(U, k):
            if all(max(x, y) % min(x, y) for x, y in combinations(A, 2)):
                yield A


def fee_frontier_by_blocks(period: int, reduced: bool = True,
                           arity: int = 3) -> dict:
    """The least fee per grid at which an antichain of at most ``arity`` parts wins.

    The closed form is wrong on a request exactly when some antichain of that
    request is cheaper bought as one grid than bought apart, and the fee at
    which that happens is :func:`block_fee` of the antichain.  So the frontier of
    a period is the minimum of that over the antichains of its universe -- and
    what a census can then ask is how many parts the minimum needs.  The answer
    is a ladder of three rungs, measured by :func:`arity_ladder_census`: a pair
    reaches ``2``, a triple reaches ``1``, and a quadruple reaches ``0``, which
    is as low as a count goes.
    """
    if arity < 2:
        raise ValueError("arity counts parts of a block, so it starts at two")
    U = universe(period) if reduced else grids(period)
    best, witness = None, ()
    pair, pair_witness = None, ()
    for A in _antichains(U, arity):
        fee = block_fee(A)
        if len(A) == 2 and (pair is None or fee < pair):
            pair, pair_witness = fee, A
        if best is None or fee < best:
            best, witness = fee, A
            if best == 0:
                break
    return {
        "period": period,
        "universe": "reduced" if reduced else "divisors",
        "universe_size": len(U),
        "arity": arity,
        "frontier": best,
        "witness": witness,
        "chain": best is None,
        "pair_frontier": pair,
        "pair_witness": pair_witness,
    }


def _exhaustive_fee(U: tuple[int, ...]) -> tuple[int | None, tuple[int, ...]]:
    """The frontier over *every* antichain of ``U``, whatever its size."""
    best, witness = None, ()
    for A in _antichains(U, len(U)):
        fee = block_fee(A)
        if best is None or fee < best:
            best, witness = fee, A
            if best == 0:
                break
    return best, witness


def fee_frontier_census(max_period: int = 2000, reduced: bool = True,
                        arity: int = 3, verify_universe: int = 14) -> dict:
    """The fee frontier over a whole range of periods, with its law and its check.

    Two questions, kept apart.  The **law** is the minimum over antichains of at
    most ``arity`` parts and it is computed for every period, in time
    polynomial in the size of the universe.  The **check** is exhaustive over
    every antichain of any size, and it is only affordable while the universe is
    small, so it runs on the periods with ``|U| <= verify_universe`` and the
    census publishes how many those are rather than pretending to all.

    Three states come out of it, and they are not the same answer.  A period
    whose universe is a chain has **no** antichain at all, so no fee ever makes
    the closed form wrong there -- and measured to 2000 those are exactly the
    ten powers of two, which is the half that can be argued (the divisors of
    ``2^k`` are totally ordered, so no antichain exists) with the converse left
    to the census.  A frontier of ``1``
    needs a triple, because an incomparable pair sits at ``-1`` or below on a
    zero fee (``lcm`` of two incomparable grids is at least twice the larger, so
    the block's price never drops under the sum of theirs).  Everything above
    ``1`` is a pair failing, and the first such period is ``6``, where ``(2, 3)``
    costs ``2`` apart and ``3`` whole: a fee of two is exactly what it takes.

    What this census cannot say is what a block *larger* than ``arity`` would
    do, and that silence is where the claim used to be too strong: the frontier
    law stopped at three parts on the strength of a horizon that reaches
    ``verify_universe`` grids, while the first block that reaches a fee of zero
    needs four parts and a universe of ``15`` -- one grid past the horizon.  The
    ladder of arities is :func:`arity_ladder_census`'s question, and it is asked
    with its own horizon rather than inherited from this one.
    """
    if max_period < 2:
        raise ValueError("the census needs at least one even period")
    if verify_universe < 2:
        raise ValueError("a universe of one grid has no antichain to verify")
    period_rows, mismatches, chains, at_zero, at_one = [], [], [], [], []
    verified, largest, witness_sizes = 0, None, {}
    pair_over, pair_misses = [], []
    for period in range(2, max_period + 1, 2):
        row = fee_frontier_by_blocks(period, reduced=reduced, arity=arity)
        period_rows.append(row)
        if row["frontier"] is None:
            chains.append(period)
        else:
            if row["frontier"] == 0:
                at_zero.append(period)
            elif row["frontier"] == 1:
                at_one.append((period, row["witness"]))
            if largest is None or row["frontier"] > largest["frontier"]:
                largest = {"period": period, "frontier": row["frontier"],
                           "witness": row["witness"]}
            if (row["pair_frontier"] is not None
                    and row["pair_frontier"] > row["frontier"]):
                pair_misses.append(period)
            if (row["pair_frontier"] is not None
                    and row["pair_frontier"] > row["frontier"]
                    and row["frontier"] >= 2):
                pair_over.append({"period": period,
                                  "frontier": row["frontier"],
                                  "pair_frontier": row["pair_frontier"],
                                  "witness": row["witness"]})
        if row["universe_size"] <= verify_universe:
            verified += 1
            truth, truth_witness = _exhaustive_fee(
                universe(period) if reduced else grids(period))
            if truth != row["frontier"]:
                mismatches.append({"period": period, "law": row["frontier"],
                                   "exhaustive": truth,
                                   "exhaustive_witness": truth_witness})
            elif truth is not None:
                size = len(truth_witness)
                witness_sizes[size] = witness_sizes.get(size, 0) + 1
    universe_name = "U(T)" if reduced else "D(T)"
    if arity < 3:
        pair_law = (
            f"NOT_APPLICABLE: at {arity} parts the law is already the pairs, so "
            f"no triple is enumerated and nothing here measures what one buys")
    else:
        pair_law = (
            f"{len(pair_misses)} periods put a triple above every pair and all "
            f"of them sit at a fee of at most one, so the third part buys "
            f"exactly one thing: the fee of one"
            if not pair_over else
            f"REFUTED at {len(pair_over)} periods: a triple reaches a fee below "
            f"every pair of its universe while the fee is already two or more")
    if verified == 0:
        verdict = (
            f"NO_MEDIDO: no period up to {max_period} has a universe of at most "
            f"{verify_universe} grids, so the law was computed but never "
            f"checked against the exhaustive frontier")
    elif mismatches:
        first = mismatches[0]
        verdict = (
            f"REFUTED: the law prices period {first['period']} at fee "
            f"{first['law']} while the exhaustive frontier is "
            f"{first['exhaustive']} -- a block of more than {arity} parts wins "
            f"earlier than every block of at most {arity}")
    else:
        shapes = sorted({w for _, w in at_one})
        if largest is None:
            widest = ("no period in this horizon has an antichain, so no fee "
                      "breaks any of them")
        else:
            widest = (f"the widest frontier is {largest['frontier']} at period "
                      f"{largest['period']} on {tuple(largest['witness'])}")
        verdict = (
            f"the law holds on all {verified} period"
            f"{'s' if verified != 1 else ''} of {universe_name} with at most "
            f"{verify_universe} grids, witnesses of size "
            f"{sorted(witness_sizes)}; {len(chains)} periods are chains where no "
            f"fee breaks the formula, {len(at_zero)} break with no fee at all"
            f"{' ' + str(at_zero) if at_zero else ''}, and {len(at_one)} break "
            f"first at a fee of 1 on {len(shapes)} shape"
            f"{'s' if len(shapes) != 1 else ''} "
            f"{', '.join(str(s) for s in shapes) if shapes else '(none)'}; "
            f"{widest}")
        if arity >= 3:
            verdict += (
                "; blocks above that arity are NOT_MEASURED here, and a fourth "
                "part reaches a fee of 0 -- `arity-ladder` is the census that "
                "asks how many parts the frontier needs")
    return {
        "periods": len(period_rows),
        "max_period": max_period,
        "universe": universe_name,
        "arity": arity,
        "verify_universe": verify_universe,
        "verified": verified,
        "mismatches": mismatches,
        "chains": tuple(chains),
        "at_zero": tuple(at_zero),
        "at_one": tuple(at_one),
        "largest": largest,
        "verified_witness_sizes": witness_sizes,
        "pair_misses": tuple(pair_misses),
        "pair_over": tuple(pair_over),
        "pair_law": pair_law,
        "rows": tuple(period_rows),
        "theorem": (
            "the frontier of a period is the minimal fee over its antichains; a "
            "pair needs 2 and attains it at (2, 3); a triple can need only 1; "
            "and a universe that is a chain never breaks. The minimum is "
            "attained by a pair unless it is 0 or 1, so the third part buys "
            "only the fee of one -- and the fourth buys zero, which this census "
            "does not measure unless arity reaches it"),
        "verdict": verdict,
    }


def arity_ladder_census(max_period: int = 1000, reduced: bool = False,
                        arities: tuple[int, ...] = (2, 3, 4, 5),
                        exhaustive_until: int = 15) -> dict:
    """How many parts the frontier needs: the ladder, rung by rung.

    A fee is a minimum over blocks, so a larger block can only lower it, and the
    question is where the lowering stops.  It stops at zero, which is also where
    it stops *earning*: no rung above one that already reaches ``0`` can improve
    on it, because a fee is a count.  So the ladder has two earning steps at
    most -- ``2 -> 3`` and ``3 -> 4`` -- and this census measures which periods
    each one moves, on what witness, and where the free rung first opens.

    The answer depends on the universe, and that is the point.  Over the grids
    the free rung opens at ``60`` on ``(12, 20, 30)``, which **three** parts
    already reach; the quadruple only earns where the triple stalls at one, and
    the first such period is ``462`` on the cofactors of ``2 * 3 * 7 * 11``.
    Over :func:`universe` -- the carriers, which is the package's default -- no
    period up to ``1000`` reaches zero at all, and the first that does is
    ``1080`` on ``(54, 90, 135)``: nine times the Giuga triple, free by scaling
    rather than by arity.  That crossing is what :func:`free_margin_census`
    measures, and the ``80`` periods between the two numbers are the price of
    reading a horizon off a fee instead of off the distance.
    Both readings are published: a free block is a reciprocal sum over the
    complementary divisors passing ``1`` by more than the share of its odd
    parts, and whether that happens at three parts or at four is a property of
    the universe, not of the fee.

    This is the census the fee frontier was missing.  Its sister
    :func:`fee_frontier_census` checks the law against every antichain, but only
    on the periods its horizon reaches, and that horizon was one grid short of
    the first four-part win: ``462`` has fifteen grids above one and the check
    stopped at fourteen, so its green was a green in the dangerous direction.
    A verification horizon is only as strong as its distance to the first
    witness, so this census names both -- ``first_free`` carries the first
    witness and ``exhaustive`` carries what ``exhaustive_until`` actually
    reached.
    """
    if max_period < 2:
        raise ValueError("the ladder needs at least one period to climb")
    rungs = tuple(sorted(set(arities)))
    if len(rungs) < 2 or rungs[0] < 2:
        raise ValueError("a ladder needs two rungs and a block starts at two")
    steps = tuple(zip(rungs, rungs[1:]))
    earnings: dict[tuple[int, int], list[int]] = {step: [] for step in steps}
    first: dict[tuple[int, int], dict] = {}
    at_floor: list[int] = []
    first_free: dict | None = None
    verified = 0
    mismatches: list[dict] = []
    exhaustive_rows: list[dict] = []
    for period in range(2, max_period + 1, 2):
        fees: dict[int, int | None] = {}
        witnesses: dict[int, tuple[int, ...] | None] = {}
        sizes: dict[int, int] = {}
        for n in rungs:
            row = fee_frontier_by_blocks(period, reduced=reduced, arity=n)
            fees[n], witnesses[n] = row["frontier"], row["witness"]
            sizes[n] = row["universe_size"]
        if all(fee is None for fee in fees.values()):
            continue
        for lo, hi in steps:
            if fees[lo] is None or fees[hi] is None or fees[hi] >= fees[lo]:
                continue
            earnings[(lo, hi)].append(period)
            first.setdefault((lo, hi), {
                "period": period,
                "universe_size": sizes[hi],
                "from_fee": fees[lo],
                "to_fee": fees[hi],
                "witness": witnesses[hi],
            })
        if fees[rungs[-1]] == 0:
            at_floor.append(period)
        if first_free is None:
            for n in rungs:
                if fees[n] == 0:
                    first_free = {
                        "period": period,
                        "parts": n,
                        "universe_size": sizes[n],
                        "witness": witnesses[n],
                        "fourth_part_needed": fees[rungs[-1]] == 0 and n == rungs[-1],
                    }
                    break
        if sizes[rungs[-1]] <= exhaustive_until:
            verified += 1
            truth, _ = _exhaustive_fee(
                universe(period) if reduced else grids(period))
            ladder = min(fee for fee in fees.values() if fee is not None)
            exhaustive_rows.append({"period": period,
                                   "universe_size": sizes[rungs[-1]],
                                   "frontier": ladder, "exhaustive": truth})
            if truth != ladder:
                mismatches.append({"period": period, "ladder": ladder,
                                   "exhaustive": truth})
    heads = " ".join(f"{lo}->{hi}: {len(earnings[(lo, hi)])}" for lo, hi in steps)
    if mismatches:
        verdict = (
            f"REFUTED: the ladder is not the frontier at {len(mismatches)} "
            f"periods, so a block above {rungs[-1]} parts wins where it cannot")
    elif verified == 0:
        verdict = (
            f"NO_MEDIDO: no period up to {max_period} has at most "
            f"{exhaustive_until} grids, so the ladder was never checked "
            f"against every antichain")
    elif first_free is None:
        verdict = (
            f"NO_PERIOD REACHES ZERO up to {max_period} over the "
            f"{'reduced' if reduced else 'divisor'} universe: the ladder earns "
            f"{heads} period(s) and no rung is free, so the floor of the fee is "
            f"one here rather than zero")
    else:
        earned_step = next((step for step in reversed(steps)
                            if earnings[step]), None)
        if earned_step is not None:
            earned = first[earned_step]
            where = (f"the last step that earns is {earned_step[0]}->"
                     f"{earned_step[1]}, at {len(earnings[earned_step])} "
                     f"period(s), first {earned['period']} on "
                     f"{earned['witness']}")
            ceiling = f"no block above {earned_step[1]} parts earns over one at {earned_step[1]}"
        else:
            where = "no step earns at all"
            ceiling = f"no block above {rungs[0]} parts earns over one at {rungs[0]}"
        verdict = (
            f"the ladder earns {heads} period(s); a free block first appears at "
            f"{first_free['period']} on {first_free['witness']}, which "
            f"{first_free['parts']} part(s) already reach, so {where}; "
            f"{ceiling} on any period up to {max_period}, and the ladder "
            f"matches the exhaustive frontier on all {verified} periods whose "
            f"universe is at most {exhaustive_until} grids")
    return {
        "max_period": max_period,
        "universe": "reduced" if reduced else "divisors",
        "arities": rungs,
        "earnings": {f"{lo}->{hi}": tuple(periods)
                     for (lo, hi), periods in earnings.items()},
        "first": {f"{lo}->{hi}": row for (lo, hi), row in first.items()},
        "first_free": first_free,
        "at_floor": tuple(at_floor),
        "exhaustive": {"until": exhaustive_until, "verified": verified,
                       "rows": tuple(exhaustive_rows),
                       "mismatches": tuple(mismatches)},
        "theorem": (
            "a fee is a minimum over blocks, so a bigger block can only lower "
            "it, and it stops at zero, which is also where it stops earning; a "
            "block is free exactly when the reciprocal sum over its "
            "complementary divisors passes 1 by more than the fraction of its "
            "parts that are odd, no pair can do that (1/2 + 1/3 < 1), and the "
            "free rung therefore needs three parts at least and nothing above "
            "four"),
        "verdict": verdict,
    }


def _best_margin(items: tuple[int, ...]) -> tuple[int | None, tuple[int, ...]]:
    """Largest :func:`free_margin` over every antichain of ``items``, by search.

    Exhaustive by construction rather than by horizon: the recursion walks the
    antichains of the universe itself, so there is no ``exhaustive_until`` to
    calibrate and no verification gap for a witness to hide in.  A block of one
    part is skipped -- a single grid is never a block -- and the bound
    ``2 * (halves + rest) - lcm`` can only overestimate a descendant, because
    the lcm never shrinks along a branch and the parity term it drops only ever
    helps.
    """
    pool = tuple(sorted(set(items)))
    size = len(pool)
    if size < 2:
        return None, ()
    apart = tuple(tuple(i != j and pool[i] % pool[j] and pool[j] % pool[i]
                        for j in range(size))
                  for i in range(size))
    best: int | None = None
    block: tuple[int, ...] = ()

    def walk(start: int, chosen: tuple[int, ...], span: int, halves: int) -> None:
        nonlocal best, block
        if len(chosen) >= 2:
            margin = 2 * halves - span + (span & 1)
            if best is None or margin > best:
                best, block = margin, tuple(pool[j] for j in chosen)
        if start == size:
            return
        rest = sum(pool[j] // 2 for j in range(start, size)
                   if all(apart[j][c] for c in chosen))
        if best is not None and 2 * (halves + rest) - span <= best:
            return
        for j in range(start, size):
            if all(apart[j][c] for c in chosen):
                walk(j + 1, chosen + (j,), lcm(span, pool[j]),
                     halves + pool[j] // 2)

    walk(0, (), 1, 0)
    return best, block


def free_margin_census(max_period: int = 1200, reduced: bool = True) -> dict:
    """The free line read in units of its own distance, not in fees.

    :func:`arity_ladder_census` asks where the floor stops being one; this asks
    the same question with the integer that decides it, and the two answers are
    not congruent because a fee is a floored quotient.  Two rows come out, and
    the first is the one a fee cannot see: ``on_line`` collects the periods
    whose best block has margin exactly ``0``, which are the periods *touching*
    the free line -- free at the next multiple, paying ``1`` today.

    Measured over the carriers to ``1200``: nineteen periods touch the line,
    every one of them a multiple of ``60`` and every one of them on
    ``(6, 10, 15)`` -- the Giuga triple, whose defect equals its single odd
    part exactly.  Exactly one period crosses, ``1080 = 18 * 60``, on
    ``(54, 90, 135) = 9 * (6, 10, 15)`` with margin ``8``: the freeness comes
    from scaling, not from a fourth part.  Over the divisors the line is
    touched at ``30`` and ``858`` -- the two smallest Giuga numbers inside the
    horizon, each on its own cofactor block -- and crossed from ``60`` on.

    The horizon has to reach the witness, and this one is chosen to: ``1080``
    sits eighty periods past the ``1000`` a fee-read census stops at, which is
    the same defect as a verification horizon one grid short of its witness.
    """
    if max_period < 2:
        raise ValueError("the census needs at least one period to measure")
    on_line: list[dict] = []
    below = 0
    first_free: dict | None = None
    best: dict | None = None
    measured = 0
    for period in range(2, max_period + 1, 2):
        pool = universe(period) if reduced else grids(period)
        margin, block = _best_margin(pool)
        if margin is None:
            continue
        measured += 1
        if margin == 0:
            on_line.append({"period": period, "block": block})
        elif margin < 0:
            below += 1
        if margin > 0 and first_free is None:
            first_free = {"period": period, "parts": len(block), "block": block,
                          "margin": margin, "universe_size": len(pool)}
        if best is None or margin > best["margin"]:
            best = {"period": period, "block": block, "margin": margin,
                    "universe_size": len(pool)}
    name = "carriers" if reduced else "divisors"
    if first_free is None:
        verdict = (f"NO_PERIOD_CROSSES: over the {name} the best margin is "
                   f"{best['margin'] if best else 'NO_MEDIDO'} up to "
                   f"{max_period}, so this horizon stops on or below the line")
    else:
        verdict = (f"the {name} cross the free line at {first_free['period']} "
                   f"on {first_free['block']} with margin "
                   f"{first_free['margin']}, after {len(on_line)} periods "
                   f"touching it at margin 0")
    return {
        "universe": name,
        "max_period": max_period,
        "periods_measured": measured,
        "on_line": tuple(on_line),
        "below": below,
        "first_free": first_free,
        "best": best,
        "verdict": verdict,
    }


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


def block_of(n: int) -> tuple[tuple[int, ...], int]:
    """The drop-one block of ``n`` and the grid that covers it whole.

    ``{n // p}`` over the distinct primes ``p`` of ``n``, keeping only the
    entries that are carriers.  The block is an antichain by construction --
    primes are distinct, so no entry divides another -- and the one grid that
    covers all of it is its least common multiple, which is ``n`` only when
    every prime power of ``n`` survives a division.  ``(54, 90, 135)`` is the
    block of ``270``; ``(12, 20, 30)`` the block of ``60``; the prime ``2``
    leaves nothing and the prime square ``4`` leaves ``{2}``, whose least common
    multiple is ``2`` and not ``4``.
    """
    block = tuple(sorted(d for d in (n // p for p in primes_of(n)) if d >= 2))
    if not block:
        return (), 1
    return block, lcm(*block)


def block_saving(n: int) -> int:
    """The points the closed form overcharges on the drop-one block of ``n``.

    ``sum(d // 2 for d in block) - lcm(block) // 2``, exact integer arithmetic.
    It is about ``(n / 2) * (sum(1 / p) - 1)`` over the primes of ``n``, so the
    criterion is a sum of reciprocal primes -- and the floors are not noise:
    ``n = 30`` has reciprocal sum ``31 / 30`` above one and still saves nothing,
    because ``6 + 10 + 15 = 31`` against ``30`` buys the same thirty points.

    A negative result is not an error: ``{2, 3}`` costs ``2`` apart and ``3``
    whole, so a small number is priced *better* by buying its parts.
    """
    block, whole = block_of(n)
    return sum(d // 2 for d in block) - whole // 2


def odd_members(n: int) -> int:
    """How many entries of the drop-one block of ``n`` are odd.

    Every entry is ``n // p``, so the parity is decided by the two-adic valuation
    of ``n`` alone: ``omega`` entries when ``n`` is odd, ``1`` when exactly one
    factor ``2`` is available, none when ``4`` divides ``n``.
    """
    return sum(1 for d in block_of(n)[0] if d % 2)


def band_slack(n: int) -> int:
    """The points the block may overcharge for free: ``odd_members(n) - [n odd]``.

    The block is charged ``d // 2`` per entry and the grid that covers it is
    charged ``lcm // 2``.  Each odd entry rounds down once and the total rounds
    down once, so the two roundings cancel exactly:

        block_saving(n) == (prime_excess(n) - band_slack(n)) // 2

    for every ``n`` with at least two distinct primes -- the unified law behind
    the three regimes of :func:`sign_regime`.  The slack is ``1`` for a period
    with exactly one factor ``2`` (the ``858`` shape), ``0`` when ``4`` divides
    the period (so an excess of ``1`` is unreachable there: ``prime_excess`` is
    even), and ``omega - 1`` for an odd period, where every entry is odd.

    So the band is not an unexplained residue; it is the region
    ``1 <= prime_excess(n) <= band_slack(n)`` of one integer predicate.  See
    :func:`excess_parity` for why the division above needs no floors.

    The value is read from the primes and not from the block, because a prime
    period has an empty block (there is nothing to drop) and would report a
    slack of ``-1``; from ``omega`` the slack is ``omega - 1 = 0`` there, which
    is the reading the law needs.
    """
    if n % 2:
        return len(primes_of(n)) - 1
    return 1 if n % 4 else 0


def prime_excess(n: int) -> int:
    """``n * (sum(1/p) - 1)`` as an exact integer: ``sum(n // p) - n``.

    Multiplying by ``n`` clears every denominator, so the sign of the block's
    excess lives in this integer and not in the reciprocal sum, and the gap to
    one is ``prime_excess(n) / n`` exactly.  ``prime_excess(30) == 1``,
    ``prime_excess(60) == 2``, ``prime_excess(858) == 1``, and an odd period of
    ``omega`` primes is feasible (a reciprocal sum above one) only from
    ``omega == 9`` on, where the first witness already exceeds ``10**8``.
    """
    return sum(n // p for p in primes_of(n)) - n


def excess_parity(n: int) -> int:
    """The parity the excess cannot escape: ``prime_excess(n) === band_slack(n)`` mod 2.

    Every term of ``M = sum(n // p) - n`` is congruent to its own parity modulo
    two, so ``M === #(p : n // p odd) - [n odd] (mod 2)``.  That count is how
    :func:`band_slack` is *defined* -- ``odd_members`` minus the period's own
    parity -- so the congruence is an identity and not a measurement:

        M === odd_entries - [n odd] === slack    (mod 2)

    What it buys, in one line each:

    * ``block_saving(n) == (M - slack) // 2`` is an *exact* division for every
      period with at least two distinct primes, because ``M - slack`` is even by
      this congruence: the two roundings of the charges cancel by construction
      and not by census.
    * ``slack == 0`` (``4`` divides ``n``) empties the band ``1 <= M <= 0``
      because ``M`` is even there -- not a separate fact.
    * ``slack == 1`` (exactly one factor ``2``) forces ``M`` odd, and the band is
      then the single point ``M == 1``: the even half of the band *is* the Giuga
      equation ``sum(1/p) - prod(1/p) == 1``.
    * ``slack == omega - 1`` (odd ``n``) forces ``M === omega - 1 (mod 2)``, so
      ``M == 1`` needs ``omega`` even.  An odd Giuga number would have an even
      number of prime factors, and the window's corner is not Giuga's corner at
      odd ``omega``.

    A prime period has ``omega == 1`` and an empty block, where the law does not
    apply -- ``block_saving`` is zero while ``M - slack`` is not -- so this raises
    instead of handing back a parity the law cannot use.
    """
    if len(primes_of(n)) < 2:
        raise ValueError(f"n must have at least two distinct primes, got {n}")
    return (odd_members(n) - (n % 2)) % 2


def band_window(product: int, excess: int, threshold: int) -> tuple[int, int]:
    """The closed interval of primes that would put ``product * q`` in the band.

    Every proper subset of a band member has a reciprocal sum at most one (see
    :func:`odd_band_census`), so a member is always a prefix with ``excess < 0``
    extended by one prime.  With ``m = -excess`` the member's excess is
    ``product - q * m``, which decreases in ``q`` -- hence the three regions are
    contiguous and the band is one interval:

        ``product - q * m`` in ``[1, threshold]``
            <=> ``q`` in ``[ceil((product - threshold) / m), (product - 1) // m]``

    A negated result is possible and means no prime can close the band here.
    The interval has width below ``threshold + 1`` integers, so at most
    ``threshold`` primes can close it above any one prefix: the band is a
    window of width ``omega - 1``, never a surface to be scanned.
    """
    if product < 1:
        raise ValueError(f"product must be positive, got {product}")
    if excess >= 0:
        raise ValueError(f"excess must be negative here, got {excess}")
    if threshold < 1:
        raise ValueError(f"threshold must be at least 1, got {threshold}")
    slack = -excess
    low = -((-(product - threshold)) // slack)
    return low, (product - 1) // slack


def carrier_block_of(n: int, carried: tuple[int, ...] | set[int]) -> tuple[tuple[int, ...], int]:
    """The drop-one block of ``n``, kept only where the universe carries it.

    ``U(T)`` is not divisor-closed: it can hold ``84``, ``140`` and ``210`` while
    missing ``60``, even though all four are the block of ``420``.  The filter
    is therefore part of the law and not a detail -- dropping an element can
    only lower the least common multiple that covers the rest, so a filtered
    block can save where the whole one cannot.
    """
    block, _ = block_of(n)
    kept = tuple(d for d in block if d in carried)
    if not kept:
        return (), 1
    return kept, lcm(*kept)


def carrier_block_saving(n: int, carried: tuple[int, ...] | set[int]) -> int:
    """The saving of the block of ``n`` after the universe keeps its part.

    Same integer rule as :func:`block_saving` -- charged apart, ``sum(d // 2)``;
    covered whole, ``lcm(d) // 2`` -- applied to the kept elements only.
    """
    kept, whole = carrier_block_of(n, carried)
    return sum(d // 2 for d in kept) - whole // 2


def antichain_saving_ceiling(period: int) -> tuple[int, tuple[int, ...]]:
    """The most any antichain with this exact least common multiple can save.

    This is the other side of the law.  ``block_saving`` prices one candidate;
    here every antichain of divisors of ``period`` whose least common multiple
    is exactly ``period`` is priced, and the best one wins.  Pairs and single
    grids are excluded because they never save: for incomparable ``a`` and ``b``
    with ``a = g * s`` and ``b = g * t``, ``lcm(a, b) = g * s * t >= a + b``
    because ``(s - 1)(t - 1) >= 1``.

    The search prunes a branch as soon as ``charged apart + everything still
    reachable`` cannot beat the best saving already found, which is why the
    ceiling is reachable for one period but not for every period at once.  The
    branches it keeps are exactly the ones that could save, so a return of
    ``(0, ())`` is a statement about every antichain and not about the search.
    """
    candidates = tuple(sorted((d for d in divisors(period) if d >= 2),
                              reverse=True))
    n = len(candidates)
    suffix = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        suffix[i] = suffix[i + 1] + candidates[i] // 2
    best = [0, ()]

    def walk(start: int, chosen: list[int], merged: int, charged: int) -> None:
        if len(chosen) >= 3 and merged == period and charged - period // 2 > best[0]:
            best[0], best[1] = charged - period // 2, tuple(chosen)
        if start >= n or charged + suffix[start] - merged // 2 <= best[0]:
            return
        for j in range(start, n):
            d = candidates[j]
            if any(d % c == 0 or c % d == 0 for c in chosen):
                continue
            chosen.append(d)
            walk(j + 1, chosen, _lcm(merged, d), charged + d // 2)
            chosen.pop()

    walk(0, [], 1, 0)
    return best[0], best[1]


def block_dominance_census(limit: int) -> dict:
    """Check both readings of the lemma the law rests on, period by period.

    The **sign** reading is the one the law needs: if no block saves, no
    antichain saves either, so a cover built out of several grids cannot beat
    the closed form.  The **ranking** reading is stronger: the block of ``n``
    is the antichain of least common multiple ``n`` that saves *most*.  The law
    only needs the first, and the second is false in general -- see
    :func:`strong_lemma_frontier` for the first period family that refutes it,
    which no period up to ``limit`` reaches.
    """
    violations = []
    sign_violations = []
    in_play = 0
    for period in range(2, limit + 1):
        ceiling, witness = antichain_saving_ceiling(period)
        block = max(0, block_saving(period))
        if ceiling or block:
            in_play += 1
        if ceiling > block:
            violations.append({
                "period": period,
                "ceiling": ceiling,
                "witness": list(witness),
                "block": block,
                "block_witness": list(block_of(period)[0]),
            })
        if ceiling > 0 and block <= 0:
            sign_violations.append({
                "period": period,
                "ceiling": ceiling,
                "witness": list(witness),
            })
    return {
        "limit": limit,
        "periods": max(0, limit - 1),
        "in_play": in_play,
        "violations": violations,
        "sign_violations": sign_violations,
        "verdict": "the block carries every saving" if not violations
                   else "the block does not carry every saving",
        "sign_verdict": "no antichain saves where no block does"
                        if not sign_violations
                        else "an antichain saves where no block does",
    }


def pair_antichain_of(period: int) -> tuple[int, ...]:
    """One grid per pair of primes of ``period``, descending.

    For distinct primes ``p, q`` of ``period`` the grid ``period // (p * q)``
    divides ``period``, and two such grids divide each other only when they come
    from the same pair -- so the family is an antichain.  Its least common
    multiple is ``period`` itself as soon as ``period`` has three primes, since
    every prime keeps its full exponent in the pairs that avoid it.

    With three primes the family is the three products of two primes --
    ``{30, 18, 12}`` for ``180`` -- and it saves nothing at all, because a
    product of two primes is a poor thing to merge: ``18 + 12`` against ``30``.
    The family only pays once the primes are many and individually small, and
    even then it takes a primorial to notice.
    """
    primes = primes_of(period)
    return tuple(sorted((period // (a * b) for a, b in combinations(primes, 2)),
                        reverse=True))


def pair_antichain_saving(period: int) -> int:
    """The points the pair antichain saves, by the rule every block is priced by."""
    grids = pair_antichain_of(period)
    if not grids:
        return 0
    return sum(d // 2 for d in grids) - lcm(*grids) // 2


def _primes_up_to(limit: int) -> list[int]:
    """Every prime at most ``limit``, by sieve."""
    if limit < 2:
        return []
    sieve = bytearray([1]) * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, isqrt(limit) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, limit + 1, i)))
    return [i for i, flag in enumerate(sieve) if flag]


def _first_primes(count: int) -> list[int]:
    """The first ``count`` primes, by sieve."""
    if count <= 0:
        return []
    bound = 16
    while True:
        sieve = bytearray([1]) * (bound + 1)
        sieve[0:2] = b"\x00\x00"
        for i in range(2, isqrt(bound) + 1):
            if sieve[i]:
                sieve[i * i::i] = bytearray(len(range(i * i, bound + 1, i)))
        found = [i for i, flag in enumerate(sieve) if flag]
        if len(found) >= count:
            return found[:count]
        bound *= 2


def strong_lemma_frontier(max_primes: int = 400) -> dict:
    """Where the ranking reading of the lemma breaks: the first primorial that does it.

    The block of ``period`` charges one grid per prime; the pair antichain
    charges one grid per *pair* of primes.  Dividing both by ``period // 2``,
    the block saves when the reciprocal sum of the primes passes ``1`` and the
    pairs save when the second elementary symmetric function of the reciprocal
    primes passes ``1``.  The first grows like ``log log`` and the second like
    its square, so the pairs overtake the block -- and they do it late: at
    ``176`` primes, the largest being ``1049``, where the excess is
    ``0.08`` per cent of a saving with ``440`` digits.

    The crossing is located with exact rational arithmetic and then confirmed
    with the same integer rule the package prices blocks by, floors included,
    on both sides of it.  So the refutation is not an artefact of dividing.
    """
    if max_primes < 3:
        raise ValueError("max_primes must be at least 3")
    primes = _first_primes(max_primes)
    total = Fraction(0)
    squares = Fraction(0)
    below: dict | None = None
    crossing: dict | None = None
    for index, prime in enumerate(primes, 1):
        total += Fraction(1, prime)
        squares += Fraction(1, prime * prime)
        pairs = (total * total - squares) / 2
        row = {"primes": index, "largest_prime": prime,
               "reciprocal_sum": float(total), "pair_sum": float(pairs)}
        if pairs > total:
            crossing = row
            break
        below = row
    if crossing is None:
        return {
            "max_primes": max_primes,
            "crossing_primes": None,
            "largest_prime": None,
            "verdict": "the block still carries the largest saving below "
                       f"{max_primes} primes",
            "crossing": None,
            "below": below,
        }
    for row in (below, crossing) if below else (crossing,):
        prime_set = primes[:row["primes"]]
        product = 1
        for p in prime_set:
            product *= p
        row["block_saving"] = (sum((product // p) // 2 for p in prime_set)
                               - product // 2)
        row["pair_saving"] = (sum((product // (a * b)) // 2
                                  for a, b in combinations(prime_set, 2))
                              - product // 2)
        row["excess"] = row["pair_saving"] - row["block_saving"]
        row["digits"] = len(str(product))
    return {
        "max_primes": max_primes,
        "crossing_primes": crossing["primes"],
        "largest_prime": crossing["largest_prime"],
        "verdict": "the drop-one block is not the best-saving antichain",
        "crossing": crossing,
        "below": below,
    }


def reciprocal_prime_sum(n: int) -> Fraction:
    """``sum(1 / p)`` over the distinct primes of ``n``, exactly.

    The whole sign reading of the lemma turns on this one number.  ``30`` gives
    ``31/30`` -- above one, and its block still saves nothing; ``858`` gives
    ``859/858``, which is one *plus its own reciprocal*; ``2310`` gives
    ``1171/770``.
    """
    if n < 2:
        raise ValueError(f"n must be at least 2, got {n}")
    total = Fraction(0)
    for prime in primes_of(n):
        total += Fraction(1, prime)
    return total


def two_adic(n: int) -> int:
    """The exponent of ``2`` in ``n``, which is ``0`` for odd ``n``."""
    if n < 1:
        raise ValueError(f"n must be a positive integer, got {n}")
    return dict(prime_factorization(n)).get(2, 0)


def sign_regime(period: int) -> dict:
    """Which of the three regimes decides the sign reading for this period.

    ``sum(1/p) <= 1`` -- **proven-zero**.  No antichain of divisors of the
    period saves a point, for every size of period: the block is the heaviest
    antichain and even it does not save.  See :func:`sign_theorem_census` for
    the statement this rests on.

    ``sum(1/p) > 1`` and the block saves -- **block-saves**.  The reading holds
    on its own witness, so nothing is left to prove.

    ``sum(1/p) > 1`` and the block saves nothing -- **band**.  The single
    remaining surface, and it is a Diophantine one: writing the block's saving
    as ``(sum(period / p) - period - odd_block) / 2`` with ``odd_block`` in
    ``{0, 1}``, the band is entered only when ``period * (sum(1/p) - 1) <= 1``,
    which for an even period means exactly ``sum(1/p) = 1 + 1/period`` with the
    period twice an odd number.  Those three members below ``2000`` are ``30``,
    ``858`` and ``1722`` -- and their equation is the Giuga equation
    ``sum(1/p) - prod(1/p) = 1``, because for a squarefree period ``1/period``
    is ``prod(1/p)``.

    All three regimes are now one predicate rather than three cases.  With
    ``m = prime_excess(period) = sum(period / p) - period`` exact and
    ``slack = band_slack(period)``, the unified law is
    ``block_saving = (m - slack) // 2``, so

        band  <=>  1 <= m <= slack

    and ``slack`` is what the parity of the period decides: ``1`` when exactly
    one factor ``2`` is present (the ``30``/``858``/``1722`` equation, which is
    the Giuga equation), ``0`` when ``4`` divides the period (so no member can
    exist there, because ``m`` is even), and ``omega - 1`` for an odd period --
    a *window* of width ``omega`` through Giuga's equation, whose corner
    ``m == 1`` is an odd Giuga number.  None is known and none exists below the
    horizon measured in :func:`odd_band_census`.
    """
    if period < 2:
        raise ValueError(f"period must be at least 2, got {period}")
    primes = primes_of(period)
    reciprocal = reciprocal_prime_sum(period)
    block = block_saving(period)
    exact_band = reciprocal == 1 + Fraction(1, period)
    twos = two_adic(period)
    if reciprocal <= 1:
        regime = "proven-zero"
    elif block > 0:
        regime = "block-saves"
    else:
        regime = "band"
    excess = prime_excess(period)
    slack = band_slack(period)
    if period % 2:
        characterisation = "giuga-window"
    elif twos == 1:
        characterisation = "giuga-equation"
    else:
        characterisation = "no-member-possible"
    return {
        "period": period,
        "primes": list(primes),
        "reciprocal_sum": str(reciprocal),
        "excess": str(period * (reciprocal - 1)),
        "excess_integer": excess,
        "slack": slack,
        "band_predicate": f"1 <= {excess} <= {slack}",
        "in_band": 1 <= excess <= slack,
        "block": block,
        "two_adic": twos,
        "band_equation": exact_band,
        "characterisation": characterisation,
        "characterised": True,
        "regime": regime,
    }


def odd_band_census(max_omega: int = 11, max_prime: int = 150,
                    slack: int | None = None) -> dict:
    """Every odd period the band can reach below a declared horizon: none, so far.

    An odd period has a wider block than an even one -- every ``period / p`` is
    odd, so ``band_slack`` is ``omega - 1`` -- and its band is therefore a
    *window* instead of an equation.  Two facts make the window searchable:

    *Minimal feasibility.*  ``M(P U {q}) = q * M(P) + L(P)`` is exact, so once a
    set of odd primes has a reciprocal sum above one, every extension of it has
    ``M >= 3 + 3 ** |P|``.  With ``|P| >= 9`` that is above ``omega - 1`` for
    every ``omega < 19686``; with ``|P| <= 8`` no odd set crosses one at all,
    because the first eight odd primes sum to ``111429982 / 111546435``.  Hence
    **every proper subset of a band member has a reciprocal sum at most one**,
    and a member is a minimal feasible set extended by a single prime.

    *The window.*  Above such a prefix the member's excess is
    ``product - q * m``, decreasing in ``q``, so the band is exactly the
    interval :func:`band_window` returns: of width ``omega``, so at most
    ``omega`` primes wide.

    Together they fix the shape: at most ``max_omega`` odd primes, each at most
    ``max_prime``, walk every prefix that has not yet crossed one and price the
    one-prime extensions of it.  ``prefixes`` reports how much was walked, so an
    empty answer is a measurement and not an empty loop.

    The boundary is named rather than hidden: ``M == 1`` on an odd period *is*
    an odd Giuga number, ``sum(1/p) - prod(1/p) == 1``, open since 1950 -- no
    horizon settles that corner.  :func:`excess_parity` narrows it: the excess of
    an odd period is congruent to ``omega - 1`` modulo two, so the corner sits at
    even ``omega``, and at odd ``omega`` the band is walked on even excesses
    only.  The rest of the band -- ``2 <= M <= omega - 1`` -- is a different
    equation, and nothing below the horizon enters either of them.

    What the walk returns is not only "no members".  ``margin`` carries the
    nearest miss -- the smallest positive excess it saw, with the primes that
    produce it, per ``omega`` -- so an empty answer is a distance with a witness
    and not a silence, and ``parity`` re-checks the congruence on every one of
    those misses.  That miss is the nearest over the whole horizon and not just
    over the prefixes walked: a set with a positive excess has every superset
    further above one (each prime multiplies the excess and adds the product),
    so the least positive excess sits on a set whose proper subsets are all
    non-positive -- and such a set is a one-prime extension of a negative
    prefix, which is exactly what the walk prices.

    ``slack`` replaces the per-member threshold ``omega - 1`` with a constant
    one.  That is how the instrument is shown not to be blind: with the
    threshold lifted the same walk must return members, or the census is
    measuring its own pruning.
    """
    if max_omega < 2:
        raise ValueError(f"max_omega must be at least 2, got {max_omega}")
    if max_prime < 3:
        raise ValueError(f"max_prime must be at least 3, got {max_prime}")
    odd = tuple(p for p in _primes_up_to(max_prime) if p % 2)
    floor_omega, floor_row = _odd_feasibility_floor(odd)
    below_floor = floor_omega is None or max_omega < floor_omega
    threshold_of = (lambda omega: omega - 1) if slack is None else (lambda _: slack)
    suffix = [Fraction(0)] * (len(odd) + 1)
    for i in range(len(odd) - 1, -1, -1):
        suffix[i] = suffix[i + 1] + Fraction(1, odd[i])
    bounds = [(s.numerator, s.denominator) for s in suffix]
    members: list[dict] = []
    nearest: tuple[int, tuple[int, ...], int, int] | None = None
    closest: dict[int, tuple[int, tuple[int, ...], int, int]] = {}
    prefixes = 0

    def note(omega: int, chosen: tuple[int, ...], prime: int, member: int,
             edge: int) -> None:
        """Record how far above the band this prefix's nearest extension lands.

        ``member`` falls as the prime rises, so the *last* positive extension of
        a prefix is its closest approach from above -- one row per prefix, no
        per-prime bookkeeping.  The members say where the band is entered; this
        says how close it is never entered, which turns "no members" from an
        empty answer into a distance with a witness.
        """
        nonlocal nearest
        seen = closest.get(omega)
        improves_global = nearest is None or member < nearest[0]
        improves_omega = seen is None or member < seen[0]
        if not (improves_global or improves_omega):
            return
        primes = (*chosen, prime)
        if improves_global:
            nearest = (member, primes, omega, edge)
        if improves_omega:
            closest[omega] = (member, primes, omega, edge)

    def walk(start: int, chosen: tuple[int, ...], product: int, excess: int) -> None:
        nonlocal prefixes
        prefixes += 1
        numerator, denominator = bounds[start]
        if -excess * denominator >= numerator * product:
            return  # the reachable reciprocals cannot lift the sum above one
        threshold = threshold_of(len(chosen) + 1)
        omega = len(chosen) + 1
        pending: tuple[int, int] | None = None
        for i in range(start, len(odd)):
            prime = odd[i]
            member = prime * excess + product
            if member >= 1:
                if member <= threshold:
                    members.append({
                        "primes": [*chosen, prime],
                        "omega": omega,
                        "product": product * prime,
                        "excess": member,
                        "slack": len(chosen),
                    })
                else:
                    pending = (prime, member)
                continue
            if pending is not None:
                note(omega, chosen, pending[0], pending[1], threshold)
                pending = None
            if omega < max_omega:
                walk(i + 1, (*chosen, prime), product * prime, member)
        if pending is not None:
            note(omega, chosen, pending[0], pending[1], threshold)

    walk(0, (), 1, -1)
    if below_floor:
        verdict = (f"NO_MEDIDO: at most {max_omega} odd primes cannot have a "
                   f"reciprocal sum above one, so this horizon cannot reach "
                   f"the band at all")
    elif members:
        verdict = (f"the band has {len(members)} odd members below the "
                   f"horizon: {[m['primes'] for m in members]}")
    else:
        verdict = (f"no odd period in the band over {prefixes} prefixes: at "
                   f"most {max_omega} odd primes, each at most {max_prime}")
    if nearest is None:
        margin = {"nearest": None, "by_omega": {},
                  "verdict": "the walk never produced a positive excess, so "
                             "there is no distance to report"}
    else:
        def row_of(entry: tuple[int, tuple[int, ...], int, int]) -> dict:
            member, primes, omega, edge = entry
            return {"primes": list(primes), "omega": omega, "excess": member,
                    "threshold": edge, "distance": member - edge}

        margin = {
            "nearest": row_of(nearest),
            "by_omega": {str(omega): row_of(entry)
                         for omega, entry in sorted(closest.items())},
            "verdict": f"the nearest miss clears the band by "
                       f"{nearest[0] - nearest[3]} at omega = {nearest[2]} "
                       f"{nearest[1]}, excess {nearest[0]} against a threshold "
                       f"of {nearest[3]}",
        }
    return {
        "max_omega": max_omega,
        "max_prime": max_prime,
        "threshold": "omega - 1" if slack is None else slack,
        "odd_primes_available": len(odd),
        "prefixes": prefixes,
        "candidates": members,
        "status": "below-floor" if below_floor else "measured",
        "floor": floor_row,
        "margin": margin,
        "parity": {
            "rule": "for an odd period, prime_excess(n) === omega - 1 (mod 2)",
            "corner": "an excess of 1 would be an odd Giuga number, and the "
                      "congruence puts that corner at even omega only: it is "
                      "unreachable at odd omega",
            "consistent": all((entry[0] - (omega - 1)) % 2 == 0
                              for omega, entry in closest.items()),
        },
        "boundary": ("an odd member with excess 1 would be an odd Giuga number, "
                     "sum(1/p) - prod(1/p) = 1, open since 1950 -- and the "
                     "parity of the excess puts that corner at even omega"),
        "verdict": verdict,
    }


def _odd_feasibility_floor(odd: tuple[int, ...]) -> tuple[int | None, dict]:
    """The smallest number of odd primes whose reciprocal sum can pass one."""
    total = Fraction(0)
    product = 1
    for count, prime in enumerate(odd, start=1):
        below = total
        total += Fraction(1, prime)
        product *= prime
        if total > 1:
            return count, {
                "omega": count,
                "primes": list(odd[:count]),
                "product": product,
                "excess": sum(product // q for q in odd[:count]) - product,
                "below": f"sum(1/p) over the first {count - 1} odd primes is "
                         f"{below}, below one",
            }
    return None, {
        "omega": None,
        "primes": list(odd),
        "product": product,
        "excess": None,
        "below": f"the {len(odd)} odd primes at most {odd[-1] if odd else 0} "
                 f"sum to {total}, below one",
    }


def worst_antichain_weight(period: int) -> tuple[int, tuple[int, ...]]:
    """The largest ``sum(a)`` over antichains of proper divisors with lcm ``period``.

    Proper means ``2 <= a <= period // 2``: the period itself is an antichain of
    one element, it weighs ``period``, and it saves nothing, so admitting it
    would make the quantity say nothing.

    This is the quantity the lemma bounds.  A saving is
    ``(sum(a) - period - odd_slack) / 2``, so once ``sum(a)`` cannot exceed the
    weight of the drop-one block -- which is ``sum(period // p)`` -- the block
    is the *heaviest* antichain, and a period whose block saves nothing has no
    saving antichain at all.  :func:`sign_theorem_census` reports where the
    maximum meets that bound, which is the sharpness of the inequality and not
    an accident of the search.  :func:`heaviest_antichain` is the same maximum
    read against that bound one period at a time, with the warrant attached.
    """
    if period < 2:
        raise ValueError(f"period must be at least 2, got {period}")
    candidates = tuple(sorted((d for d in divisors(period)
                               if 2 <= d <= period // 2), reverse=True))
    n = len(candidates)
    suffix = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        suffix[i] = suffix[i + 1] + candidates[i]
    best = [0, ()]

    def walk(start: int, chosen: list[int], merged: int, total: int) -> None:
        if merged == period and total > best[0]:
            best[0], best[1] = total, tuple(chosen)
        if start >= n or total + suffix[start] <= best[0]:
            return
        for j in range(start, n):
            d = candidates[j]
            if any(d % c == 0 or c % d == 0 for c in chosen):
                continue
            chosen.append(d)
            walk(j + 1, chosen, _lcm(merged, d), total + d)
            chosen.pop()

    walk(0, [], 1, 0)
    return best[0], best[1]


def heaviest_antichain(period: int) -> dict:
    """The heaviest antichain of ``period``, read against the block's own weight.

    The drop-one block ``{period / p}`` is an antichain whose least common
    multiple is ``period``, and its weight is ``sum(period / p)``.  So that
    number is a *lower* bound on the maximum for every period -- the block is
    always a candidate -- and the lemma of :func:`sign_theorem_census` supplies
    the upper bound exactly where ``sum(1/p) <= 1``.  Two lines, and the
    inequality in the docstring of :func:`worst_antichain_weight` becomes an
    equality:

        max sum(a)  ==  sum(period / p)   for every period with sum(1/p) <= 1,

    attained by the block and by nothing heavier.  That is the extremal reading
    of the lemma, and it is what the closed form needs: with the maximum equal
    to the block's weight, "no antichain saves where no block does" is not a
    search that came out empty, it is an identity.

    The band is the one region where the upper bound is not proved, so there the
    equality is a *measurement* and this function labels it as such.  What it
    measures there is the tightest form of the sign reading.  A saving needs
    ``weight >= period + slack + 2`` while the block weighs ``period + m`` with
    ``m = prime_excess(period)``, so

        miss = slack - m + 2,

    which is ``2`` exactly when ``m == slack`` -- the top of the band, the last
    value that still saves nothing.  ``weight`` and ``period + slack`` share a
    parity, so no antichain can miss by less than two points, and two points of
    weight *is* one unit of saving: the tightest band member is one saving away
    and cannot be closer.  ``miss`` is non-positive exactly when the block
    already saves, i.e. exactly outside the band.

    That makes the band's two instruments one statement: ``m == slack`` is the
    condition for tightness here and the threshold :func:`odd_band_census`
    searches for, so an empty odd band is also the statement that no odd period
    is anywhere near one unit of saving.

    And the equality is not a theorem about all periods: at ``176`` primes the
    pair family weighs more than the block, as :func:`strong_lemma_frontier`
    refutes with exact integers.  So the warrant is per region and the region is
    part of the answer -- a gap of zero measured here is not a claim about a
    period with a hundred and seventy-six prime factors.
    """
    primes = primes_of(period)
    if len(primes) < 2:
        raise ValueError(f"period must have at least two distinct primes, "
                         f"got {period}")
    weight, witness = worst_antichain_weight(period)
    excess = prime_excess(period)
    bound = period + excess
    drop = tuple(sorted((period // p for p in primes), reverse=True))
    regime = sign_regime(period)["regime"]
    slack = band_slack(period)
    miss = slack - excess + 2
    return {
        "period": period,
        "primes": list(primes),
        "regime": regime,
        "warrant": "theorem" if regime == "proven-zero" else "measured",
        "excess": excess,
        "bound": bound,
        "weight": weight,
        "gap": bound - weight,
        "witness": list(witness),
        "drop_one": list(drop),
        "block_attains": tuple(witness) == drop,
        "slack": slack,
        "miss": miss,
        "tightest_miss": miss == 2,
    }


def extremal_census(limit: int = 2000) -> dict:
    """The extremal reading over a population, with the warrant of each region.

    Three regions, three warrants, one number -- ``gap = sum(period / p) -
    max sum(a)``:

        region          gap      warrant
        proven-zero      0       theorem: the block attains the bound and the
                                 lemma bounds every antichain
        block-saves      0       measured: the proof does not reach here
        band             0       measured: and `miss = slack - m + 2` is a
                                 minimum only at `m == slack`, so a band member
                                 is tight iff it sits at the top of the window

    A negative gap is the only outcome that refutes the package: the block would
    not be the heaviest antichain, an antichain could save where the block saves
    nothing, and the sign reading would fall.  That is what this census reports
    as a violation, and it is the code it returns.

    Measured below ``200 000`` (``181 879`` periods, ``173 537`` proved, ``8 338``
    block-saves, ``4`` band members: ``30``, ``858``, ``1722``, ``66198`` -- the
    four Giuga numbers below ``200 000``): the gap is zero on every one of them
    and the heaviest antichain is the drop-one block on every one of them.

    The horizon is declared and not implied.  The equality is *false* in general
    -- see :func:`strong_lemma_frontier`, where the pairs outweigh the block from
    ``176`` primes on -- so this census measures a region and names it; it does
    not certify the whole family of periods, which no finite population could.
    """
    if limit < 6:
        raise ValueError(f"limit must be at least 6, got {limit}")
    regions: dict[str, dict] = {}
    band_rows: list[dict] = []
    violations: list[dict] = []
    tightest: dict | None = None
    periods = 0
    for period in range(2, limit + 1):
        if len(primes_of(period)) < 2:
            continue
        row = heaviest_antichain(period)
        periods += 1
        bucket = regions.setdefault(row["regime"], {
            "rows": 0, "gap_zero": 0, "block_attains": 0,
            "tightest_miss": 0, "warrant": row["warrant"], "max_gap": None,
        })
        bucket["rows"] += 1
        bucket["gap_zero"] += row["gap"] == 0
        bucket["block_attains"] += row["block_attains"]
        bucket["tightest_miss"] += row["tightest_miss"]
        if bucket["max_gap"] is None or row["gap"] > bucket["max_gap"]:
            bucket["max_gap"] = row["gap"]
        if row["gap"] < 0:
            violations.append(row)
        if row["regime"] != "block-saves" and (
                tightest is None or (row["miss"], row["period"]) < (
                    tightest["miss"], tightest["period"])):
            tightest = row
        if row["regime"] == "band":
            band_rows.append(row)
    attained = sum(bucket["block_attains"] for bucket in regions.values())
    return {
        "limit": limit,
        "periods": periods,
        "regions": {name: dict(row) for name, row in sorted(regions.items())},
        "band_rows": band_rows,
        "violations": violations,
        "tightest": tightest,
        "block_attains_everywhere": attained == periods,
        "theorem": "max sum(a) == sum(period / p) for every period with "
                   "sum(1/p) <= 1, attained by the drop-one block",
        "frontier": "sum(period / p) is not an upper bound in general: the pair "
                    "family outweighs the block from 176 primes on, largest "
                    "1049 (strong_lemma_frontier)",
        "verdict": "the block is the heaviest antichain wherever the gap is zero, "
                   "by theorem where the reciprocal sum is at most one and by "
                   "measurement elsewhere in this population"
                   if not violations else
                   "an antichain outweighs the block, so the block is not the "
                   "heaviest: the sign reading falls",
    }


def sign_theorem_census(limit: int = 2000) -> dict:
    """The sign reading of the lemma, proved where it is provable and measured elsewhere.

    The lemma the closed form rests on says: if no block saves, no antichain
    saves.  Half of it is a theorem.

    *Lemma.*  If ``sum(1/p) <= 1`` over the primes of ``period`` and ``C`` is an
    antichain of ``period``-smooth integers above one, then
    ``sum(1/c for c in C) <= sum(1/p)``.  By induction on the largest prime
    ``q``: split ``C`` by the exact power of ``q`` it carries.  Dividing by
    ``q^k`` preserves comparability, so each slice is an antichain of
    ``Q``-smooth integers with ``Q`` the primes below ``q``, and the hypothesis
    gives ``sum(slice) <= q**-k * sum(1/p for p in Q)``; summing the geometric
    series gives ``sum(1/p for p in Q) * q / (q - 1)``, and
    ``sum(1/p for p in Q) <= 1 - 1/q`` is exactly what ``sum(1/p) <= 1`` leaves
    at the largest prime.  Equality holds for the drop-one block, so the lemma
    is sharp.

    *Theorem.*  Writing ``c = period / a`` and ``sum(a) = period * sum(1/c)``,
    the lemma gives ``sum(a) <= sum(period / p)``; a saving is
    ``(sum(a) - period - odd_slack) / 2`` with ``odd_slack`` at most one, so
    ``sum(a) <= period`` would already leave it non-positive, and the lemma
    forces the same conclusion for every size of period.  Integer arithmetic
    then does the rest: a quantity strictly below ``1/2`` and integral is at
    most zero.

    The census verifies the two halves against each other -- the maximum weight
    of :func:`worst_antichain_weight` against the block's weight, and the
    exhaustive ceiling of :func:`antichain_saving_ceiling` against the
    theorem's zero -- and reports the band on its own, because the band is the
    one place where the proof does not reach and a measurement is all there is.

    *The law is a congruence, not a census.*  :func:`excess_parity` shows that
    ``prime_excess(n)`` and ``band_slack(n)`` share their parity for every period
    with two distinct primes, so ``(excess - slack) // 2`` divides exactly and
    the three regimes really are one predicate.  The census re-derives the
    congruence on every period it visits, which makes it a check of a proof
    rather than the proof itself -- the same shape as the lemma: proved once,
    re-measured everywhere.
    """
    if limit < 2:
        raise ValueError(f"limit must be at least 2, got {limit}")
    proven = band = block_saves = tight = 0
    law_periods = 0
    band_rows: list[dict] = []
    lemma_violations: list[dict] = []
    sign_violations: list[dict] = []
    law_violations: list[dict] = []
    mismatches: list[dict] = []
    for period in range(2, limit + 1):
        row = sign_regime(period)
        regime = row["regime"]
        if len(row["primes"]) >= 2:
            law_periods += 1
            if (row["excess_integer"] - row["slack"]) % 2:
                law_violations.append({"period": period,
                                       "excess": row["excess_integer"],
                                       "slack": row["slack"]})
        if regime == "block-saves":
            block_saves += 1
            continue
        if regime == "band":
            band += 1
            if row["characterised"] is False:
                mismatches.append({"period": period,
                                   "reciprocal_sum": row["reciprocal_sum"],
                                   "two_adic": row["two_adic"],
                                   "block": row["block"]})
            ceiling, witness = antichain_saving_ceiling(period)
            row = {**row, "ceiling": ceiling, "witness": list(witness)}
            band_rows.append(row)
            if ceiling > 0:
                sign_violations.append({"period": period, "ceiling": ceiling,
                                        "witness": list(witness),
                                        "regime": "band"})
            continue
        proven += 1
        weight, heavy = worst_antichain_weight(period)
        bound = sum(period // p for p in primes_of(period))
        if weight > bound:
            lemma_violations.append({"period": period, "weight": weight,
                                     "bound": bound, "witness": list(heavy)})
        elif weight == bound:
            tight += 1
        ceiling, witness = antichain_saving_ceiling(period)
        if ceiling > 0:
            sign_violations.append({"period": period, "ceiling": ceiling,
                                    "witness": list(witness),
                                    "regime": "proven-zero"})
    return {
        "limit": limit,
        "periods": max(0, limit - 1),
        "proven_zero": proven,
        "tight": tight,
        "band": band,
        "block_saves": block_saves,
        "band_rows": band_rows,
        "lemma_violations": lemma_violations,
        "sign_violations": sign_violations,
        "characterisation_mismatches": mismatches,
        "law_periods": law_periods,
        "law_violations": law_violations,
        "law": "prime_excess(n) === band_slack(n) (mod 2), an identity, so "
               "block_saving = (excess - slack) // 2 divides exactly",
        "law_verdict": "the unified law is a congruence and holds on every "
                       "period of this census"
                       if not law_violations
                       else "the congruence failed, so the law needs floors",
        "lemma": "sum(1/p) <= 1 implies sum(1/c) <= sum(1/p) over antichains "
                 "of period-smooth integers above one",
        "verdict": "the sign reading is proved where the reciprocal sum is at "
                   "most one and measured where the block is empty-handed"
                   if not lemma_violations else "the lemma is false",
        "sign_verdict": "no antichain saves where the theorem forbids it"
                        if not sign_violations
                        else "an antichain saves where the theorem forbids it",
    }


def block_catalogue(limit: int) -> tuple[tuple[int, int, tuple[int, ...]], ...]:
    """Every ``n <= limit`` whose block saves and that no divisor of it improves.

    The divisibility-minimal set is the useful one: if ``d`` divides ``n`` and
    ``d`` already saves, then ``n`` is a scaled copy of a smaller block.  This
    is the design table -- a period divisible by one of these has a request the
    closed form misprices.
    """
    saving = {n: block_saving(n) for n in range(2, limit + 1)}
    out = []
    for n in range(2, limit + 1):
        if saving[n] <= 0:
            continue
        if any(saving.get(d, 0) > 0 for d in divisors(n) if d != n):
            continue
        out.append((n, saving[n], block_of(n)[0]))
    return tuple(out)


def _block_break(period: int, reduced: bool) -> tuple[int, int, tuple[int, ...], int]:
    """``(saving, n, kept, whole)`` of the biggest saving of any block of a period.

    ``saving`` is zero when no block of the period saves, which is the reading
    the sign lemma argues about; ``n`` and ``kept`` are the witness when it does.
    """
    carried = set(universe(period)) if reduced else set(divisors(period))
    best = (0, 0, (), 0)
    for n in divisors(period):
        if n < 2:
            continue
        kept, whole = carrier_block_of(n, carried)
        saving = sum(d // 2 for d in kept) - whole // 2
        if saving > best[0]:
            best = (saving, n, kept, whole)
    return best


def frontier_by_blocks(max_period: int, reduced: bool = False) -> dict:
    """The frontier over either universe, by one law with the universe as a filter.

    The law is the same sentence in both universes: a period is mispriced when
    some divisor ``n`` of it has a drop-one block whose kept part is cheaper
    bought as one grid than bought apart.  Only the filter changes -- ``D(T)``
    keeps everything, ``U(T)`` keeps the carriers -- and that single difference
    is the whole eighteen-fold gap between the two frontiers.

    Over ``D(T)`` the answer is a divisibility statement: the first period is
    ``60``, because ``{12, 20, 30}`` costs ``31`` apart and ``30`` as the grid
    ``60``.  Measured against the exhaustive search on every even period to
    2000, the law agrees on all 1000 with the same witness and the same saving;
    it also agrees on all 125 breaks to 3000, where no witness of the maximum is
    anything other than a drop-one block.

    Over ``U(T)`` the first period is ``1080``, where ``{54, 90, 135}`` costs
    ``139`` apart and ``135`` as the grid ``270``.  ``840`` is a multiple of
    ``60`` and still does not break: it carries ``20`` and ``12`` but not
    ``30``, so the block of ``60`` survives there with two elements, and two
    incomparable grids never save.  Measured against the exhaustive search on
    every even period to 8000, the law agrees on all 16 breaks -- and on the
    *value* of every saving, not only on break against no-break.
    """
    periods = range(2, max_period + 1, 2)
    breaks = []
    for period in periods:
        best = _block_break(period, reduced)
        if best[0] > 0:
            breaks.append((period,) + best)
    scanned = max_period // 2
    universe_name = "U(T)" if reduced else "D(T)"
    if not breaks:
        return {"period": None, "saving": 0, "witness": (), "bound": 0,
                "cover": 0, "grid": 0, "breaks": 0, "scanned": scanned,
                "universe": universe_name, "law": "carrier blocks"}
    period, saving, n, kept, whole = breaks[0]
    return {
        "period": period,
        "saving": saving,
        "witness": kept,
        "bound": sum(d // 2 for d in kept),
        "cover": whole // 2,
        "grid": whole,
        "breaks": len(breaks),
        "scanned": scanned,
        "universe": universe_name,
        "law": "carrier blocks",
    }


def exactness_frontier(max_period: int = 8000, reduced: bool = True) -> dict:
    """The first even period whose universe the closed form misprices.

    The frontier is a property of the **universe**, and the two answers differ
    by a factor of eighteen: ``60`` over every divisor and ``1080`` over the
    reduced carriers.  Both are decided by the same law -- the carrier block of
    a divisor -- and the only difference between the two universes is which
    elements survive the filter, so neither answer is searched for.  What
    carries each universe, and how the law was checked against the exhaustive
    search there, is stated in :func:`frontier_by_blocks`.

    The witness is certified by the cover, not by a DP: ``{lcm(witness)}`` is a
    legal design that covers the block, so its price is an upper bound the DP
    would only confirm.  At ``1080`` the DP was ``2^23`` masks and at ``2520``
    ``2^32`` -- the certificate has to survive where the search cannot run, and
    :func:`price_exact` now runs there too by taking the maximals only.
    """
    return frontier_by_blocks(max_period, reduced=reduced)


def reciprocal_filter_census(max_period: int = 400) -> dict:
    """Cross the sign lemma's hypothesis with both universes, period by period.

    The lemma is proved from ``lcm(A)`` alone, so it transfers to any set of
    grids a schedule may use -- and a theorem that transfers is only worth the
    row that exercises it.  Every even period up to ``max_period`` is decided
    two ways here: by ``prime_excess`` (the hypothesis ``sum(1/p) <= 1``) and by
    the biggest saving of any of its blocks over ``U(T)`` and over ``D(T)``.

    The two columns separate two different reasons for having no saving.  Below
    one, *no* universe may break, and neither does: that is the proof reaching
    the census.  Above one the proof says nothing, and the period is decided by
    which divisors survive the filter -- ``60`` breaks over every divisor and
    not over the carriers, and ``30`` breaks nowhere, because the band's saving
    is a point the floors absorb.

    A period whose excess is above one and whose two columns disagree is the
    whole finding: same period, same reciprocal sum, same law, other universe.
    """
    rows = []
    for period in range(2, max_period + 1, 2):
        excess = prime_excess(period)
        saving_u, _, kept_u, _ = _block_break(period, True)
        saving_d, _, kept_d, _ = _block_break(period, False)
        rows.append({
            "period": period,
            "excess": excess,
            "above_one": excess > 0,
            "reduced": saving_u,
            "reduced_witness": kept_u,
            "divisors": saving_d,
            "divisors_witness": kept_d,
        })
    proved = [row for row in rows if not row["above_one"]]
    band = [row for row in rows if row["above_one"]]
    clean = [row for row in proved
             if row["reduced"] == 0 and row["divisors"] == 0]
    filtered = [row for row in band if row["reduced"] == 0]
    split = [row for row in band
             if row["reduced"] == 0 and row["divisors"] > 0]
    widest = max(rows, key=lambda row: row["divisors"], default=None)
    if not rows:
        verdict = "no period was censused"
    elif len(clean) < len(proved):
        verdict = "the proof does not reach the census: a period below one breaks"
    elif split:
        verdict = ("the filter decides, not the reciprocal sum: "
                   f"{len(split)} of {len(band)} periods above one break over "
                   "every divisor and none over the carriers")
    else:
        verdict = "no period here separates the two universes"
    return {
        "max_period": max_period,
        "periods": len(rows),
        "proved": len(proved),
        "proved_without_break": len(clean),
        "band": len(band),
        "band_clean_reduced": len(filtered),
        "band_breaking_divisors": len(split),
        "widest": widest,
        "rows": rows,
        "verdict": verdict,
    }
