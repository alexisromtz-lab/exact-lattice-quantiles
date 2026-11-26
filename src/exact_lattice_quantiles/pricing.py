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
from itertools import combinations
from math import gcd, isqrt, lcm

from .arithmetic import divisors, prime_factorization, primes_of
from .lattice import cell_widths, grids, mask_table, universe

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

    size = mask_table(len(U), "the cheapest-cover table", period=period)
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
        cost = M // 2
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
    for mask in range(1, len(table)):
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
    full = len(best) - 1
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
    an accident of the search.
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
