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

from .arithmetic import divisors, primes_of
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
    best = [_INF] * (1 << k)
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
        for s in range(1 << k):
            v = best[s]
            if v < _INF and v + cost < nxt[s | mask]:
                nxt[s | mask] = v + cost
        best = nxt
    return best[(1 << k) - 1]


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
        carried = set(universe(period)) if reduced else set(divisors(period))
        best = (0, 0, (), 0)
        for n in divisors(period):
            if n < 2:
                continue
            kept, whole = carrier_block_of(n, carried)
            saving = sum(d // 2 for d in kept) - whole // 2
            if saving > best[0]:
                best = (saving, n, kept, whole)
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
