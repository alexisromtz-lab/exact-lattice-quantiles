# exact-lattice-quantiles

[![tests](https://github.com/alexisromtz-lab/exact-lattice-quantiles/actions/workflows/tests.yml/badge.svg)](https://github.com/alexisromtz-lab/exact-lattice-quantiles/actions/workflows/tests.yml)

Exact, provable phase margin for multi-rate lattice sampling — in rational
arithmetic, with zero dependencies and zero floats.

A periodic system with period `T` is sampled by **grids**: the divisors of `T`
above one. Grid `M` cuts the phase into `M / 2` section classes, and every
divisor `d` of `T` that divides `M` is seated in a block of them. The
**universe** `U(T)` is the set of *carriers* of the even grids — for each grid,
the divisor whose block straddles the median class. A **design** is a set of
grids; its **margin** is how far a shared section may drift before two
divisors are confused, measured as a fraction of the phase.

Every question about such a system has an exact answer, and this package
returns it as a `fractions.Fraction`:

| question | what it is | function |
| --- | --- | --- |
| What is the best margin this period can ever have? | a closed-form constant of `T`, not a search | `instance_ceiling` |
| What is the cheapest design reaching margin `v`? | exact hypergraph coverage | `minimum_cost_design` |
| Which divisors does that search actually have to tell apart? | the kernel of the coverage matrix | `coverage_kernel` |
| Which margin costs `c` points? | the whole trade-off curve | `ladder` |
| How many shared sections hold margin `v`? | set cover over the divisor universe | `minimum_probes` |
| What is the best margin `k` sections can hold? | maximin over probe positions | `best_probes` |
| Is the margin I asked for even possible? | a certificate, not an error message | `InfeasibleAboveInstanceCeiling` |
| How much common drift does the margin survive? | a curve, not a constant | `drift_tolerance` |
| How much common drift still covers every divisor? | the exact coverage radius | `drift_radius` |
| What does keeping margin `v` under drift `delta` cost? | one purchase, `v + delta` | `guaranteed_margin` |

## Install

```
pip install git+https://github.com/alexisromtz-lab/exact-lattice-quantiles
```

Python 3.11+. No dependencies, at runtime or for tests.

## Quick start

```python
from fractions import Fraction
from exact_lattice_quantiles import (instance_ceiling, period_of,
                                     solve_lattice_schedule)

instance_ceiling(180).margin          # Fraction(1, 3)  — tied by divisor 6
period_of((15, 20, 30, 45, 60))       # 180 — the period the lengths generate
instance_ceiling(360).margin          # Fraction(4, 15) — tied by divisor 90

schedule = solve_lattice_schedule(180, max_cost=64, k_probes=4, min_margin="1/3")
schedule.design                       # (4, 6, 10, 12, 15, 18, 20, 45)
schedule.cost                         # 64
schedule.witness_margin               # Fraction(1, 3)

solve_lattice_schedule(180, max_cost=64, k_probes=4, min_margin="2/5")
# InfeasibleAboveInstanceCeiling: rc 1 — no schedule exists at any cost
```

Or from the shell:

```
$ exact-lattice-quantiles period 15 20 30 45 60   # which period a clinic runs over
$ exact-lattice-quantiles ceiling 180
$ exact-lattice-quantiles solve 180 --max-cost 64 --probes 4 --min-margin 1/3
$ exact-lattice-quantiles ladder 180
$ exact-lattice-quantiles designs 180
$ exact-lattice-quantiles census 2 400
```

## Five claims, each reproducible by a command

### 1. The ceiling belongs to the instance

`1/3` is not a law of the problem. Of the 200 even periods from 2 to 400,
**167 exceed it**, **30 hit it exactly** (180 among them) and **3 fall below**:

```
252 -> 2/7      360 -> 4/15      396 -> 10/33
```

`360 -> 4/15` is the sharpest counterexample to a fixed threshold. 90 and 360
have the same prime support, `{2, 3, 5}`, so no prime-level argument separates
them — what separates them is the radical `rad = 30`:

```
U( 90) = (2, 3, 5, 9, 15, 45)                    no multiple of 30 inside
U(360) = (2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 18,
          20, 24, 36, 40, 45, 90)               90 = 3 * 30 inside
```

Once a multiple of `rad(T)` is a carrier, the ceiling collapses to
`psi(rad)/rad` — here `8/30 = 4/15`, tied by divisor 90 — and no design of any
size can recover the `8/15` that `T = 90` still has:

```
T =  60 -> 1/3     T = 120 -> 1/3     T = 180 -> 1/3     T = 240 -> 1/3
T =  90 -> 8/15    T = 360 -> 4/15
```

Above 400 the same law holds — `720 -> 4/15`, and `840`, `1260`, `2520 ->
8/35` — and it is read in closed form:

    nu(T) = min_{d in U(T)} psi(d) / d,   psi(2) = 2,  psi(d) = phi(d) for d > 2
          = P(T) * prod_{q in F} q / (q - 1),   P(T) = prod_{p | T} (1 - 1/p)

where `F` is the set of primes the universe cannot afford to drop. When `F` is
empty the ceiling is exactly `P(T)`, and `Ceiling.radical_bound_hit` says so.
Two independent routes (`ceiling_by_product` and `ceiling_by_scan`) are tested
against each other on every even period to 400; `T = 2` is the single
degenerate case, where `psi(2) = 2` has no mass to divide.

**The subject is the even periods, and the boundary is named.** A section is a
class of the phase `M // 2`, so only an even grid carries a divisor; an odd `M`
has no middle class. Every divisor of an odd `T` is odd, so `U(T)` is empty for
every odd `T` — measured, `bool(universe(T)) == (T % 2 == 0)` on all 400
periods from 1 to 400. There is no universe to take a minimum over, hence no
ceiling and no schedule, and the package says so instead of letting `min()`
report on an empty iterable:

```
>>> instance_ceiling(9)
OddPeriodUnsupported: period 9 is odd: U(T) is empty because every grid of an
odd period is odd and only an even grid carries a divisor, so there is no
ceiling and no schedule -- the subject is the even periods
```

### 2. The margin/cost curve is exact, not a sample

For `T = 180` — 12 divisors in the universe, 17 candidate grids — the entire
curve, with no threshold guessed anywhere:

| cost (points) | design margin | cheapest design |
| --- | --- | --- |
| 47 | 1/6 | `{12, 18, 20, 45}` |
| 52 | 8/45 | |
| 59 | 1/5 | |
| 64 | 1/3 | `{4, 6, 10, 12, 15, 18, 20, 45}` |

"Margin at least `v`" is a coverage condition — the grids whose cells are at
least `v` wide must cover the whole universe — so its minimum cost is a DP over
`2^|universe|` masks, checked against brute force. `designs 180` counts all
`2^17 = 131072` designs: `96768` cover the universe, and `47` points is the
cheapest covering design.

A lower ceiling is not a cheaper instance. `T = 360` reaches its lower ceiling
`4/15` for **141** points, and its ladder is strictly worse everywhere:
`95 -> 1/12`, `100 -> 4/45`, `107 -> 1/10`, `112 -> 1/6`, `129 -> 1/5`,
`139 -> 1/4`, `141 -> 4/15`.

### 3. Sections are a second budget, and it does not trade linearly

`max_cost` (grid points) and `k_probes` (shared sections) are **not
interchangeable**, and the package keeps them separate: a design can be
affordable and still not be servable.

On `T = 180`:

| design | cost | margin | sections needed | best margin with `k = 4` |
| --- | --- | --- | --- | --- |
| `{12, 18, 20, 45}` | 47 | 1/6 | 8 | 1/15 |
| `{4, 6, 10, 12, 15, 18, 20, 45}` | 64 | 1/3 | 4 | 1/3 |

The cheap design is the expensive one to serve: `1/6` needs eight sections
while the maximum margin needs four, because a design whose blocks are only
just wide enough forces each section to serve a single divisor. At `k = 4` on
the cheap design no set of sections can hold more than `1/15`, and no `k` at
all reaches past `1/6`:

```python
solve_lattice_schedule(180, max_cost=47, k_probes=4, min_margin="1/6")
# verdict = "PROBES_INSUFFICIENT", probes_needed = 8, witness_margin = 1/15
```

### 4. The price is a minimum, not a formula

The tempting formula is "buy one grid per maximal divisor", so the price would
be `sum(d // 2 for d in maximal(needed))`. It is always a *legal* design and it
is the exact optimum over `U(T)` -- but it is not the price over the grids.

A block of grids can be bought whole for less than the sum of its parts:
`{12, 20, 30}` costs `31` bought separately and `30` as the single grid `60`,
their least common multiple. Two incomparable divisors never save
(`lcm(a, b) >= a + b`), so the formula is exact on every request whose every
frontier has no triple, and wrong everywhere else.

```
$ exact-lattice-quantiles price 180
180 over U(T) (12 elements, 4095 requests)
  one grid per maximal    closed form, 117 frontiers
  exact minimum           cheapest cover 47 points, one grid alone 90
  requests overcharged    0 (0) on 0 frontiers
180 over every divisor (17 elements, 131071 requests)
  one grid per maximal    closed form, 173 frontiers
  exact minimum           cheapest cover 90 points, one grid alone 90
  requests overcharged    9984 (9984/131071) on 10 frontiers
```

The boundary is measured, not assumed. Over the reduced universe the formula is
the exact price for **all 200 even periods up to 400**, and it is still never a
single grid that does the work: in 192 of the 200 periods the cheapest cover
beats every single grid.

```
$ exact-lattice-quantiles prices
even periods 2..400 over reduced: 200 censused
  closed form is the price   200
  closed form overcharges    0
  cheaper than one grid      192
  costliest cover            160 points at period 384
  skipped (unaffordable)     0
```

200/200 over `2..400` is a bounded census, not a law. The shortcut breaks the
moment a block of three grids enters `U(T)`: `{6m, 10m, 15m}` has lcm `30m`,
so it costs `15m + floor(m/2)` bought apart and `15m` bought whole, and the
deficit is `floor(m/2)` -- zero at `m = 1`, which is why the frontier is not
period 2. The first mispriced period is `T = 1080`, with the witness
`{54, 90, 135}` costing `139` apart and `135` as the single grid `270`:

```
$ exact-lattice-quantiles frontier
even periods 2..8000 over U(T): 16 of 4000 mispriced, first at T = 1080
  witness                    (54, 90, 135)
  closed form                139 points
  one grid 270               135 points
  saving                     4 points
  decided by                 carrier blocks
```

Below that boundary the formula *is* the price over `U(T)`, and the theorem is
sharp: `16` breaks, `16` triples, no break without a triple and no triple
without a break. Divisibility by `30` is necessary and not sufficient -- `840`
holds `12` and `20` but not `30`, so it is exact.

#### The law behind the boundary, and why the universe decides it

The frontier is a property of the **universe**, and the two answers are eighteen
times apart. One law decides both, and the only thing that changes between them
is a filter.

For `N` let `B(N) = {N // p}` over the distinct primes `p` of `N`: an antichain
by construction, whose least common multiple is `N` when every prime power of
`N` survives the division. Keep from it only the elements the universe actually
carries, `K(N, V) = B(N) ∩ V`, and price what survives by the same integer rule:

    saving(N, V) = sum(b // 2 for b in K(N, V)) - lcm(K(N, V)) // 2

Then the closed form is wrong on the universe `V(T)` **if and only if some `N`
dividing `T` has `saving(N, V) > 0`**. Over `D(T)` the filter is the identity
and the criterion is a plain divisibility statement; over `U(T)` it is the same
statement after the carriers have had their say.

That quantity is about `(N / 2) * (sum(1 / p) - 1)` over the primes of `N`, so
the criterion is a sum of reciprocal primes -- and the floors are not noise:
`N = 30` has reciprocal sum `31 / 30`, above one, and still saves nothing,
because `6 + 10 + 15 = 31` buys the same thirty points as the single grid `30`.

Measured against the exhaustive block search: **0 disagreements** on every even
period to 2000 in *both* universes, and on every even period to 20000 in `U(T)`
-- and not only on break against no-break, the law names the same saving, point
for point (`52` breaks to 20000, all of them equal in value). It runs in
`O(tau(T)^2)` where the search needs the antichains.

```
$ exact-lattice-quantiles frontier --universe all --max-period 2000
even periods 2..2000 over D(T): 82 of 1000 mispriced, first at T = 60
  witness                    (12, 20, 30)
  closed form                31 points
  one grid 60                30 points
  saving                     1 points
  decided by                 carrier blocks

$ exact-lattice-quantiles blocks --limit 250
divisibility-minimal blocks up to 250: 4
      60  saves   1  (12, 20, 30)
      90  saves   1  (18, 30, 45)
     150  saves   2  (30, 50, 75)
     210  saves  18  (30, 42, 70, 105)
```

| universe | frontier | breaks to 8000 | witness |
|---|---|---|---|
| `D(T)`, every divisor | **60** | 336 | `{12, 20, 30}` → `31` / `30` |
| `U(T)`, reduced carriers | **1080** | 16 | `{54, 90, 135}` → `139` / `135` |

Reducing to carriers buys **18x** of delay and changes the witness, not the law.
`U(2520)` keeps `84`, `140` and `210` and loses `60`, so the block of `420`
survives there with three elements: `434` apart against `420` as one grid, a
saving of `7`. `U(840)` keeps `12` and `20` and loses `30`: two elements, and a
pair of incomparable grids never saves -- which is why `840` is a multiple of
`60` and does not break. A frontier quoted without its universe is a number
about the wrong set.

#### The proof transfers to any universe, and the census says where it lands

The sign lemma is proved from `lcm(A)` alone, so no set of grids escapes it: if
`sum(1/p) <= 1` over the primes of `T`, no block saves -- over `U(T)`, over
`D(T)`, over any filter a schedule may apply. A theorem that transfers is worth
exactly the row that exercises it, so `universes` crosses the hypothesis with
both columns on every even period:

```
$ exact-lattice-quantiles universes --max-period 400
reciprocal sum against universe, even periods to 400
  periods                                200
  sum 1/p at most one                    187  no block saves in either universe
  above one                               13  the proof says nothing here
    clean over U(T)                       13  the filter folds the band
    breaking over D(T)                    12  same period, same sum, other universe
  widest gap                           T=390  excess 43  U(T) 0  D(T) 21
  verdict: the filter decides, not the reciprocal sum: 12 of 13 periods above one break over every divisor and none over the carriers
```

`187` is not a claim, it is a count that the census recomputes: it is the number
of even periods to `400` whose reciprocal sum is at most one, and every one of
them is clean in *both* universes. The thirteen above one are where the proof is
silent and the universe answers instead -- `60` breaks over every divisor and
not over the carriers, `180` saves `3` there against `0`, `390` saves `21`
against `0` -- except `30`, which breaks in neither, because its saving is
exactly the point the floors absorb. That row is the band seen from the other
side: same period, same reciprocal sum `31/30`, and no authority but the filter.

**What the law is not.** It is a statement about the *maximum*, not about every
saving set. In `D(1050)` the antichain `{42, 150, 350, 525}` saves `8` points
and `42` is no `1050 / p`: its reciprocals are `1/2 + 1/3 + 1/7 + 1/25 > 1`, so
the drop-one shape is not the only way to save. The maximum still comes from the
block of `1050`, worth `92`.

#### The one lemma that would make the law a theorem

Half of the law needs nothing. If some divisor has a block that saves, that
block **is** a cover, so the closed form that prices its parts is overcharging,
and the period is mispriced. The other half -- *no* saving block means the
closed form is optimal -- is a statement about covers that are not blocks.

Take any cover of a request by several grids and fold it: the grids of a group
are covered by the single grid named by their least common multiple `L`, so the
group costs `L // 2`. Two grids in a group can never be cheaper together than
apart -- for incomparable `a = g * s` and `b = g * t`, `lcm(a, b) / 2 = g * s *
t / 2 >= (a + b) / 2` because `(s - 1)(t - 1) >= 1`. For three or more, the chain
closes if and only if

> **no antichain with least common multiple `L` saves a point where the block
> of `L` saves none.**

With it, a group never beats paying its parts, every cover costs at least the
closed-form price, and the form is optimal exactly when no block saves.

#### The lemma has two readings, and the stronger one is false

Read as a **ranking** -- the block saves *most* among the antichains sharing its
least common multiple -- the lemma is false, and the refutation is short enough
to quote. Subtract the single grid `L // 2` from both savings and divide by it:
the block pays `e1 - 1`, and the family holding one grid per *pair* of primes,
`{L // (pq)}`, pays `e2 - 1`, where `e1` is the reciprocal sum of the prime set
of `L` and `e2` the reciprocal sum over its pairs. Both scale with `L`, so which
wins depends on the prime set alone. `e1` grows like `log log` and `e2` like its
square, so `e2` has to overtake -- but not soon, because a product of two primes
is a poor thing to merge: at `180` the pair family is `{30, 18, 12}` and it
*loses* `60` points. The overtake happens at 176 primes:

```
$ exact-lattice-quantiles strong-lemma
ranking lemma: no antichain saves more than the block of its own least common multiple
  crossing at 176 primes, the largest being 1049
  at 175 primes (largest 1039, 437-digit period)  the block wins by 1.704e+432
  at 176 primes (largest 1049, 440-digit period)  the pairs win by 1.209e+436 (+0.082 %)
  verdict  the drop-one block is not the best-saving antichain
```

The crossing is located with exact rationals and then confirmed with the same
integer rule the package prices blocks by, floors included, on both sides of it,
so the defeat is not an artefact of dividing. It is also out of reach of the
census below: no period under `2000` puts more than four primes in play.

What the law actually uses is the **sign** reading -- no antichain saves where
no block does -- and `dominance` reports it on a line of its own, so the two
readings can never be confused again:

```
$ exact-lattice-quantiles dominance --max-period 2000
antichains with a least common multiple up to 2000: 82 periods where a saving is in play
  every saving is carried by the drop-one block
  no antichain saves where no block does (the reading the law needs)
```

`82` is the same `82` the search finds over `D(T)`, and the ceiling it computes
is the block's saving period by period, witness included --
[tests/test_block_dominance.py](tests/test_block_dominance.py) fails the day an
antichain saves where no block does, and
[tests/test_strong_lemma_frontier.py](tests/test_strong_lemma_frontier.py)
pins the day the ranking dies. The sign reading is no longer a search either:
`sign-theorem` decides it from the prime factorisation alone.

#### The sign reading is an inequality, and it is the one-line kind

Write `e1(P) = sum(1/p)` over the primes of `L`. Subtracting the single grid
`L // 2` from both sides, an antichain `A` saving `S(A)` satisfies

    S(A) = (sum_A a - #odd(A)) / 2 - (L - [L odd]) / 2,

so the whole question is how heavy an antichain can be. Induction on the largest
prime `q` of `L` answers it. Split `A` by `v_q`: dividing a member by `q^k`
preserves divisibility, so each slice is again an antichain, now `q`-smooth
except for the exact power `q^k` charged to it, and the induction hypothesis
applies to the slice without it. The geometric series closes as

    sum_A 1/a  <=  e1(P),   every antichain of proper divisors with lcm = L,

and the block `{L // p}` meets it with equality -- the lemma is sharp, not
slack. Therefore `sum_A a <= L * e1`, hence

> **`e1(P) <= 1` implies `S(A) <= 1/2` for every antichain, and savings are
> integers, so `S(A) <= 0`: no antichain saves a point where no block does.**

The proof never looks at the universe, only at `lcm(A)`, so it transfers to any
set of grids a schedule may use. What it does not cover is the band `e1 > 1`
where the block still saves nothing. There the block's own saving is
`(L * e1 - L - (#odd_B - [L odd])) / 2`, and `#odd_B` -- the odd members of
`{L // p}` -- is what makes a band possible at all. With `v2(L) >= 2` every
member is even, so `#odd_B = 0` and the block saves `L * (e1 - 1) / 2 >= 1`,
which confines the band to `v2(L) = 1`: there the single odd member `L // 2`
buys back exactly one point, and the band is `L * (e1 - 1) = 1`. For squarefree
`L` that is `sum(1/p) - prod(1/p) = 1` -- Giuga's equation -- and the census
below finds its only three members under `2000`:

```
$ exact-lattice-quantiles sign-theorem
sign lemma: if no block saves, no antichain saves
  periods 2..2000                   1999  censused
  sum 1/p at most one               1914  proved: no antichain saves any point
    the block is the heaviest       1581  of those, the lemma is sharp
  the block saves already             82  nothing left to prove
  band, sum 1/p above one              3  measured only, the proof does not reach
        30  sum 1/p = 31/30      period*(sum-1) = 1   block   0  ceiling 0  Giuga equation
       858  sum 1/p = 859/858    period*(sum-1) = 1   block   0  ceiling 0  Giuga equation
      1722  sum 1/p = 1723/1722  period*(sum-1) = 1   block   0  ceiling 0  Giuga equation
  verdict  no antichain saves where the theorem forbids it
```

#### The block is not only allowed, it is the maximum

The lemma above is stated as an inequality: every antichain weighs at most
`e1(P)`. But the drop-one block `{L // p}` *is* one of those antichains, its
least common multiple is `L`, and its weight is exactly `sum(L // p)`. So

    max_A sum_A a  ==  sum(L // p)   whenever  e1(P) <= 1,

where the left side is a maximum over every antichain of proper divisors of `L`.
The lower bound is the block; the upper bound is the lemma; there is nothing
left between them. That matters because it turns the sign reading from a search
that came out empty into an identity: a period whose block saves nothing has no
saving antichain *because nothing else can be heavier than the block*, not
because a walk stopped early.

The equality has three regions and this package refuses to report one number
without the other two. Where `e1(P) <= 1` it is a theorem; where the block
already saves it is measured and harmless; in the band it is measured, and the
measurement is tight, with an identity behind it. A saving needs
`sum_A a >= L + slack + 2` while the block weighs `L + m`, so

    miss  =  slack - m + 2,

which is `2` exactly when `m == slack` -- the top of the band, the last value
that still saves nothing. `sum_A a` and `L + slack` share a parity, so no
antichain can miss by fewer than two points, and two points of weight *is* one
unit of saving: the tightest band member is one saving away and cannot be
closer. `miss` is non-positive exactly when the block already saves, so the
field doubles as the region test. That also makes the band's two instruments one
statement, because `m == slack` is the condition for tightness here and the
threshold `odd-band` searches for: an empty odd band is the same sentence as *no
odd period is anywhere near one unit of saving*.

```
$ exact-lattice-quantiles extremal-gap
extremal reading: the heaviest antichain is the drop-one block
  periods 2..2000                   1666  with at least two distinct primes
  proven-zero                       1581  gap 0 on 1581, warrant theorem
  block-saves                         82  gap 0 on 82, warrant measured
  band                                 3  gap 0 on 3, warrant measured
  the block is the maximum          1666  of those periods, on every one
  band members                    30 858 1722  gap 0, one unit of saving away
  closest to a saving                 30  weight 31 on a bound of 31, 2 points short
  theorem  max sum(a) == sum(period / p) for every period with sum(1/p) <= 1, attained by the drop-one block
  frontier sum(period / p) is not an upper bound in general: the pair family outweighs the block from 176 primes on, largest 1049 (strong_lemma_frontier)
  verdict  the block is the heaviest antichain wherever the gap is zero, by theorem where the reciprocal sum is at most one and by measurement elsewhere in this population
```

A horizon that cannot reach the band exits `4`: the only region whose warrant is
a measurement would be missing from the table, and an empty measured column is
not a green one.

```
$ exact-lattice-quantiles extremal-gap --limit 20 ; echo "rc=$?"
NO_MEDIDO: no band member below 20, so the one region whose warrant is a measurement was not measured
extremal reading: the heaviest antichain is the drop-one block
  periods 2..20                        7  with at least two distinct primes
  proven-zero                          7  gap 0 on 7, warrant theorem
  the block is the maximum             7  of those periods, on every one
  closest to a saving                  6  weight 5 on a bound of 5, 4 points short
  ...
rc=4
```

The horizon is declared, not implied. The equality is **false** in general: the
pair family `{L // (p q)}` is also an antichain with least common multiple `L`,
and at `176` primes -- the largest `1049`, a `440`-digit period -- it outweighs
the block, for the same reason the ranking reading dies there. That failure is
*outside* the region the sign reading needs: it happens where the block already
saves `1.4e+439` points to begin with, so the closed form never consults the
bound there. Below `200 000` the census visits `181 879` periods with at least
two distinct primes -- `173 537` proved, `8 338` block-saves, `4` band members
(`30`, `858`, `1722`, `66198`, the four Giuga numbers below `200 000`) -- and the
gap is zero on every one of them, with the drop-one block as the witness on
every one of them. `--limit` moves the horizon; the report carries the frontier
with the number so the two cannot be cited apart.

#### The three regimes are one integer predicate

`#odd_B` is decided by the two-adic valuation alone -- `w` when `L` is odd, one
when `v2(L) = 1`, none when `4` divides `L` -- so the whole sign reading is

    block_saving(L) = (M - slack(L)) // 2,      M = L * (e1 - 1) in Z,
    band  <=>  1 <= M <= slack(L),

with `slack = 0` when `4 | L`, `1` when `v2(L) = 1`, and `w - 1` when `L` is
odd.

*The law is a congruence, not a census.* `M = sum_{p|L} L/p - L`, read modulo
two, counts exactly the primes whose cofactor is odd, minus whether `L` is odd --
which is `slack` by definition:

    M  ===  #{p : L/p odd}  -  [L odd]  (mod 2)  =  slack(L)  (mod 2),

so `M - slack` is always even, the division in `block_saving` never drops a
unit, and the statement holds for **every** `L` with `w >= 2` instead of for the
`2000` periods it was measured over. `excess_parity(L)` is that residue;
`tests/test_band_parity.py` holds it against `prime_excess` on every period to
`2000` and against `block_saving` on the same population. The `w = 1` corner is
where it fails (`L = 4` has `M = -2` against `slack = 0`), so the function
refuses a prime period rather than returning a number the law cannot use.

One corollary reaches the open problem. On an odd period `slack = w - 1`, hence
`M === w - 1 (mod 2)`: an odd Giuga number would need `M = 1`, and `M = 1`
forces an **even** number of prime factors. No odd Giuga number is known; this
does not shorten the search, but it turns the corner from a point into a parity
class. It is also why the even band is a single point -- `v2(L) = 1` gives
`slack = 1`, so the band is `M = 1` and nothing else.

The odd row is the one the proof could not close, because `M = 1` there *is* an
odd Giuga number -- `sum(1/p) - prod(1/p) = 1` -- open since 1950. It is no
longer an unexplained residue either: the band is a window, and the window has a
floor. Eight odd primes are not enough to cross one (`111429982/111546435`, short
by `116453/111546435`), so the first candidate sits at `w = 9`:

```
$ exact-lattice-quantiles odd-band --max-omega 10 --max-prime 100
odd band: 1 <= period * (sum 1/p - 1) <= omega - 1
  odd primes available                    24  at most 100
  prefixes walked                      13866  minimal feasible sets only
  threshold                         omega - 1  omega - 1 by law
  feasibility floor                 omega = 9  product 3234846615 clears the window by 108169298
  members in the band                      0
  nearest miss                         99604  clears the band by 99596 at omega = 9
  3 5 7 11 13 17 19 43 47              99604  the closest approach from above
  omega 9                               99604  by 99596 over 8
  omega 10                            6669349  by 6669340 over 9
  parity of the excess              consistent  for an odd period, prime_excess(n) === omega - 1 (mod 2)
  boundary: an odd member with excess 1 would be an odd Giuga number, sum(1/p) - prod(1/p) = 1, open since 1950 -- and the parity of the excess puts that corner at even omega
  verdict: no odd period in the band over 13866 prefixes: at most 10 odd primes, each at most 100
```

The floor and the margin are different optima, and the walk now answers both.
The floor minimises the **product** that crosses one; the margin minimises the
**excess** and reports the nearest miss with the primes that produce it --
`{3, 5, 7, 11, 13, 17, 19, 43, 47}`, `L = 9801536745`, excess `99604` against a
slack of `8`: the band is missed by `99596`, eight thousand times closer than the
floor witness whose own product is three times smaller. `margin` carries that row
for every `w`, and `parity` re-checks the congruence on each one, so an empty
answer comes back as a distance with a witness instead of a silence.

`13866` is what verifies, not what it would take to search: the walk enumerates
only *minimally feasible* sets, because `M(P + q) = q * M(P) + L(P)` with
`M({}) = -1` means every proper subset of a member has `M <= 0` -- otherwise the
member would need an excess above `w - 1` for `w < 19686`, which the recursion
forbids. So the search of an open question stays finite, and it is finite at
`w = 12` and primes to `200` too: `8410793` prefixes in `31.5 s`, no member --
and at `w = 13` with primes to `250`: `102422466` prefixes in `700.7 s`, still
none.

An empty answer from a pruned walk proves nothing unless the pruning is blind
-proofed, so the same walk runs with the threshold lifted to `10**9`: it returns
the floor witness itself, `{3, 5, 7, 11, 13, 17, 19, 23}`, excess `108169298`,
in `0.03 s`. `--slack` is that control, and `tests/test_odd_band.py` also
enumerates the band exhaustively by combinations over a smaller horizon and gets
the same zero. The instrument returns its own refutation when one exists.

What is *not* claimed: that the odd band is empty. It is empty below a horizon
(`w <= 12`, primes `<= 200`), and a member above it is an odd Giuga number. That
is why `sign_regime` publishes `characterisation: "giuga-window"` for odd
periods and `"giuga-equation"` only for the three that solve it -- and why the
version stays below `1.0.0`.

`1581` of the `1914` proved periods are *sharp*: the heaviest antichain weighs
exactly the block, so the bound is tight and a proof cannot be pushed further
without extra structure. In `858` the heaviest is `429 + 286 + 78 + 66 = 859`
against `L * e1 = 859`, and the odd-count slack absorbs the extra point -- the
sharpness is what keeps the band from being empty by accident.

The reciprocal sum is not the criterion either. `N = 30` has `1/2 + 1/3 + 1/5 >
1` and saves nothing, because `7 + 5 + 3 = 15` buys exactly the fifteen points
the single grid `30` costs; `60` and `90` are the same shape with a square
factor and save one point each.

#### Why the price is answerable where the table is not

`price_exact` runs its DP over the **maximal elements of the request**, not over
the universe: a grid that covers `b` covers every divisor of `b` for free, so a
request costs exactly what its maximals cost, and the state space is `2^k` with
`k` the number of maximals. Verified against the `2^|U|` table on all 8 388 607
requests over the divisors of 360, with 0 discrepancies -- and it is what makes
`price_exact(1080, (54, 90, 135))` return `135` in under a millisecond where the
full table is `2^23` masks, and `2520` is `2^32` — a census the package refuses
by name (`RC 5`, `MaskTableTooLarge`) instead of trying to allocate. The
contrast is the whole argument: `price_form_census(2520)` enumerates all `2^32`
requests and is refused, while `price_exact(2520, U(2520))` returns `649` in
under a millisecond, because those 32 divisors have 6 maximals and the request
costs what its frontier costs.

So the library never prices with the formula: `price_exact` is a subset DP over
the coverage masks, and `price_form_census` publishes how often the shortcut
agrees. Restricting the universe to `U(T)` is what makes the shortcut look
right, and that restriction is a modelling decision, not a measurement --
`price_form_census(period, reduced=False)` is the row that shows it.

#### The formula also assumes no fixed cost per grid

Every census above prices a grid at `grid // 2`. That is a hypothesis, not a
fact about billing, and the formula does not state it. Charge a fixed cost
`setup` per grid bought and the closed form -- which buys one grid per maximal
divisor and never buys a block -- stops being the price at a fee of **one**:

```
$ exact-lattice-quantiles setup-fee
setup fee: the closed form assumes a grid costs grid // 2 and nothing else
  period                180           reduced universe, |U| 12, 4095 requests
  fee                   discrepancies pairs  witness
    0                   0                 0  -
    1                   32                0  (6, 10, 15) by 2
    2                   146               2  (6, 10, 15) by 4
    3                   237               4  (6, 10, 15) by 6
    4                   1697              8  (6, 10, 15) by 8
  frontier              1             first fee per grid that breaks the closed form
  pairs                 2             first fee at which a pair is what fails
  theorem  advantage(B, setup) = setup * (len(B) - 1) + sum(b // 2) - lcm(B) // 2; an incomparable pair sits at -1 or below on a zero fee, so no fee below 2 lets a pair win
  verdict  the closed form stops being the price at a fee of 1 per grid: 32 of 4095 requests, witness (6, 10, 15) costs 16 exact against 18 by the formula
```

One identity decides the whole column. Folding a block `B` into its least
common multiple refunds one fee per folded grid and pays the rest with the
parts' weight against the lcm's:

```
advantage(B, setup) = setup * (len(B) - 1) + sum(b // 2) - lcm(B) // 2
```

Two consequences, and the second is why the frontier has to be measured over
every request of the instance rather than argued from pairs:

- **The fee favours the big block.** A triple folds two fees, so the requests
  that break first are triples, not pairs. The witness `(6, 10, 15)` has
  `advantage = 0` at a zero fee -- it *tied* -- and the two refunded fees turn
  the tie into a 2-point saving. That is why its gap grows by exactly
  `len(B) - 1 = 2` per unit of fee, `2 4 6 8` in the column above.
- **No pair can win below a fee of two.** For an incomparable pair
  `lcm(a, b) >= a + b` makes `sum(b // 2) - lcm // 2` at most `-1`, so the
  advantage stays at `0` or below for a fee of one, and `(2, 3)` attains the
  bound exactly. A census that only looked at pairs would report the closed
  form surviving one fee further than it does.

On the divisors universe the formula is already broken at a zero fee, so the
frontier is `0` and the witness is the block that breaks it there:

```
$ exact-lattice-quantiles setup-fee --universe divisors --max-setup 1
    0                   9984              0  (36, 60, 90) by 3
    1                   10016             0  (36, 60, 90) by 5
  frontier              0             first fee per grid that breaks the closed form
  pairs                 none          first fee at which a pair is what fails
```

And when the horizon never reaches a break, the command refuses to report a
frontier it did not find -- `RC 4`, with the reason on stderr:

```
$ exact-lattice-quantiles setup-fee --max-setup 0 ; echo $?
  verdict  NO_MEDIDO: no fee up to 0 breaks the closed form on period 180; the frontier is beyond the horizon, not absent
4
```

A note on what this does *not* say: it does not say the exact price is wrong.
`price_exact(period, request, setup=s)` prices every fee exactly, and the
identity above is what a design should be compared against. It says the
**formula** has a domain, the domain is `setup = 0`, and the package now
publishes where it ends instead of assuming nobody bills a setup.

#### The frontier is a law, and only a triple can be priced at one fee

The `setup-fee` command walks the fee upwards on **one** period. The same
question has a closed answer, and it is cheaper to compute than the walk: a
request can only break the formula through an antichain -- the formula prices
the *maximals* of a request, so a part that divides another is never priced
apart -- and the fee at which an antichain `A` breaks it is `block_fee(A)`,
which is `block_advantage` inverted in one integer division:

```
$ exact-lattice-quantiles fee-census
fee census: the fee that breaks the closed form, period by period
  horizon                 2000            U(T) universe, blocks of at most 3 parts
  periods                 1000            909 verified against every antichain (|U| <= 14)
  mismatches              0               periods where no block of that size attains the frontier
  chains                  10              no antichain at all: no fee breaks these
  at zero                 1               already broken at no fee (1080,)
  at one                  51              break first at a fee of 1
  widest frontier         499             period 1994 on (2, 997)
  witness sizes           2: 889, 3: 10   on the verified periods
  pairs miss              52              a triple beats every pair
  pair over               0               of those, above a fee of one: the law refuted if not 0
  pair law 52 periods put a triple above every pair and all of them sit at a fee of at most one, so the third part buys exactly one thing: the fee of one
  verdict  the law holds on all 909 periods of U(T) with at most 14 grids, witnesses of size [2, 3]; 10 periods are chains where no fee breaks the formula, 1 break with no fee at all [1080], and 51 break first at a fee of 1 on 2 shapes (6, 10, 15), (6, 14, 21); the widest frontier is 499 at period 1994 on (2, 997)
```

Four statements come out of it, and each one is a place the census could have
failed:

- **A fee of one needs a triple, and there are only two shapes.** In the whole
  horizon only `(6, 10, 15)` and `(6, 14, 21)` break first at a fee of one --
  `{2p, 2q, pq}` with `p = 3` and `q` equal to `5` or `7`. Both were already the
  witnesses of the walk on period 180.
- **A fee of two is what a pair costs, and `(2, 3)` attains it.** That is the
  pair theorem of the section above, used here as the floor of the census
  instead of as a footnote: no period breaks below two without a triple.
- **Ten periods of a thousand are chains no fee can break**, and below 2000 they
  are exactly the powers of two. That half is argued rather than searched: the
  divisors of `2^k` are totally ordered, so `U(2^k)` has no antichain at all.
  The census answers `chain` -- not a zero, which would claim a frontier of no
  fee.
- **`1080` breaks with no fee at all**, which is the break the carrier-block law
  found from the other side. Two instruments, one number.

The law is a minimum over antichains of at most three parts, which is cheap; the
*check* is exhaustive over antichains of any size, which is not. So the census
publishes `909` verified periods rather than pretending to `1 000`, and the check
is the reason to trust the law outside them. Restricting the law to pairs is
enough to refute it, and the refutation names exactly the periods a triple is
needed for:

```
$ exact-lattice-quantiles fee-census --max-period 200 --arity 2 ; echo $?
REFUTED: the block law misprices 5 of 100 periods up to 200
fee census: the fee that breaks the closed form, period by period
  horizon                 200             U(T) universe, blocks of at most 2 parts
  periods                 100             100 verified against every antichain (|U| <= 14)
  mismatches              5               periods where no block of that size attains the frontier
    period 60             1               law says 2 by (6, 10, 15)
    period 84             1               law says 2 by (6, 14, 21)
    period 120            1               law says 2 by (6, 10, 15)
    period 168            1               law says 2 by (6, 14, 21)
    period 180            1               law says 2 by (6, 10, 15)
  chains                  7               no antichain at all: no fee breaks these
  at zero                 0               none in this horizon
  at one                  0               break first at a fee of 1
  widest frontier         49              period 194 on (2, 97)
  witness sizes           2: 88           on the verified periods
  pair law NOT_APPLICABLE: at 2 parts the law is already the pairs, so no triple is enumerated and nothing here measures what one buys
  verdict  REFUTED: the law prices period 60 at fee 2 while the exhaustive frontier is 1 -- a block of more than 2 parts wins earlier than every block of at most 2
1
```

`RC 1` is the package disagreeing with itself, and it exits on the periods where
the third part is what wins.

And the third part wins for exactly one reason, which is the sharper statement of
the law: **a triple only ever buys a fee of one.** Enumerate pairs and you get the
frontier everywhere the fee reaches two; the `52` periods where you do not are
the `51` at a fee of one plus the `1` at no fee at all, and nothing else:

```
$ exact-lattice-quantiles fee-census --max-period 400
fee census: the fee that breaks the closed form, period by period
  horizon                 400             U(T) universe, blocks of at most 3 parts
  periods                 200             198 verified against every antichain (|U| <= 14)
  mismatches              0               periods where no block of that size attains the frontier
  chains                  8               no antichain at all: no fee breaks these
  at zero                 0               none in this horizon
  at one                  10              break first at a fee of 1
  widest frontier         100             period 398 on (2, 199)
  witness sizes           2: 181, 3: 9    on the verified periods
  pairs miss              10              a triple beats every pair
  pair over               0               of those, above a fee of one: the law refuted if not 0
  pair law 10 periods put a triple above every pair and all of them sit at a fee of at most one, so the third part buys exactly one thing: the fee of one
  verdict  the law holds on all 198 periods of U(T) with at most 14 grids, witnesses of size [2, 3]; 8 periods are chains where no fee breaks the formula, 0 break with no fee at all, and 10 break first at a fee of 1 on 2 shapes (6, 10, 15), (6, 14, 21); the widest frontier is 100 at period 398 on (2, 199)
```

That is what makes the census affordable instead of merely correct. The law over
antichains of at most three parts costs `C(|U|, 3)`; the law over pairs costs
`C(|U|, 2)`, and the periods it cannot answer are the ones whose answer is `0` or
`1` -- which a triple supplies. So the enumeration is a search for the two shapes
of the fee of one, not a search for the frontier.

### 5. A margin survives exactly zero drift

The witness margin is a knife edge: the divisor that fixes it sits on the edge
of its block, so **any** common drift loses it. That is a property of the
optimum, not a slack in the implementation — the sections are placed to just
seat every divisor, and the last one seated is the one that leaves no room.

`drift` prints both halves of the answer. Two quantities, not one:

- **coverage drift** — the largest common `delta` the sections still *cover*
  every divisor under. A section may be pushed into a block it never captured,
  so this is the weaker question.
- **margin kept** — the largest common `delta` the sections still hold the
  witness margin under. This is the guarantee.

```
$ exact-lattice-quantiles drift 180 --max-cost 64 --probes 4
period           180
design           [12, 18, 20, 45]
cost             47 points (budget 64)
sections (4)  ['3/10', '17/30', '67/90', '29/30']
witness margin   1/15
coverage drift   1/15  (captured, margin not kept)
margin kept      drift tolerated
  1/15            0
  1/20            1/60
  1/30            1/30
  1/60            1/20
  0               1/15
```

Read the first row and the last: the margin is held under **zero** drift, and
the full margin is only available at zero drift. In between the curve is the
line `W - v` — the margin you can *demand* costs exactly itself in tolerance,
so margins and drift tolerance are the same currency spent twice. The design
that carries `1/6` on the same period behaves identically
(`drift 180 --max-cost 47 --probes 8`: `1/6 -> 0`, `1/12 -> 1/12`).

Coverage is a genuinely weaker measure, and the gap is not small. This design
holds `1/6` and still covers under a drift of `7/30`:

```python
from exact_lattice_quantiles import (best_probes, drift_radius, universe,
                                    witness_margin)

G = (2, 3, 4, 6, 9, 12, 15, 20, 36, 45, 180)
_, Q = best_probes(G, universe(180), 4)
min(witness_margin(G, Q, universe(180)).values())          # 1/6
min(drift_radius(G, Q, universe(180)).values())            # 7/30
```

`7/30` is the exact radius of the connected component of the union of
admissible drift spans that contains zero. An earlier version of
`drift_radius` answered `1/12` here: it only counted sections that captured a
divisor at `delta = 0`, and drift can capture. The under-report is asserted in
`tests/test_coverage.py` against the old formula, so it cannot come back.

## What the exit codes mean

The distinction is the point: an infeasible margin cannot be bought, an
over-budget one can.

| RC | exception | meaning |
| --- | --- | --- |
| 0 | — | a schedule exists and is returned |
| 1 | `InfeasibleAboveInstanceCeiling` | above `nu(T)`: no design at **any** cost |
| 2 | `InfeasibleWithinBudget` | exists, but `cheapest > max_cost` |
| 3 | `ValueError` (CLI) | malformed request — bad fraction, bad flag, odd period |
| 4 | — (CLI) | the verdict is not `OK`: no schedule was delivered |
| 5 | `MaskTableTooLarge` | the answer exists; its table is above `TABLE_BUDGET` |

`RC 3` covers two different refusals and names which. A malformed request is a
`ValueError` about the *syntax*; an odd period is `OddPeriodUnsupported` (a
`ValueError` subclass, so a caller that caught the old one is unaffected) about
the *subject*: `U(T)` is empty, so there is nothing to compute. Before 0.9.1
the second arrived as CPython's `ValueError: min() iterable argument is empty`
— the arithmetic reporting on itself, with the period nowhere in the message.

`RC 5` is not `RC 3`, and the distance between them is the whole distinction: an
odd period has no answer to give, a period above the table budget has an answer
and no room to compute it here. Six call sites size a bitmask table by `2^n` in
the number of divisors -- the cheapest-design DP, the set-cover DP, the design
census, the maximal-divisor table, the cheapest-cover table and the frontier
price table -- and before 0.11.0 the only bound on that allocation was the
machine, so the process was killed with `RC 137` and the caller learned nothing:
not the period, not the site, not that the schedule exists. From 0.11.0
`mask_table` is the only way a site sizes a table, and it refuses by name before
a single list is built:

```
$ exact-lattice-quantiles solve 4032 --min-margin 1/6     # RC 5
period 4032: the cheapest-design DP would allocate 2^26 masks over 12 rows
= 805306368 cells, above the budget of 536870912; the answer exists and this
implementation will not allocate the table
```

#### The cheapest-design table was a property of the formulation

From 0.12.0 the design DP does not count divisors, it counts **classes of
divisors that are not already implied**. Write `cols[d]` for the grids whose
block for `d` is at least `v` wide. A design holds `d` exactly when it contains
a grid of `cols[d]`, so if `cols[e]` is contained in `cols[d]`, every design
that holds `e` also holds `d`: the constraint on `d` is implied and dropping it
cannot move the optimum. The divisors that survive that test are the **kernel**
of the coverage matrix, and the DP runs over `2^|K|` states instead of
`2^|U|`, where `K` is the kernel:

| period | divisors | kernel at `1/6` | cells | cost at `1/6` | time |
| --- | --- | --- | --- | --- | --- |
| 720 | 22 | 14 | `12 * 2^14` = 197 K | 271 | 0.0 s |
| 1440 | 27 | 19 | `12 * 2^19` = 6.3 M | 590 | 1.3 s |
| 1680 | 27 | 19 | `12 * 2^19` = 6.3 M | 652 | 1.3 s |
| 1800 | 24 | 13 | `12 * 2^13` = 98 K | 499 | 0.0 s |
| 2520 | 32 | 20 | `12 * 2^20` = 12.6 M | 888 | 2.9 s |
| 2880 | 32 | 24 | `12 * 2^24` = 201 M | 1227 | 51 s |

`T = 1440`, `1680` and `2520` were refused by 0.11.0 -- `2520` asked for
`2^32` masks -- and are now computed. The budget did not move: the table shrank,
from `2^32` to `2^20` on the period that was four thousand times over it. The
kernel is a property of the *threshold*, not of the period: `2520` keeps 18 of
its 32 divisors at `1/126`, 20 at `1/6` and 23 at its own ceiling of `8/35`,
because near the ceiling almost nothing is implied any more. Above the ceiling
every column empties, all of them are equal, and the kernel is one empty column.

The equivalence is pinned twice, not argued: against an enumeration of all
`2^12` designs at every one of the 22 widths of `T = 180`, and against the
previous `2^|U|` DP on **2377 instances** -- every even period from 2 to 400
against every width of the instance *and* the coverage-only threshold `0`, with
0 discrepancies. Threshold `0` is the case that a defaulted read of a
non-existent pair gets wrong: a pair that does not exist is not a wide block,
and reading it as zero width makes every grid cover every divisor and the
cheapest design collapse to a single grid. `margin 0` is the coverage price,
`47` points for `T = 180`.

Three sites still refuse, and they refuse for a reason the kernel does not
touch:

| site | period | cells | how a caller meets it |
| --- | --- | --- | --- |
| design | 4032 | `12 * 2^26` = 805 M | through the API: at `1/6` its kernel keeps 26 of 33 divisors |
| probes | 4032 | `2 * 2^33` = 17 G | `solve` refuses here first, before the design DP |
| census | 720 | `29 * 2^22` = 122 M | `designs 720` |

So the design DP is no longer the binding refusal anywhere a caller can reach
it: the set-cover table grows with the universe and not with the kernel, so for
`|U| >= 29` it refuses before the design DP is sized, and below that the kernel
fits. Measured on every even period up to 12000 with `|U| <= 28`: the kernel
never exceeds 25 divisors at `1/6` or `1/5`, and 25 is exactly what the budget
allows.

`TABLE_BUDGET = 2^29` is a **cell count** -- slots touched -- and not seconds,
because wall time on a shared machine moved by 4x between two runs of the same
instance. Peak RSS read with `/usr/bin/time -l`; `getrusage` inside the process
reported `0.54 GB` right after a `1.2 GB` allocation and is not an instrument
here. And the refusal eats nothing that was ever answered: `T = 360` at a `1/6`
margin is still `112` points, and that is a test, not a promise.

`PROBES_INSUFFICIENT` is a verdict on the returned `Schedule`, not an
exception: the design was found, the sections were not enough, and the object
says by how much (`probes_needed`). The CLI prints the object and exits `RC 4`
— a caller in a shell must not read "we have a schedule" from a schedule that
was not delivered.

`UNDERCOVERED` is the third verdict and the one that is not about the margin at
all: the sections placed do not land inside every divisor of the universe, so
the margins that would be reported are ``0`` for a divisor nobody holds — an
absent measurement wearing the value of a measured zero. From 0.9.0 such a
schedule carries `verdict = "UNDERCOVERED"` and the divisors it missed in
`uncovered`, and **both `witness_margin` and `robust_margin` are `None`**, in
the object and in `as_dict()`. The CLI prints `NO_MEDIDO` for them and exits
`RC 4`. Before 0.9.0 the same object reported `0` and, when the demanded margin
was `0`, called itself `OK` with an empty section set; a caller that reads a
margin must now read `None` first.

## Four margins, and which one to quote

Four quantities, and they answer four different questions:

- `design_margin` — what the grid set guarantees before any section is placed.
- `witness_margin` — what the chosen sections actually hold, with each divisor
  leaning on its own grid/section pair. This is what a schedule delivers, and
  it never exceeds `design_margin`: the four quantities nest.
- `drift_radius` (`Schedule.robust_margin`) — the largest common drift `delta`
  the sections still **cover** every divisor under. Exact for coverage; it is
  *not* optimistic, it is a weaker question: a section may be pushed into a
  block it never captured, so this can exceed the witness margin (measured
  above: `7/30` against `1/6`).
- `drift_tolerance` — the largest common drift `delta` the sections still hold
  a demanded margin under. `drift_tolerance(..., 0) == drift_radius(...)`, and
  `drift_tolerance(..., witness_margin)` is `0` on every design measured.

Coverage and slack are not the same predicate, and neither is the other's
bound.

The hierarchy `design_margin >= witness_margin >= drift_radius` holds at all
200 even periods in `2..400`, and so does `design_margin <= instance_ceiling(T)`.
It did not always: before 0.3.0 the block map divided by the number of sections
(`M // 2`) while `cell_width` divided by half the period (`M / 2`), and the two
agree only on even grids, so every odd grid inflated the witness by `M / (M-1)`
and 96 periods reported sections holding more than their own design guaranteed
([issue #6](https://github.com/alexisromtz-lab/exact-lattice-quantiles/issues/6)).
One normalisation in `_interval_map` closed it; the census that pins it is
[tests/test_solver_contract.py](tests/test_solver_contract.py).

Two figures changed with the repair, both in the direction of the geometry:
`k = 4` on the cheap design holds `1/15` rather than the inflated `2/33`, and
the coverage radius of the design above is `7/30` rather than `37/210`.

## The price of a margin is one-for-one, all the way down

A **design** buys margin with drift tolerance at a price of exactly one for one,
across its entire range:

```
tolerance(v) = radius - v        for every 0 <= v <= radius
```

On `T = 180` at `k = 4` the sections cover a common drift of `1/3` at zero
margin, and every further tenth of margin costs a tenth of tolerance:

```
$ exact-lattice-quantiles drift 180 --max-cost 64 --probes 4 --min-margin 1/3
coverage drift   1/3  (captured, margin not kept)
tolerance(v)     1/3 - v for every v up to 1/3  (one-for-one, no cliff)
narrowest        divisor regime ends at 1/30  (pessimistic: the design pays nothing there)
margin kept      drift tolerated
  1/3             0
  1/4             1/12
  1/6             1/6
  1/12            1/4
  0               1/3
```

All five rows lie on `1/3 - v`. That ``narrowest`` line is the interesting one.
Asked *per divisor*, `drift_law` reports a critical margin and a cliff after it:

```
tolerance(v) = radius - v        for every v < critical_margin
tolerance(v) < radius - v        at  v = critical_margin
```

and the divisors of this design do break at `1/30`. The design does not,
because the tolerance of a design is the **minimum** over its divisors: a
minimum of lines of slope `-1` is a line of slope `-1`, and a divisor that falls
under its own line only wins that minimum when it also holds the radius. The
divisor that breaks first is not the divisor that binds — so a per-divisor
reading is pessimistic about the design, and a design has no cliff at all.

Three measurements, all reproducible:

- per divisor, strictly below its edge, the line holds in **1680 of 1680**
  readings over the 200 even periods at `k = 4` and `k = 6`; at the edge the
  tolerance is under the line in **369 of 369**, never on it;
- per design, the line holds on every interval between consecutive knots of the
  span geometry in **3076 of 3076** designs — 200 even periods, budgets
  `{16, 24, 32, 47, 64, 96}`, `k` in `{2, 3, 4, 5, 6, 7, 8}`;
- a design's line rests on one inequality per divisor — a drop of `cliff_d` from
  `radius_d` only undercuts it while `cliff_d > radius_d - radius` at a critical
  margin inside the range — and that holds in the same **3076 of 3076**. It is
  not slack: in **1234** readings the divisor *does* have a cliff to jump and
  lands exactly **on** the line. The smallest case is `T = 12`, design `(4, 6)`,
  divisor `4`: it breaks at `1/6` and drops `1/3` from a radius of `2/3` — exactly
  the design's excess — so its tolerance meets the line at `v = 1/6` and no lower;
- for **78 of 240** designs at `k = 4` and `k = 6` the first divisor break is
  strictly smaller than the radius, and the design pays nothing at it.

[tests/test_drift_law.py](tests/test_drift_law.py) pins the per-divisor half —
the second check fails if `_linear_regime` is replaced by "the line holds
everywhere" — and [tests/test_design_linearity.py](tests/test_design_linearity.py)
pins the design half.

## Drift is not a second budget

The price above is paid in one currency. A design that holds a margin `v` under
a common drift `delta` is quoted as one number, `v + delta`, and 0.17.0 reads
the exchange in both directions:

```python
guaranteed_margin(G, Q, U, delta)    # the margin the sections still keep
margin_under_drift(G, Q, U, delta)   # the same, per divisor
```

`guaranteed_margin` is the inverse of `drift_tolerance`, and it is how the
fourth margin above should be read: drift buys nothing on its own, it is spent
from the witness margin. On `T = 180` the cheapest design for `1/30` holds
`1/15` with no drift at all, so a `1/60` drift leaves it half of that:

```
$ exact-lattice-quantiles drift 180 --max-cost 64 --probes 4 --min-margin 1/30 --delta 1/60
period           180
design           [12, 18, 20, 45]
cost             47 points (budget 64)
sections (4)  ['1/30', '17/30', '7/10', '83/90']
witness margin   1/15
coverage drift   1/15  (captured, margin not kept)
tolerance(v)     1/15 - v for every v up to 1/15  (one-for-one, no cliff)
narrowest        divisor regime ends at 2/45  (pessimistic: the design pays nothing there)
margin kept      drift tolerated
  1/15            0
  1/20            1/60
  1/30            1/30
  1/60            1/20
  0               1/15
under drift      1/60  ->  keeps margin 1/20
budget           one currency: keeping 1/20 under that drift prices the same as keeping 1/15 with none
```

Read off the table instead of the law: `1/20` kept at `1/60` drift is the same
row as `1/15` kept at no drift, one budget spent either way.

Two measurements, one over the 60 even periods `2..120`, both budgets
`{the cover price, 64}`, `k` in `{2, 4}`, `delta` in `{0, 1/120, 1/60, 1/30}` —
**641 readings**:

- `guaranteed_margin(design, Q, U, delta) == witness_margin - delta` in
  **641 of 641**. No exception, no cliff: the law is not a fit, it is the
  arithmetic of a design bought at its cover price;
- the inverse round-trips in the same **641 of 641**:
  `drift_tolerance(..., guaranteed_margin(..., delta)) == delta`.

The exchange rate is the **witness margin**, not the coverage radius, and the
two part company as soon as a design is not bought at its cheapest cost. In
**354 of 1800** random section sets the radius exceeds the witness margin, and
what the radius adds is not margin:

```
T = 12, grids (4, 6), sections (1/6, 1/4, 1/2, 7/12, 11/12)
witness margin   1/6
coverage drift   1/3
tolerance(1/12)  1/4        <- the line is r - v here
tolerance(1/6)   0          <- divisor 2 breaks exactly at the witness margin
guaranteed_margin(delta = 1/100)   = 1/6, not 1/3 - 1/100
guaranteed_margin(delta = 1/4)     = 1/12 = 1/3 - 1/4
```

The sixth above `1/6` is reach the sections have and margin they never held. So
`drift_radius` is what the sections **cover**, never what they **keep**, and a
caller who quotes it as a guarantee quotes a number the schedule cannot hold.

That the radius never falls below the witness margin is a theorem, not a
census. `drift_radius` is the largest margin whose tolerance is still
non-negative, so it cannot sit below the margin the design already holds: a
shrunk span holds a *strict* seed at zero only when its pair's slack exceeds
`v`, and at `v == w` the divisor that attains `w` has no pair above it, so
`tolerance(w) == 0` and stays `0` for every larger `v`. Two consequences, and
the sweep in [tests/test_drift_budget.py](tests/test_drift_budget.py) checks
both over **186 designs** of `T` in `{12, 18, 20}` — of which **49** do carry a
radius above their witness margin:

- `radius >= witness_margin` in **186 of 186**. An exhaustive census says the
  same over every design with up to 3 sections of `T` in `{12, 18, 20, 24, 30}`:
  **0 of 951508** below, **58712** above — the gap is the common case, not an
  accident of one period;
- `guaranteed_margin(delta) <= witness_margin` for every `delta > 0`, so a
  service level `keep` survives a drift `delta` only while `keep + delta` is
  inside the instance ceiling. The drift is paid out of the ceiling.

Which makes a service level one command and not two budgets:

```
$ exact-lattice-quantiles plan 180 --max-cost 64 --probes 5 --keep 1/30 --delta 1/60
period           180
service level    keep 1/30 under a drift of 1/60  ->  buy margin 1/20
design           [12, 18, 20, 45]
cost             47 points (budget 64)
sections (5)  ['1/18', '1/2', '53/90', '13/18', '8/9']
witness margin   1/9
under drift      1/60  ->  keeps margin 17/180
verdict          OK
```

`plan` is the plain solve at `keep + delta` and then re-reads its own answer
with `guaranteed_margin` before printing: `0` the level is met, `1` the level
plus the drift is above the instance ceiling (no schedule at any cost), `2` the
budget is short, `4` `NO_MEDIDO` — the sections do not realise the level asked
for, and a shortfall is reported instead of a number.

One reading rule, because it changes the answer: the value returned is a
**supremum**. Where a divisor breaks, the margin it names is approached, not
attained — `drift_tolerance` at that value is `0`. The guarantee is therefore
strict: every margin strictly below the answer tolerates the drift. Reading an
attained value instead would understate the answer by an arbitrary amount; on
the design above at `delta = 1/6` it reports `1/12` where the sections keep
`1/6`.

[tests/test_drift_budget.py](tests/test_drift_budget.py) pins the law, the round-trip,
the gap, and the theorem that closes it.

## A clinic is a period

A chair offering 15, 20, 30, 45 and 60 minute services repeats its booking
grammar every `lcm(15, 20, 30, 45, 60) = 180` minutes, and the candidate cadences
of that grammar are the divisors of 180 above one — 17 of them. The lattice of
`T = 180` **is** the lattice of that chair, so the arithmetic above answers a
scheduling question with no scheduler in it:

```
$ PYTHONPATH=src python examples/clinic_period.py
durations            15, 20, 30, 45, 60 minutes
period               180 minutes  (lcm of the durations, 17 candidate cadences)
instance ceiling     1/3 = 60 minute(s)  (critical divisor 6)
cheapest design      47 points: cadences [12, 18, 20, 45]
at the ceiling       64 points: cadences [4, 6, 10, 12, 15, 18, 20, 45]  ->  holds 1/3 = 60 minute(s)
ladder of period 180
     47 points -> margin 1/6 = 30 minute(s)
     52 points -> margin 8/45 = 32 minute(s)
     59 points -> margin 1/5 = 36 minute(s)
     64 points -> margin 1/3 = 60 minute(s)
service level        keep 1/30 = 6 minute(s) under a drift of 1/60 = 3 minute(s)
  bought as          margin 1/20, 47 points, 5 shared sections
  still keeps        17/180 = 17 minute(s)
  ceiling check      keep + drift = 1/20 <= 1/3: True
```

The same bridge from a schedule to a period is one command, so the ceiling can
be asked for the lengths a system actually runs:

```
$ exact-lattice-quantiles period 15 20 30 45 60
durations         15, 20, 30, 45, 60
period            180  (lcm of the durations)
grids             17 divisors of 180 above one
universe U(T)     12 carriers of the even grids
ceiling nu(T)     1/3  = 60 of the period
critical divisor  6
radical bound hit False
services per period  15: 12, 20: 9, 30: 6, 45: 4, 60: 3
```

An odd period is refused by name rather than answered: every divisor of an odd
`T` is odd, so `U(T)` is empty and there is no ceiling to report — `rc 3`, never
a number. That is the same boundary `ceiling` already keeps.

Every number in that transcript is computed by
[examples/clinic_period.py](examples/clinic_period.py) and pinned by
[tests/test_examples.py](tests/test_examples.py); none is quoted by hand.

Read it as a reading of the arithmetic, not as a shipped scheduler. What the
library prices is the margin that survives a delay hitting *every* service at
once — a late opening, a cleaning cycle running long — because that is the
failure an average slot length cannot see: padding each service separately pays
for the same delay once per service, and a shared buffer collapses the moment
they all slip together. The clinic schedules itself; what it cannot compute by
eye is where its own ceiling is, and `keep + drift <= 1/3` is that answer for a
chair with those five lengths.

## Exactness

Everything the library computes and returns is `int` or
`fractions.Fraction`. No comparison is decided by rounding and no optimum can
be `1/3 - 1e-16`. `as_dict()` serialises every rational as `"n/d"`; `--json` on
the CLI does the same. The single float in the package is the human-readable
percentage the CLI prints next to the exact fraction, and the two reciprocal
sums `strong_lemma_frontier` reports for a human to read beside the exact
integers that decide. No comparison is ever made on either: the crossing is
located with `Fraction` and confirmed with `int`, floors included.

## Verification

```
python -m unittest discover -s tests -t .     # 397 tests
exact-lattice-quantiles --self-test           # 102 checks, RC 0
```

The suite anchors each number quoted above to the function that produces it,
including the cases where the first implementation was wrong: filtering the
scan to even grids halves the ceiling of every odd divisor; `T = 2` needs its
degenerate branch; a design that does not seat every divisor has margin `0`
rather than raising; and the coverage DP is compared against brute force.

### A claim kept outside the package can be recomputed

Every number above gets quoted somewhere else — in notes, tickets, drafts — and
nothing recomputes those quotes. So the quote and the library drift apart while
both look fine, and the only witness to the boundary is somebody remembering
it. `verify-acta` turns that boundary into a gate: give it a **manifest**, a
JSON list of anchors each stating one quantity for one period, and it
recomputes every one in exact arithmetic and reports **per field**, so a
divergence names the quantity instead of the file.

```
exact-lattice-quantiles verify-acta \
  src/exact_lattice_quantiles/manifests/acta_margen_secciones_20261125.json
# 90 checks over 30 anchors: 90 reproduce, 0 divergent        RC 0
```

An anchor is one of three kinds:

| kind | required | optional |
| --- | --- | --- |
| `ceiling` | `period`, `margin` | `critical_divisor`, `grids`, `universe_size`, `cost_cover`, `cost_ceiling` |
| `design` | `period`, `max_cost`, `probes`, `margin` | `cost`, `design`, `probes_needed`, `robust_margin`, `verdict` |
| `ladder` | `period`, `design`, `probes`, `margin` | — |
| `frontier` | `period` | `max_period`, `reduced`, `witness`, `bound`, `cover`, `saving`, `breaks`, `scanned` |

Exit code `1` means the library and the manifest disagree: the manifest is
stale, or the change was not intended — either way a human decides. Exit code
`3` means the manifest cannot be judged at all: not JSON, no `anchors` list, an
empty list, an unknown kind, a missing field, or a **float** where an exact
value belongs. A float is refused because it cannot express `1/3`, and an empty
manifest is refused because verifying nothing is not a green light.

### The other edge: the manifest against the acta it quotes

A manifest could be perfectly consistent with this package and still cite a
source that no longer says what it said. `source` closes that edge — the file,
its blob, its line count, and the tables the manifest quotes:

```json
"source": {
  "path": "acta_margen_secciones.py",
  "blob": "a6975e65920af7f5c66dc238e5938aeba10c11e4",
  "lines": 1606,
  "tables": [
    {"name": "FAMILIA_P54", "line": 240,
     "columns": ["grids", "universe_size", "margin", "critical_divisor",
                 "cost_cover", "cost_ceiling"]}
  ]
}
```

```
exact-lattice-quantiles verify-acta <manifest> --source <acta.py>
# 169 checks over 30 anchors: 169 reproduce, 0 divergent        RC 0
```

The acta is **parsed, not imported**: the reader walks the module's AST and
evaluates literals only — integers, strings, tuples, dicts, and `Fraction(...)`
/ `F(...)`. It never executes a line of a simulation just to read a table, and
it refuses a float for the same reason the manifest refuses one. A table
computed at run time therefore has to be written down before it can be quoted.

Three verdicts, and they are not the same verdict:

- a **blob that moved** while every quoted cell still agrees is `RC 0` with a
  printed `MOVIDO`: the numbers hold, the citation of the source is stale. A
  stale address is not a false claim.
- a **cell that disagrees** is `RC 1`, naming `table[period].column` — the same
  shape as a package divergence, so one reader covers both edges.
- a source that **cannot be measured** (no git, missing table, a float, no
  `source` block) is `RC 4 NO_MEDIDO`, which never degrades to `0`.

A cell no anchor declares is reported as `sin ancla` rather than passed: an
unmeasured cell is a blind spot with a name.

### The tables the boundary cannot see, and the one that looked like a defect

Twelve tables in that acta are keyed by period and written as literals, so the
reader can evaluate all twelve. Seven of them are quoted. Declaring the rest is
not mechanical, and one of them is a trap worth writing down.

Closed by declaring them: `TECHOS_FUERA_P55`, `U_FUERA_P55` and
`D_ATADO_FUERA_P55` carry the acta's out-of-sample claim — four highly composite
periods above 400 (`720, 840, 1260, 2520`) where the factorised ceiling must
still match the exact sweep. Anchoring them puts that claim on **both** edges at
once: this package recomputes the ceiling for those periods, and the acta's own
cells are re-read and compared. `RAD30_P55` is the six periods that share
`rad(T) = 30` and do not share a ceiling, which is the point of quoting them.

Partly anchored: `ESCALERA_47` and `ESCALERA_64`. The rows that do reproduce are
anchored, as a fourth kind, `ladder` — a row states a period, the design it is
measured on, the number of shared sections `k`, and the margin claimed. Unlike
`design`, the design is an **input** here: a ladder row is a claim about a
fixed design, and letting the design float would ask the question `design`
already answers. A ladder anchor also refuses a design that carries no block for
some divisor of the universe, because that row would report the margin of a
divisor nobody holds. Anchored: `k = 5, 7, 8` of `ESCALERA_47` and all three
rows of `ESCALERA_64`.

The rows that are **not** quoted look comparable — same period, same cost, a
margin per `k` — and they are not, for a reason that took measuring: the two
sides count the section space in different currencies. `ESCALERA_47`'s `k = 3`
row is not a margin at all — that design needs four sections even with no slack,
so at three there is no schedule whose margin could be quoted.

A block of `a` classes measures `2a / M` of the section circle — the period
currency, which is what every margin here is quoted in. The acta's P52 ladder
builds its cells as `a / (M // 2)` instead. The two agree whenever `M` is even
and differ on an odd grid by exactly `M / (M - 1)`. On the design that carries
`T = 180` the odd grid is 45, so its cells are stretched by `45 / 44`: on
`M = 45` every one of the five cells differs, while on 12, 18 and 20 all are
identical.

It is not that the acta uses one design and this package uses every design — at
cost 47 both pick `(12, 18, 20, 45)`, the same design, with the same cost. It is
that a "shared section `q`" in the stretched reading is not the same phase in
every grid, so sections that are shared there are not shared. Swapping only the
currency in the acta's own `escalera_global` — nothing else — turns its ladder
into this package's, row for row: `k = 4` goes `2/33 → 1/15`, `k = 6` goes
`4/33 → 2/15`, and `k = 3, 5, 7, 8` are unchanged, because those rows never
needed a section in the stretched corner. Its sibling `ESCALERA_64` — a design
with two odd grids, 15 and 45 — reproduces here at all three of its rows
(`1/6`, `3/10`, `1/3`), which is what makes the divergence a property of the
two rows and not of the table: the stretcher only bites where an optimum sits
in the stretched corner. The acta's own `asientos` computes
`2a / M`, its own `identidad_w` asserts it, and the P54/P55 ceiling law lives in
it — so the ladder is the one table in the acta that is in the other currency.
The consequence is that its `k = 4` value is not an upper bound: the same cost
buys `1/15 = 6.67 %` here against the `2/33 = 6.06 %` quoted there.

192 of the 200 even universes in `2..400` carry an odd grid — only a power of
two does not — so the two rows that move are recorded in the manifest's
`not_reproducible` list with that reason rather than anchored — together with
`PASOS_P55` (an instrument measurement), `CAIDAS_P55` (a census over 116
universes, not a quantity of one period) and `FUERA_P55` (a set per period,
which no column kind can hold yet).
`TestTheCurrencyHasAnIndependentReferee`
pins the currency against a referee that counts phases from `M / gcd(j, M)` and
never consults the seat map, so the two readings cannot quietly converge again.

## Citation

See `CITATION.cff`. Author ORCID: [0009-0001-9116-075X](https://orcid.org/0009-0001-9116-075X).

## License

MIT.

## Field notes

This repository belongs to an audit corpus whose measurements are published as
field notes. Two of them come from the same work:

- [Maths 2.0: why AI needs auditors, not just code generators](https://dev.to/alexisromtz7/maths-20-por-que-la-ia-necesita-auditores-no-solo-generadores-de-codigo-2f29)
- [Mide el jacobiano antes de optimizar: cuántas restricciones tiene de verdad tu objetivo](https://dev.to/alexisromtz7/mide-el-jacobiano-antes-de-optimizar-cuantas-restricciones-tiene-de-verdad-tu-objetivo-3oe1)
