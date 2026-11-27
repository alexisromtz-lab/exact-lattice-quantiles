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

This census prices the frontier with antichains of at most three parts, which is
cheap; the *check* is exhaustive over antichains of any size, which is not. So the
census publishes `909` verified periods rather than pretending to `1 000`, and the
check is the reason to trust the law outside them. Restricting the law to pairs is
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

And the third part is the third part in **either** universe, not only in the
reduced one. The full divisor lattice `D(T)` stores twelve periods at a fee of no
fee against the single `1080` of `U(T)` -- unbounded grids break the closed form
earlier and more often -- but where the closed form does hold, the pairs hold it:

```
$ exact-lattice-quantiles fee-census --universe divisors --max-period 400 --verify-universe 16
fee census: the fee that breaks the closed form, period by period
  horizon                 400             D(T) universe, blocks of at most 3 parts
  periods                 200             192 verified against every antichain (|U| <= 16)
  mismatches              0               periods where no block of that size attains the frontier
  chains                  8               no antichain at all: no fee breaks these
  at zero                 12              already broken at no fee (60, 90, 120, 150, 180, 210, 240)
  at one                  9               break first at a fee of 1
  widest frontier         100             period 398 on (2, 199)
  witness sizes           2: 169, 3: 15   on the verified periods
  pairs miss              21              a triple beats every pair
  pair over               0               of those, above a fee of one: the law refuted if not 0
  pair law 21 periods put a triple above every pair and all of them sit at a fee of at most one, so the third part buys exactly one thing: the fee of one
  verdict  the law holds on all 192 periods of D(T) with at most 16 grids, witnesses of size [2, 3]; 8 periods are chains where no fee breaks the formula, 12 break with no fee at all [60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390], and 9 break first at a fee of 1 on 2 shapes (6, 10, 15), (6, 14, 21); the widest frontier is 100 at period 398 on (2, 199)
```

`21 = 12 + 9` the same way `52 = 51 + 1` did, on a different universe and with
different witnesses: the fee of one still needs `(6, 10, 15)` or `(6, 14, 21)`
and nothing above it ever needs a triple. What changes between universes is not
the pair law but where the frontier reaches zero -- and that is published as data,
not assumed by the law.

#### The fourth part is free, and the free rung is Giuga's

The frontier is a minimum over blocks, so a bigger block can only lower it, and
until v0.19.0 this package stopped the claim at three parts. That stop was wrong,
and it was wrong in the direction that matters: the census verified every period
whose universe fits in **14** grids, and the first block that reaches a fee of
zero needs **15**. A verification horizon one grid short of the first witness is
a green in the dangerous direction, so v0.20.0 publishes the horizon and the
witness side by side:

```
$ exact-lattice-quantiles arity-ladder --max-period 1000
arity ladder: how many parts the fee frontier needs
  horizon               1000              divisors universe, rungs (2, 3, 4, 5)
    2->3                52                first 30 on (6, 10, 15) (2 -> 1, 7 grids)
    3->4                7                 first 462 on (42, 66, 154, 231) (1 -> 0, 15 grids)
    4->5                0                 never earns
  first free block      60                (12, 20, 30), already free with 3 part(s)
  periods at zero       38                (60, 90, 120, 150, 180, 210, 240)
  verified              442               against every antichain (|U| <= 15)
  mismatches            0                 the ladder against the exhaustive frontier
  verdict  the ladder earns 2->3: 52 3->4: 7 4->5: 0 period(s); a free block first appears at 60 on (12, 20, 30), which 3 part(s) already reach, so the last step that earns is 3->4, at 7 period(s), first 462 on (42, 66, 154, 231); no block above 4 parts earns over one at 4 on any period up to 1000, and the ladder matches the exhaustive frontier on all 442 periods whose universe is at most 15 grids
```

Three statements come out of it, and the first one corrects the package:

- **The free rung opens with three parts, at `60`.** Not with four. `(12, 20, 30)`
  costs nothing at all, and it is the cofactor block of the primes `2, 3, 5`.
- **The fourth part is free where the third stalls at one.** `462 = 2 * 3 * 7 * 11`
  has `(42, 66, 154, 231)` free while every triple there costs `1`, and the seven
  periods up to `1000` where that happens are `462, 546, 714, 798, 858, 924, 966`.
- **Nothing above four earns, because zero is the floor.** A fee is a count: once
  a rung reaches `0`, no higher rung can improve on it. That is why the ladder has
  two earning steps and no third.
- **Over the carriers the freeness arrives by scaling, at `1080`.** Not by a
  fourth part and not inside the popular horizon of `1000`: `(54, 90, 135)` is
  nine times the Giuga triple and costs nothing while the triple pays `1`. A fee
  cannot see that distance, which is why the instrument that measures it is the
  margin and not the fee.
- **The crossing survives a second reader that shares no bound.** Grouping blocks
  by their lcm instead of walking antichains -- same maximum, no objective bound
  to be wrong -- reads the same `1080` and splits on no period in either universe,
  and the plausible shortcut of blocking only the maximal grids would have read
  `-20` on `U(60)` where the truth is `0`.

The reason a block can be free at all is an inversion worth stating plainly: the
block is cheaper than buying each part, *and* its parts are dearer than the one
grid that covers them. With `L = lcm(A)` and `delta = sum - L`,

```
free  <=>  delta > k  <=>  2B < 0  <=>  sum over the complements  1/(L/a) > 1 + k/L
```

so a free block is a set of divisors whose **reciprocals pass one by the share of
its odd parts** -- and no less. Reading the last link as `> 1` was the package's
own published error until v0.21.0, and its own README witness falsifies it:
`(6, 10, 15)` has complements `(5, 3, 2)`, so `1/5 + 1/3 + 1/2 = 31/30 > 1`, and
it still pays `1`, because one odd part moves the threshold to
`1 + 1/30 = 31/30` -- the block sits *on* the line, not past it. The smallest
set of primes that does is `1/2 + 1/3 + 1/5 = 31/30`, which is why `60` is the
first period: there the cofactors `(30, 20, 12)` are all even, `k = 0`, and a
defect of `2` already frees them -- whereas at `30` the cofactors `(15, 10, 6)`
contain an odd part, the defect is `1`, and the same prime set costs `1`. One
exponent of `2` is the whole difference.

And **no pair can ever be free**, which is what puts `2` at the bottom of the
ladder instead of `1`: for incomparable `a < b`, `lcm(a, b) >= 2b > a + b`, so the
defect is at most `-1` and `B >= 1` in either parity -- a fee of at least two,
measured over the `8 062` incomparable pairs of both universes up to `600`.

The far end of the same statement is Giuga's equation, which the package met in
v0.10.0 from the odd band and meets again here from the free rung. For the
cofactor block of a squarefree `L`, `delta = L * (sum 1/p - 1)`, so

- `delta = 1` is `sum 1/p - 1/L = 1`, **Giuga's equation**: `30`, `858`, `1722`.
  A block one defect short of free, and it pays `1`, not `0`.
- `delta = -1` is `sum 1/p + 1/L = 1`, the **primary pseudoperfect** periods,
  whose smallest members with two primes are `6` and `42`.
- `delta >= 2` is past Giuga, and it is free: `{2, 3, 7, 11}` gives `462` with a
  defect of `31`. `{2, 3, 7}` alone gives `41/42 < 1`, which is exactly why the
  quadruple is needed *there* and not at `60`.

`cofactor_block((2, 3, 7, 41))` and `cofactor_block((2, 3, 7, 11))` are the two
sides of that boundary -- `1722` pays one, `462` pays nothing -- and `block_free`
is the predicate that separates them.

The ladder is not the same over the other universe, and the census does not
average the two:

```
$ exact-lattice-quantiles arity-ladder --max-period 1000 --universe reduced
arity ladder: how many parts the fee frontier needs
  horizon               1000              reduced universe, rungs (2, 3, 4, 5)
    2->3                25                first 60 on (6, 10, 15) (2 -> 1, 8 grids)
    3->4                0                 never earns
    4->5                0                 never earns
  first free block      none in horizon   every period here keeps a fee of at least one
  periods at zero       0                 none in this horizon
  verified              468               against every antichain (|U| <= 15)
  mismatches            0                 the ladder against the exhaustive frontier
  verdict  NO_PERIOD REACHES ZERO up to 1000 over the reduced universe: the ladder earns 2->3: 25 3->4: 0 4->5: 0 period(s) and no rung is free, so the floor of the fee is one here rather than zero
```

Over the carriers `U(T)` no period up to `1000` reaches zero: the floor of the fee
is **one**, and the census says so with a sentence rather than with a zero, which
would be a frontier of no fee. Over the divisors `D(T)` it is zero from `60` on.
Two universes, two floors, both published -- the free rung is a property of the
lattice you price, not of the fee. And `1000` is not where the carriers stop: the
floor of one is a statement about this horizon, and the next section measures how
far past it the witness sits.

Exit codes are the ladder's own: `RC 0` when the ladder *is* the frontier on every
period it could check exhaustively, `RC 1` when it is not (the package refuting
itself), and `RC 4` when it checked nowhere -- because a ladder nobody verified is
not a ladder measured.

#### The fee cannot see the distance: the carriers cross at `1080`

The transcript above reads a *fee*, and a fee is a quotient. `block_fee` divides a
block's margin by twice one less part, so on a triple the single fee `1` is
charged for margin `0` -- free at the next multiple -- and for margin `-2` alike.
A horizon chosen on fees therefore cannot tell a period sitting on the free line
from one below it, and that is how the carriers were published with a floor of
one when their first free block was eighty periods away. The margin itself is an
integer with no such blind spot: `free_margin(A) = 2 * block_advantage(A)`, which
is `delta - k` when the lcm is even and `delta - (n - 1)` when it is odd, so a
block is free exactly when its margin is positive.

```
$ exact-lattice-quantiles free-margin
free margin: the distance to a free block, in integers
  horizon                 1200            carriers universe, 590 periods measured
  first crossing          1080            (54, 90, 135), margin 8, 23 grids
  touching the line       19              margin exactly 0, free at the next multiple: 60, 120, 180, 240, 300, 360, 420, 480
  strictly below          570             every measured period with a negative margin
  widest margin           8               at 1080 on (54, 90, 135)
  verdict  the carriers cross the free line at 1080 on (54, 90, 135) with margin 8, after 19 periods touching it at margin 0
```

Nineteen periods *touch* the line, every one of them a multiple of `60`, and every
one of them on the same block `(6, 10, 15)` -- the Giuga triple, whose defect
equals its single odd part exactly. The first period to cross is `1080 = 18 * 60`,
on `(54, 90, 135) = 9 * (6, 10, 15)`: **the carriers reach free by scaling, not
by adding a part.** The block is the same three parts multiplied by nine, and it is
in `U(1080)` because `1080` is divisible by all three of them.

The mechanism is one line of algebra. The defect and the lcm both scale,
`delta(kA) = k * delta(A)` and `lcm(kA) = k * lcm(A)`, while the count of odd
parts does not: `kA` is all even for even `k` and keeps `A`'s odd parts for odd
`k`. So `margin(kA) = k * delta(A) - (k odd ? #{odd parts of A} : 0)`, which for
the Giuga triple is `k` for even `k` and `k - 1` for odd `k` -- free for every
`k >= 2`. Which `k` a carrier universe admits is then the only question, and up
to `4000` the answer is exactly `9`, `14`, `18` and `22`: the triple itself costs
`1`, and every scaling that lands in the lattice lands free.

Over the divisors the same instrument reads `30` and `858` -- the two smallest
Giuga numbers inside the horizon, each on its own cofactor block, each paying
exactly `1` -- and crosses from `60`:

```
$ exact-lattice-quantiles free-margin --universe divisors
free margin: the distance to a free block, in integers
  horizon                 1200            divisors universe, 590 periods measured
  first crossing          60              (12, 20, 30), margin 2, 11 grids
  touching the line       2               margin exactly 0, free at the next multiple: 30, 858
  strictly below          542             every measured period with a negative margin
  widest margin           184             at 1050 on (150, 210, 350, 525)
  verdict  the divisors cross the free line at 60 on (12, 20, 30) with margin 2, after 2 periods touching it at margin 0
```

Two universes again, and this time the numbers say what the fees only hinted at:
the divisor lattice *leaves* the line at `60`, while the carrier lattice sits on it
nineteen times and only `1080` gets past inside this horizon.

`free-margin` exits `RC 0` when its own witness reproduces -- the block is in the
universe, `block_free` agrees, and the margin is the one claimed -- `RC 1` when it
does not, and `RC 4` when the horizon had no pair of incomparable grids at all, so
nothing was measured. There is no `--exhaustive-until` here on purpose: the search
walks every antichain of the universe itself, so its horizon is the horizon and
not a calibration of one.

#### Read twice: the same crossing by lcm, without a bound

Every number above comes out of one instrument, and an instrument's number is a
claim about the instrument. The search is a depth-first walk whose cost is paid
by a bound -- `2 * (halves + rest) - span <= best` -- and that bound is where a
bug would live: too small and it silently drops the crossing, and the census
reports a floor the universe does not have. So the question is asked again by a
reader that shares nothing with it. A margin is fixed by the block
(`2 * sum(a // 2) - lcm`), so the blocks partition by their lcm, and every member
of a block of lcm `L` must divide `L`: the candidates are known before any choice
is made, each class is small enough to enumerate whole, and **there is no
objective bound at all** -- the only cut left is structural. A wrong bound hides
one block; it cannot hide a class.

```
$ exact-lattice-quantiles cross-read
cross-read: the free margin by search and by lcm class
  horizon                 1200            carriers universe, 590 periods read twice
  by search               pruned DFS      walk over antichains, cut by an objective bound
  by lcm                  class by class  every block of lcm L, enumerated whole, no bound to be wrong
  disagreements           0               none
  witnesses               0               blocks that fail to reproduce their own margin
  first crossing by lcm   1080            (54, 90, 135), margin 8, 3 parts
  on the line by lcm      19              60, 120, 180, 240, 300, 360, 420, 480
  crossed by lcm          1               1080 (8)
  verdict  the two readers agree on all 590 periods to 1200 over the carriers, and every witness reproduces its own margin
```

The classes are built from the closure of `lcm` over the universe and not from
pairs, because three grids can reach an lcm no pair of them reaches:
`lcm(2, 3, 5) = 30` is the lcm of no pair of `{2, 3, 5}`, so a horizon generated
pair by pair would skip the Giuga triple while agreeing on everything around it.

Agreement is only half the claim; the other half is that the second reader
*would* have caught the bug. Take the plausible simplification -- only maximal
grids are worth blocking -- and it reads `-20` on `U(60)`, where the truth is `0`
on `(6, 10, 15)`: the line at `60` lives strictly below the maxima, and a reader
that looks only at the top of the lattice cannot see it. A test states exactly
that, so the cross-read is not a page of green: it is an instrument whose
disagreement is a measured possibility.

Read far enough and the crossings are more than one. To `3000` there are three,
and every one of them is the Giuga triple scaled:

| period | block | `k` | margin |
|---|---|---|---|
| `1080` | `(54, 90, 135)` | 9 | 8 (`k - 1`, odd `k`) |
| `2160` | `(108, 180, 270)` | 18 | 18 (`k`, even `k`) |
| `2520` | `(84, 140, 210)` | 14 | 14 (`k`, even `k`) |

Over the divisors the second reader puts a name on the line. Extended to `2400`
the divisor line is touched at `30`, `858` and `1722` and **nowhere else**, on
`(6, 10, 15)`, `(66, 78, 286, 429)` and `(42, 246, 574, 861)` -- the cofactor
blocks of `2*3*5`, `2*3*11*13` and `2*3*7*41`, which are the three smallest Giuga
numbers:

```
$ exact-lattice-quantiles cross-read --universe divisors --max-period 2400
cross-read: the free margin by search and by lcm class
  horizon                 2400            divisors universe, 1189 periods read twice
  by search               pruned DFS      walk over antichains, cut by an objective bound
  by lcm                  class by class  every block of lcm L, enumerated whole, no bound to be wrong
  disagreements           0               none
  witnesses               0               blocks that fail to reproduce their own margin
  first crossing by lcm   60              (12, 20, 30), margin 2, 3 parts
  on the line by lcm      3               30, 858, 1722
  crossed by lcm          99              60 (2), 90 (2), 120 (4), 150 (4), 180 (6), 210 (36)
  verdict  the two readers agree on all 1189 periods to 2400 over the divisors, and every witness reproduces its own margin
```

That row is a prediction that was written before the sweep: if the zeros of the
divisor line are Giuga's equation `sum 1/p - 1/L = 1`, the third Giuga number is
on the line and nothing else is. It is, and nothing else does.

`cross-read` exits `RC 0` when the two readers agree and every witness reproduces
its own margin, `RC 1` when they split or a witness contradicts itself -- the
horizon then is the period, not the census -- and `RC 4` when no period had two
incomparable grids to read at all.

#### The universe is a filter over divisors, and the scales are a law

`universe(T)` builds every even grid of the period and unions their medians: a
*forwards* definition, and forwards the answer costs `T / 2` seats per grid. The
same set can be read *backwards*, and backwards it is a property of the divisor
alone -- `a` is in `U(T)` exactly when some even grid `q * a` dividing `T`
carries `a`. That is one filter over the divisors of `T`, no grid built, and
`carrier-law` reads both ways and compares, period by period:

```
$ exact-lattice-quantiles carrier-law
carrier-law: U(T) read by the divisor, and the divisors no grid carries
  periods read twice      600             every even period to 1200
  splits                  0               none, the two routes agree
  carriers measured       922             divisors to 1000 that some even grid carries, by window and not by limit
  largest cofactor        4               no carried divisor to this horizon needs more
  no grid carries         77              30, 42, 60, 66, 78, 102, 114, 120
  cross-check             20000           grids read to this cofactor, against the window
    window short at       nowhere         the theorem holds
    limit short at        nowhere         the grids reached every one
  verdict  the divisors and the grids name the same universe on all 600 even periods to 1200
  verdict  the window and the grids name the same 77 uncarried divisors to 1000, and the law needs no limit to do it
```

Two facts, and the second one is the useful one. The routes agree on all **2000**
even periods to `4000`, so the law is the same set and costs `tau(T / a)` mass
reads per divisor instead of `T / 2` seats per grid. And a divisor that is no
even grid's median is in **no** universe at all, for any period: up to `1000`
there are `77` such divisors, and every one of them is even. The three periods
where the *divisor* universe touches the free line -- `30`, `858`, `1722` -- are
all in that list, as is `60`, where the carrier universe touches. The line is
where a divisor is not carried.

#### A grid limit was never needed: it was a window

Reading the image backwards said *which* divisors are uncarried but still swept
the grids forwards to find out, and a sweep needs a horizon -- the v0.23.0
report literally named `20 * max_divisor` and left the rest unmeasured. The
horizon can be replaced by a theorem, because the median of `M = q * a` is a
running sum of the mass and two sums bracket it:

- *everything below `a` falls short of the middle*, so `cum(a-) <= M / 2`; every
  divisor of `a` but `a` itself divides `M` and lies below `a`, so
  `cum(a-) >= a - mass(a)`, hence `q >= 2 (a - mass(a)) / a`;
- *the middle is reached at or before `a`*, so `M / 2 <= cum(a)`, and `cum(a)`
  sums a subset of the integers up to `a`, hence `q <= 2 Phi(a) / a`.

Neither bound mentions `T`, a grid or a cofactor beyond `q`, so whatever carries
`a` carries it inside **one window read off `a` alone** -- and the window is
tight: at `a = 2` it is `[1, 2]` and the ceiling is attained, so no argument
reading only `a` lowers it. The floor is a ceiling and not a strict successor
because of the tie: when `4` divides `M` and the mass lands exactly on the
middle, the divisor straddling it carries too, which is where `q = 1` survives
for every power of two. Read strictly, that tie drops `4`, `8`, `16`, ... from
their own image -- how the window was first written, and wrong.

What the theorem buys is that the second list stops depending on a horizon:

```
$ exact-lattice-quantiles carrier-law --max-divisor 2000
  carriers measured       1838            divisors to 2000 that some even grid carries, by window and not by limit
  largest cofactor        5               no carried divisor to this horizon needs more
```

**And that `5` is a horizon, not a bound.** The least witness -- the smallest
cofactor that carries a divisor -- climbs a staircase, and every step is a
divisor that is the median of exactly *one* grid:

| `a` | 2 | 6 | 84 | 90 | 1848 | 4620 |
|---|---|---|---|---|---|---|
| least witness | 1 | 2 | 3 | 4 | 5 | 6 |

`carrier_cofactors(4620)` is `(6,)` inside a window of `[2, 2808]`: one carrier
grid, `27720`, and that grid is outside the table a `horizon = 20000` read ever
builds. So the transcript above answers *how far the read went*, never *how far
the object goes* -- which is why the line says "to this horizon". No even divisor
below `4620` has witness `6`, and none below `300000` has `7`; the staircase is
measured, and a measurement is exactly what a bound is not. The proof of the
window does not depend on it: the window is read off `a` alone and holds whether
the witness is `2` or `6`.

Two readers now disagree in *two distinct ways*, and `carrier-law` keeps them
apart, because only one of them is the theorem failing. A divisor the grids carry
and the window calls absent is the window missing a cofactor: `RC 1`, and no
wider sweep repairs it. A divisor the window carries and the grids miss is the
grid limit being short: the measurement had not looked far enough, and the law
was right. `tests/test_fee_law.py::TheCarrierWindow` mechanises both directions --
it cripples the window by one cofactor and asserts the census says *the theorem
failing*, not *the limit short*.

And the complement is **not** characterized by `a / phi(a)`. `210` has ratio
`4.375` and `630` carries it at cofactor `3`; `30` sits at `3.75` and no grid
carries it at all. So the uncarried divisors are a measurement and not a closed
form, and the package says so instead of guessing one.

That is what turns the scales of the Giuga triple from a sweep into a law. The
margin is not the question, because it is closed: `margin(kA) = k - [k odd]`,
positive for every `k >= 2`. What decides a crossing is whether the three scaled
parts can land in the same `U(T)` -- a membership question, answered by the
carrier grids of `6k`, `10k` and `15k` and one lcm, with no period swept:

```
$ exact-lattice-quantiles scalar-law --max-period 4001
scalar-law: every scale of (6, 10, 15) that clears the free line
  horizon                 4001            2000 even periods read twice
  crossing periods        5               1080 (k=9, margin 8), 2160 (k=18, margin 18), 2520 (k=14, margin 14), 3240 (k=27, margin 26), 3960 (k=22, margin 22)
  reader splits           0               none, both readers agree
    measured 1080         8               (54, 90, 135)
    measured 2160         18              (108, 180, 270)
    measured 2520         14              (84, 140, 210)
    measured 3240         26              (162, 270, 405)
    measured 3960         22              (132, 220, 330)
  verdict  the law and the census name the same 5 crossing periods to 4001, out of 2000 even periods read twice
```

`3240` is the row that matters. Every census published before v0.23.0 stopped at
`3000` and so saw three crossings; the law reads `k = 27` off the grids of
`162`, `270` and `405` and names a period thirty past that horizon -- and the
universe of `3240` then confirms it with both readers, at margin `26`. A law
that only reproduced the rows already measured would be a summary; one that
names a row no sweep had reached is a forecast. `scalar-law` exits `RC 1` the
moment the two readers split or the census finds a crossing the law did not
name, so a sixth scaled crossing cannot arrive quietly.

The other direction is measured too, past the horizon the two readers were built
for: `cross-read --universe divisors --max-period 10000` reads **4987** periods
twice in three minutes fifty and finds `0` disagreements, the line touched only
at `30`, `858` and `1722`, and `418` crossed periods with `60` at margin `2` on
`(12, 20, 30)` as the first. The fourth Giuga number is `66198`, so nothing new
was owed on the line -- what was owed was the certificate, and it is now a
transcript instead of a bound.

#### Three routes to one universe, and the tie is what tells them apart

Everything above rests on `U(T)`, and until v0.25.0 the package read it twice in
a way that was not independence. `universe` walks the grids and takes the
carriers; `carrier_universe_by_law` walks the divisors and takes the carriers;
both end at `carriers`, and `carriers` is one loop that accumulates `mass`
(`euler_phi`) over `divisors`. A bias in that loop would have travelled with
both readings and come out looking like agreement: the double reading was
measuring the same route twice, which is exactly the shared blind spot the
package had been naming without closing.

So there is a third route, built to share nothing. Class `j` of grid `M` meets
the divisor `M // gcd(j, M)` -- the phase count the width contract already
trusts -- so reading the **median of the class values** names the carrier
without ever touching `mass`, `divisors`, `seats` or a factorization; the grids
themselves come from trial division instead of from the prime walk. The only
thing this route shares with the other two is the definition of a carrier.

```
$ exact-lattice-quantiles median-oracle --max-period 1200
median-oracle: U(T) read from the classes, a third time
  periods read            600             every even period to 1200
  grids read              600             even grids whose carrier is read both ways
  widest universe         23              carriers in the largest U(T) read
  grid splits             0               trial division and the factorization agree
  carrier splits          0               the median of the class values is the mass carrier
  universe splits         0               the same set, period by period
  verdict  the median of the class values names every carrier the mass names, over 600 even grids of 600 periods to 1200, so U(T) is a property of the lattice and not of the route that reads it
```

Two things fall out, and the second is the reason to have built it.

- **`mass` is a theorem about the class map, not an assumption of it.** A
  divisor `d` appears among the class values of `M` exactly `mass(d) / 2` times,
  which is why the accumulated mass seats a block of `d` at the middle of the
  grid. The class route derives that identity on every divisor of every even
  grid to `360` and then never consults it.
- **The tie branch is load-bearing, and the universe hides the damage.** Read
  the median as a single middle value -- the natural wrong reading, since a
  median usually is one -- and carrier `2` of grid `4`, carrier `4` of grid `8`,
  and every `4 | M` whose two middle classes differ, drop out. No universe below
  `36` changes, because other grids of those periods still carry the same
  divisor: **the union absorbs the mistake**. To `1200` that mutant parts from
  the mass route on `29` grids, the first at `4`, but on only `33` universes,
  the first at `36` -- so the universe is the claim and the grid is the
  instrument, and a census that compared only `U(T)` would have reported
  agreement on a wrong reading. `tests/test_median_oracle.py` patches the median
  that way and requires the split to be reported, grid by grid, with `RC 1`.

Measured, `0` splits of any kind to `1200` (600 grids over 600 periods, widest
universe `23`), to `4001` (2000 and 2000, widest `33`, `1.0 s`) and to `6000`
(3000 and 3000, widest `42`, `2.3 s`) -- every even grid below each horizon, with
the grids both ways and the universes period by period. The verb exits `RC 1` on
the first disagreement and `RC 3` on a horizon below `2`; `RC 4` is not
reachable here, because a grid exists for every even period and so the census
can never come back empty.

#### A necessary condition that is not a decision procedure, and why that is arithmetic

The carrier cofactor window is a theorem: a divisor `a` is carried by the grid
`q * a` only if `q` lies in `[⌈2(a - mass(a))/a⌉, ⌊2Φ(a)/a⌋]`, both ends read off
`a` alone. It is tempting to read it backwards -- *inside the window, therefore
carried* -- and the backwards reading is false at `a = 30`.

What is worth noticing first is that emptiness was never going to be the
evidence. `Σ_{d | a} φ(d) = a` puts the floor at most `2` (because
`a - mass(a) < a`) and the ceiling at least `2` (because `Φ(a)` is a sum that
already contains `Σ_{d | a} φ(d) = a`), so **`q = 2` is in every window the
package ever builds**, for every `a ≥ 2`. A census reporting `0` empty windows
is reporting arithmetic, not a near miss.

The informative count is the **gap**: divisors whose window is nonempty and none
of whose cofactors inside it carries them. `carrier_window_census` measures it
against `carriers`, with no grid limit anywhere, and the first member is `30` --
`carriers(30) = (15,)`, so `30` is a carrier of nothing, while `carrier_cofactors(30) = ()`
sits inside `[2, 18]`.

```
$ exact-lattice-quantiles carrier-law --max-period 600 --max-divisor 600
carrier-law: U(T) read by the divisor, and the divisors no grid carries
  periods read twice      300             every even period to 600
  splits                  0               none, the two routes agree
  carriers measured       554             divisors to 600 that some even grid carries, by window and not by limit
  largest cofactor        4               no carried divisor to this horizon needs more
  no grid carries         45              30, 42, 60, 66, 78, 102, 114, 120
  window empty            never           sum(phi(d) for d | a) = a keeps q = 2 inside every window
  in window, uncarried    45              30, 42, 60, 66, 78, 102, 114, 120
    gap all sixes         yes             every one a multiple of six, no counterexample to 600
  widest window           599             [1, 365] is 364 cofactors wide
  ...
  verdict  the window is necessary and not sufficient: 45 of the 599 divisors
           to 600 sit in a nonempty window that carries them nowhere, first 30,
           and every one of them is a multiple of six, so no bound read off a
           alone decides membership
```

Measured to `600`: `45` gaps of `599` divisors, all of them multiples of six,
`0` empty windows, widest window `[1, 365]` at `599`. To `3000`: `245` of `2999`
(`8.2` per cent), first six gaps `30, 42, 60, 66, 78, 102`, again every gap a
multiple of six and no window empty, widest `[1, 1824]` at `2999`, in `18 s`.

Two consequences, and they point in opposite directions.

* The window cannot replace the search, so `carrier_cofactors` keeps calling
  `carriers`: a period census that trusted the window alone would hold `245`
  divisors that no period holds.
* The gap is not the evens alone but the **multiples of six**, at every horizon
  measured, with no counterexample. That is a census claim and not a theorem, and
  it is stated with its horizon because it has one: `tests/test_carrier_gap.py`
  pins it to `600` and would report the first counterexample rather than defend
  the pattern.

The `largest cofactor` line is a horizon fact too, and its own staircase keeps
climbing: the least cofactor a carried divisor needs is `1, 2, 3, 4, 5, 6` at
`a = 2, 6, 84, 90, 1848, 4620`, so a census that read `4` to `600` as a ceiling
would be wrong at `1848` and the one that read `5` to `3000` would be wrong at
`4620`. No cofactor bound is assumed anywhere: `_least_carrier_cofactor` walks
the window until it finds the hit or reaches the end of it.

#### The fee of a block is not a search: it is half its defect

`block_advantage` is affine in the fee with slope `n - 1`, so inverting it is one
integer division -- and the numerator is not a coincidence. Write `delta` for
the **defect** of a block, `sum - lcm`, the distance between buying every part
and buying the one grid that covers them all, and `k` for the number of its odd
parts. Then

- if the lcm is even, `floor(lcm / 2) = lcm / 2` and every odd part contributes
  `(p - 1) / 2`, so twice the advantage at a zero fee is `k - delta`;
- if the lcm is odd, every part is odd, so it is `n - 1 - delta`.

Either way `2 * B` is an integer and

```
fee(block) = 0                          if 2 * B < 0
fee(block) = 2 * B // (2 * (n - 1)) + 1  otherwise
```

No antichain, no walk, no search: the fee of a block is **the defect, halved and
divided by one less part**. Two shapes follow, and both are the ones the census
above kept finding:

```
$ exact-lattice-quantiles block-fee 6 10 15
block fee: (6, 10, 15)
  lcm                     30              one grid covering every part
  sum                     31              one grid per part
  defect                  1               sum - lcm: zero is a perfect block
  odd parts               1               k - defect, the lcm is even
  twice B                 0               twice the advantage at a zero fee
  fee                     1               floor(twice B / 4) + 1
  by advantage            1               the same fee, inverted the slow way
  verdict  breaks the closed form from a fee of 1 up

$ exact-lattice-quantiles block-fee 6 21 28
block fee: (6, 21, 28)
  lcm                     84              one grid covering every part
  sum                     55              one grid per part
  defect                  -29             sum - lcm: zero is a perfect block
  odd parts               1               k - defect, the lcm is even
  twice B                 30              twice the advantage at a zero fee
  fee                     8               floor(twice B / 4) + 1
  by advantage            8               the same fee, inverted the slow way
  verdict  breaks the closed form from a fee of 8 up
```

The two lines matter together: a block **one short of perfect** is the cheapest
thing that can break the formula, and the further its lcm runs above its sum the
more fee it needs before it is worth folding. That is why the fee of one needs a
triple and why there are only two of them below 2000 -- for a triple with an even
lcm the fee is one exactly when the defect sits within three of the odd-part
count, a window two blocks wide in the whole horizon.

The closed form is not a replacement for the instrument, it is a second reading
of it, and the package keeps both: `block_fee` inverts the advantage the slow
way, `block_fee_closed` divides once, and the two are compared on **every** block
of two and three parts in both universes up to `T = 400` -- `42 287` blocks, `0`
mismatches, which is the shape of a proof rather than a census. `block-fee`
prints both rows side by side so a disagreement would be visible rather than
argued about, and a block of one grid is refused with `RC 3` instead of being
priced at zero.

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
python -m unittest discover -s tests -t .     # 479 tests
exact-lattice-quantiles --self-test           # 127 checks, RC 0
exact-lattice-quantiles median-oracle --max-period 6000    # 0 splits, RC 0
exact-lattice-quantiles carrier-law --max-divisor 3000     # 245 gaps of 2999, RC 0
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

Version 0.26.0 is archived in Software Heritage, and the identifiers below are
the archive's answer and not a string this repository composed:

    swh:1:rev:2eddead349935f9f478959cabf35e2801e4e7d89
    swh:1:dir:dacb487a3f69a825d38395541cf12b4ffb29a74f

The revision is the merge commit of `#34`; the directory is what the archive
resolved from it, and both the revision and the directory API answer `200`. They
are quoted only after `/api/1/revision/<sha>/` returned them: an identifier that
resolves is evidence, and one that does not is decoration. No DOI is published
here, because a DOI has not been minted -- the deposit needs a token this
repository does not carry, and a written DOI that resolves to nothing would be
the same failure this package exists to measure.

## License

MIT.

## Field notes

This repository belongs to an audit corpus whose measurements are published as
field notes. Two of them come from the same work:

- [Maths 2.0: why AI needs auditors, not just code generators](https://dev.to/alexisromtz7/maths-20-por-que-la-ia-necesita-auditores-no-solo-generadores-de-codigo-2f29)
- [Mide el jacobiano antes de optimizar: cuántas restricciones tiene de verdad tu objetivo](https://dev.to/alexisromtz7/mide-el-jacobiano-antes-de-optimizar-cuantas-restricciones-tiene-de-verdad-tu-objetivo-3oe1)
