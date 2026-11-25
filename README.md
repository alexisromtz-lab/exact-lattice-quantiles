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
| Which margin costs `c` points? | the whole trade-off curve | `ladder` |
| How many shared sections hold margin `v`? | set cover over the divisor universe | `minimum_probes` |
| What is the best margin `k` sections can hold? | maximin over probe positions | `best_probes` |
| Is the margin I asked for even possible? | a certificate, not an error message | `InfeasibleAboveInstanceCeiling` |
| How much common drift does the margin survive? | a curve, not a constant | `drift_tolerance` |
| How much common drift still covers every divisor? | the exact coverage radius | `drift_radius` |

## Install

```
pip install git+https://github.com/alexisromtz-lab/exact-lattice-quantiles
```

Python 3.11+. No dependencies, at runtime or for tests.

## Quick start

```python
from fractions import Fraction
from exact_lattice_quantiles import instance_ceiling, solve_lattice_schedule

instance_ceiling(180).margin          # Fraction(1, 3)  — tied by divisor 6
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
pins the day the ranking dies. What neither says is *why* the sign reading
holds. The honest statement is the one in the code: `antichain_saving_ceiling`
is a search, and a proof would replace it with a one-line inequality.

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
full table is `2^23` masks, and `2520` is `2^32`.

So the library never prices with the formula: `price_exact` is a subset DP over
the coverage masks, and `price_form_census` publishes how often the shortcut
agrees. Restricting the universe to `U(T)` is what makes the shortcut look
right, and that restriction is a modelling decision, not a measurement --
`price_form_census(period, reduced=False)` is the row that shows it.

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

## What the three exit codes mean

The distinction is the point: an infeasible margin cannot be bought, an
over-budget one can.

| RC | exception | meaning |
| --- | --- | --- |
| 0 | — | a schedule exists and is returned |
| 1 | `InfeasibleAboveInstanceCeiling` | above `nu(T)`: no design at **any** cost |
| 2 | `InfeasibleWithinBudget` | exists, but `cheapest > max_cost` |
| 3 | `ValueError` (CLI) | malformed request — bad fraction, bad flag |

`PROBES_INSUFFICIENT` is a verdict on the returned `Schedule`, not an
exception: the design was found, the sections were not enough, and the object
says by how much (`probes_needed`, `witness_margin`).

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
  span geometry in **1712 of 1712** designs — 200 even periods, budgets
  `{16, 24, 47, 64, 96}`, `k` in `{2, 4, 6, 8}`;
- for **78 of 240** designs at `k = 4` and `k = 6` the first divisor break is
  strictly smaller than the radius, and the design pays nothing at it.

[tests/test_drift_law.py](tests/test_drift_law.py) pins the per-divisor half —
the second check fails if `_linear_regime` is replaced by "the line holds
everywhere" — and [tests/test_design_linearity.py](tests/test_design_linearity.py)
pins the design half.

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
python -m unittest discover -s tests -t .     # 184 tests
exact-lattice-quantiles --self-test           # 38 checks, RC 0
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
# 84 checks over 24 anchors: 84 reproduce, 0 divergent        RC 0
```

An anchor is one of three kinds:

| kind | required | optional |
| --- | --- | --- |
| `ceiling` | `period`, `margin` | `critical_divisor`, `grids`, `universe_size`, `cost_cover`, `cost_ceiling` |
| `design` | `period`, `max_cost`, `probes`, `margin` | `cost`, `design`, `probes_needed`, `robust_margin`, `verdict` |
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
# 163 checks over 24 anchors: 163 reproduce, 0 divergent        RC 0
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
reader can evaluate all twelve. Four of them were quoted. Declaring the rest is
not mechanical, and one of them is a trap worth writing down.

Closed by declaring them: `TECHOS_FUERA_P55`, `U_FUERA_P55` and
`D_ATADO_FUERA_P55` carry the acta's out-of-sample claim — four highly composite
periods above 400 (`720, 840, 1260, 2520`) where the factorised ceiling must
still match the exact sweep. Anchoring them puts that claim on **both** edges at
once: this package recomputes the ceiling for those periods, and the acta's own
cells are re-read and compared. `RAD30_P55` is the six periods that share
`rad(T) = 30` and do not share a ceiling, which is the point of quoting them.

Deliberately **not** quoted: `ESCALERA_47` and `ESCALERA_64`. They look
comparable — same period, same cost, a margin per `k` — and they are not. The
acta maximises the margin in `k` shared sections over **one fixed design** of
that cost. `solve_lattice_schedule` maximises over **every** design within the
budget. Same numbers on the surface, two different subjects underneath, and the
values do differ: at `k = 4` and `k = 6` on cost 47, and at `k = 2` on cost 64,
the package reaches a margin the acta's fixed design does not. Quoting both as
one field would have produced a divergence that looks like a bug and is a
mismatch of questions. They are recorded in the manifest's `not_reproducible`
list with that reason instead, along with `PASOS_P55` (an instrument
measurement), `CAIDAS_P55` (a census over 116 universes, not a quantity of one
period) and `FUERA_P55` (a set per period, which no column kind can hold yet).

## Citation

See `CITATION.cff`. Author ORCID: [0009-0001-9116-075X](https://orcid.org/0009-0001-9116-075X).

## License

MIT.

## Field notes

This repository belongs to an audit corpus whose measurements are published as
field notes. Two of them come from the same work:

- [Maths 2.0: why AI needs auditors, not just code generators](https://dev.to/alexisromtz7/maths-20-por-que-la-ia-necesita-auditores-no-solo-generadores-de-codigo-2f29)
- [Mide el jacobiano antes de optimizar: cuántas restricciones tiene de verdad tu objetivo](https://dev.to/alexisromtz7/mide-el-jacobiano-antes-de-optimizar-cuantas-restricciones-tiene-de-verdad-tu-objetivo-3oe1)
