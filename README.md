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
  decided by                 search
```

Below that boundary the formula *is* the price over `U(T)`, and the theorem is
sharp: `16` breaks, `16` triples, no break without a triple and no triple
without a break. Divisibility by `30` is necessary and not sufficient -- `840`
holds `12` and `20` but not `30`, so it is exact.

#### The law behind the boundary, and why the universe decides it

The frontier is a property of the **universe**, and the two answers are eighteen
times apart. Over every divisor, a period is mispriced exactly when one of its
divisors has a **drop-one block**: for `N` let `B(N) = {N // p}` over the
distinct primes `p` of `N`, an antichain by construction whose least common
multiple is `N`. Then

    saving(N) = sum(b // 2 for b in B(N)) - lcm(B(N)) // 2

and the closed form is wrong on `D(T)` **if and only if some `N` dividing `T`
has `saving(N) > 0`**. That quantity is about
`(N / 2) * (sum(1 / p) - 1)` over the primes of `N`, so the criterion is a sum
of reciprocal primes -- and the floors are not noise: `N = 30` has reciprocal
sum `31 / 30`, above one, and still saves nothing, because `6 + 10 + 15 = 31`
buys the same thirty points as the single grid `30`. Checked against the
exhaustive block search on all 1000 even periods to 2000: **0 disagreements**,
same witness and same saving, in `O(tau(T))` instead of a search.

```
$ exact-lattice-quantiles frontier --universe all --max-period 2000
even periods 2..2000 over D(T): 82 of 1000 mispriced, first at T = 60
  witness                    (12, 20, 30)
  closed form                31 points
  one grid 60                30 points
  saving                     1 points
  decided by                 blocks

$ exact-lattice-quantiles blocks --limit 250
divisibility-minimal blocks up to 250: 4
      60  saves   1  (12, 20, 30)
      90  saves   1  (18, 30, 45)
     150  saves   2  (30, 50, 75)
     210  saves  18  (30, 42, 70, 105)
```

| universe | frontier | decided by | witness |
|---|---|---|---|
| `D(T)`, every divisor | **60** | `saving(N) > 0` for some divisor, `O(tau(T))` | `{12, 20, 30}` → `31` / `30` |
| `U(T)`, reduced carriers | **1080** | search; the law is necessary, not sufficient | `{54, 90, 135}` → `139` / `135` |

Two readings of one theory. Divisibility by `60` is sufficient over `D(T)` and
49 of the 82 breaks to 2000 are no multiple of `60`, so the rule reads the
divisors and not the multiples. Reducing to carriers buys **18x** of delay and
destroys the sufficiency: `840` keeps `12` and `20` and loses `30`, and of the
16 breaks in `U(T)` below 8000, nine have witnesses that are no divisor's whole
block. A frontier quoted without its universe is a number about the wrong set.

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

## Exactness

Everything the library computes and returns is `int` or
`fractions.Fraction`. No comparison is decided by rounding and no optimum can
be `1/3 - 1e-16`. `as_dict()` serialises every rational as `"n/d"`; `--json` on
the CLI does the same. The single float in the package is the human-readable
percentage the CLI prints next to the exact fraction — it is never compared,
returned or stored.

## Verification

```
python -m unittest discover -s tests -t .     # 124 tests
exact-lattice-quantiles --self-test           # 13 checks, RC 0
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
# 69 checks over 21 anchors: 69 reproduce, 0 divergent        RC 0
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

## Citation

See `CITATION.cff`. Author ORCID: [0009-0001-9116-075X](https://orcid.org/0009-0001-9116-075X).

## License

MIT.
