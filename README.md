# exact-lattice-quantiles

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

## Four claims, each reproducible by a command

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
| `{12, 18, 20, 45}` | 47 | 1/6 | 8 | 2/33 |
| `{4, 6, 10, 12, 15, 18, 20, 45}` | 64 | 1/3 | 4 | 1/3 |

The cheap design is the expensive one to serve: `1/6` needs eight sections
while the maximum margin needs four, because a design whose blocks are only
just wide enough forces each section to serve a single divisor. At `k = 4` on
the cheap design no set of sections can hold more than `2/33`, and no `k` at
all reaches past `1/6`:

```python
solve_lattice_schedule(180, max_cost=47, k_probes=4, min_margin="1/6")
# verdict = "PROBES_INSUFFICIENT", probes_needed = 8, witness_margin = 2/33
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

So the library never prices with the formula: `price_exact` is a subset DP over
the coverage masks, and `price_form_census` publishes how often the shortcut
agrees. Restricting the universe to `U(T)` is what makes the shortcut look
right, and that restriction is a modelling decision, not a measurement --
`price_form_census(period, reduced=False)` is the row that shows it.

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

## Two margins, and which one to quote

`Schedule` carries both, and they are not the same quantity:

- `witness_margin` — what the chosen sections actually hold, with each divisor
  leaning on its own grid/section pair. **Quote this when a schedule is being
  promised.**
- `robust_margin` — the common drift `delta` absorbed by every section at once.
  It is **at least** the witness margin, because it lets two different pairs
  absorb the two drift directions and does not check that the union of
  admissible `delta` stays an interval. It is therefore *optimistic*; quote it
  when the question is how much clock error the design was built around.

`design_margin` is the third view: what the grid set guarantees before any
section is placed. The hierarchy `design_margin >= witness_margin` and
`robust_margin >= witness_margin` is asserted in the test suite.

## Exactness

Everything the library computes and returns is `int` or
`fractions.Fraction`. No comparison is decided by rounding and no optimum can
be `1/3 - 1e-16`. `as_dict()` serialises every rational as `"n/d"`; `--json` on
the CLI does the same. The single float in the package is the human-readable
percentage the CLI prints next to the exact fraction — it is never compared,
returned or stored.

## Verification

```
python -m unittest discover -s tests -t .     # 58 tests
exact-lattice-quantiles --self-test           # RC 0
```

The suite anchors each number quoted above to the function that produces it,
including the cases where the first implementation was wrong: filtering the
scan to even grids halves the ceiling of every odd divisor; `T = 2` needs its
degenerate branch; a design that does not seat every divisor has margin `0`
rather than raising; and the coverage DP is compared against brute force.

## Citation

See `CITATION.cff`. Author ORCID: [0009-0001-9116-075X](https://orcid.org/0009-0001-9116-075X).

## License

MIT.
