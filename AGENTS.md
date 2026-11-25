# Working notes for this repository

Public, standalone package. Nothing here depends on any private project, and
nothing private may be copied into it: no local paths, no credentials, no
internal vocabulary. If a claim needs a private document to be understood, it
does not belong here.

## Rules that decide whether a change is acceptable

1. **Exactness.** Every number the library returns is an `int` or a
   `fractions.Fraction`. Floating point inside a computation is a defect, not a
   shortcut — it makes an optimum a rounding artefact. The one permitted
   conversion is the CLI's percentage column, which is display only.
2. **Zero dependencies.** Runtime and tests use the standard library only.
   Adding a dependency needs a reason that survives review.
3. **A number in the README is a claim.** Every figure quoted in prose must be
   reproducible by a command in this repository (`ladder`, `census`, `ceiling`,
   `solve`, or a test). If it cannot be produced, delete it — do not soften it.
4. **A ceiling above the instance's is not a bug to fix by searching.**
   `InfeasibleAboveInstanceCeiling` (rc 1) and `InfeasibleWithinBudget` (rc 2)
   stay distinct: one cannot be bought, the other can.
5. **Two budgets are never merged.** Grid cost and shared sections are separate
   inputs; a design that cannot be served is reported
   (`PROBES_INSUFFICIENT`), never silently downgraded to a weaker margin.
6. **Say which margin you mean.** `design_margin`, `witness_margin`,
   `robust_margin` (== `drift_radius`) and `drift_tolerance` are four
   quantities. `robust_margin` is exact for *coverage* and is not a guarantee:
   it can exceed the witness margin because drift may push a section into a
   block it never captured. The margin is what `drift_tolerance` protects, and
   it is protected under zero drift only.

## Commands

```
python -m unittest discover -s tests -t .        # full suite
PYTHONPATH=src python -m exact_lattice_quantiles.cli --self-test
PYTHONPATH=src python -m exact_lattice_quantiles.cli ladder 180
PYTHONPATH=src python -m exact_lattice_quantiles.cli drift 180 --max-cost 64 --probes 4
PYTHONPATH=src python -m exact_lattice_quantiles.cli verify-acta \
  src/exact_lattice_quantiles/manifests/acta_margen_secciones_20261125.json
```

A quote of this package kept somewhere else is recomputed by `verify-acta`
against a manifest of anchors: one `Check` per field, so a divergence names the
quantity and not the file. Editing a manifest is a claim that the new number is
the intended one — never a way to turn a red gate green.

The suite runs in about half a minute; the slow tests are the ones that
compare the factorized ceiling against the full scan on every even period to
400 and the coverage DP against brute force. Keep those: they are the reason
the closed forms are trusted.

## Adding a result

A new closed form enters with three things in the same change: the function,
the test that pins it to numbers already measured, and the reason the naive
version was wrong (the docstring names the defect it avoids — the even-grid
filter that halved odd divisors, or the `T = 2` mass degeneracy, are the two
that shipped). If the naive version was never wrong, there is nothing to
document and the function is probably unnecessary.

## Style

Standard library, type annotations on public functions, no framework. Module
docstrings explain *why* a formula is what it is; line comments only where the
code would otherwise look wrong.
