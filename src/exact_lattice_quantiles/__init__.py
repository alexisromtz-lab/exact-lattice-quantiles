"""Exact lattice quantiles: provable phase margin for multi-rate sampling.

A periodic system with period ``T`` is sampled by grids -- the divisors of
``T`` above one -- and each grid cuts the phase into ``T / 2`` section classes.
A design is a set of grids; its margin is how far a section can drift before
two divisors are confused.

Two facts make this a solved problem instead of a tuning exercise.

**The ceiling belongs to the instance.**  ``nu(T) = min_{d in U(T)} mass(d)/d``
in closed form, no search.  ``1/3`` is not a law: 167 of the 200 even periods
below 400 exceed it, 30 hit it and 3 fall below.

**The bound is a coverage condition, not a sample.**  "margin at least ``v``"
is "the grids whose blocks are at least ``v`` wide cover the universe"; its
minimum cost is a DP, and when that hypergraph has no cover, ``v`` is
impossible at any price.

Quick start::

    from exact_lattice_quantiles import instance_ceiling, solve_lattice_schedule

    ceiling = instance_ceiling(180)          # 1/3, tied by divisor 6
    schedule = solve_lattice_schedule(180, max_cost=64, k_probes=4, min_margin="1/3")
    schedule.cost                            # 64
    schedule.design                          # (4, 6, 10, 12, 15, 18, 20, 45)
    schedule.witness_margin                  # 1/3

Every number is an ``int`` or a ``fractions.Fraction``.  No floats, no
dependencies.
"""

from __future__ import annotations

from .arithmetic import (divisors, dropped_primes, euler_phi, mass,
                         prime_factorization, prime_product, primes_of, radical,
                         totient_sum)
from .ceiling import (Ceiling, InfeasibleAboveInstanceCeiling,
                      OddPeriodUnsupported, ceiling_by_product,
                      ceiling_by_scan, ceiling_census, instance_ceiling)
from .coverage import (DesignLaw, DriftLaw, best_probes, cover_masks,
                       coverage_kernel, covered_divisors, design_law,
                       design_margin, drift_law, drift_radius, drift_tolerance,
                       guaranteed_margin, margin_ladder, margin_under_drift,
                       minimum_cost_design, minimum_probes, witness_margin)
from .lattice import (TABLE_BUDGET, MaskTableTooLarge, carriers, cell_intervals,
                      cell_width, cell_widths, class_values, grids,
                      grids_by_trial, mask_table, median_carriers,
                      median_oracle_census, period_of, seats, universe,
                      universe_by_median, universe_from_factorization, widths)
from .pricing import (PriceFormCensus, antichain_saving_ceiling,
                      arity_ladder_census,
                      best_block_saving, block_advantage, block_catalogue,
                      block_defect, block_dominance_census, block_fee,
                      block_fee_closed, block_free, block_of, block_saving,
                      band_slack, band_window, carrier_block_of,
                      carrier_block_saving, carrier_cofactor_window,
                      carrier_cofactors, carrier_image, carrier_image_by_law,
                      carrier_image_census,
                      carrier_law_census, carrier_universe_by_law,
                      in_carrier_universe,
                      cofactor_block, crossing_scalars, exactness_frontier,
                      free_margin, free_margin_census,
                      free_margin_cross,
                      excess_parity, extremal_census, fee_frontier_by_blocks,
                      fee_frontier_census, frontier_by_blocks,
                      heaviest_antichain, maximal_divisors,
                      odd_band_census,
                      odd_members, pair_antichain_of, pair_antichain_saving,
                      price_exact, price_form_census,
                      price_form_sweep, price_one_per_maximal, prime_excess,
                      reciprocal_filter_census, reciprocal_prime_sum,
                      scaling_census, scaling_margin,
                      setup_fee_frontier, sign_regime, sign_theorem_census,
                      strong_lemma_frontier, two_adic,
                      worst_antichain_weight)
from .solve import (OK, PROBES_INSUFFICIENT, UNDERCOVERED,
                    InfeasibleWithinBudget, Schedule, ladder,
                    solve_lattice_schedule)

__version__ = "0.27.0"

__all__ = [
    "Ceiling",
    "antichain_saving_ceiling",
    "block_dominance_census",
    "pair_antichain_of",
    "pair_antichain_saving",
    "reciprocal_prime_sum",
    "sign_regime",
    "sign_theorem_census",
    "heaviest_antichain",
    "extremal_census",
    "block_advantage",
    "block_free",
    "cofactor_block",
    "free_margin",
    "free_margin_census",
    "free_margin_cross",
    "carrier_cofactors",
    "carrier_universe_by_law",
    "in_carrier_universe",
    "carrier_law_census",
    "carrier_cofactor_window",
    "carrier_image",
    "carrier_image_by_law",
    "carrier_image_census",
    "scaling_margin",
    "crossing_scalars",
    "scaling_census",
    "arity_ladder_census",
    "block_fee",
    "block_fee_closed",
    "block_defect",
    "setup_fee_frontier",
    "fee_frontier_by_blocks",
    "fee_frontier_census",
    "strong_lemma_frontier",
    "two_adic",
    "worst_antichain_weight",
    "InfeasibleAboveInstanceCeiling",
    "OddPeriodUnsupported",
    "MaskTableTooLarge",
    "TABLE_BUDGET",
    "mask_table",
    "InfeasibleWithinBudget",
    "OK",
    "PROBES_INSUFFICIENT",
    "UNDERCOVERED",
    "PriceFormCensus",
    "Schedule",
    "__version__",
    "best_block_saving",
    "covered_divisors",
    "block_catalogue",
    "block_of",
    "block_saving",
    "band_slack",
    "band_window",
    "odd_band_census",
    "excess_parity",
    "odd_members",
    "prime_excess",
    "reciprocal_filter_census",
    "carrier_block_of",
    "carrier_block_saving",
    "frontier_by_blocks",
    "best_probes",
    "carriers",
    "class_values",
    "median_carriers",
    "median_oracle_census",
    "universe_by_median",
    "grids_by_trial",
    "ceiling_by_product",
    "ceiling_by_scan",
    "ceiling_census",
    "cell_intervals",
    "cell_width",
    "cell_widths",
    "coverage_kernel",
    "design_margin",
    "design_law",
    "DesignLaw",
    "DriftLaw",
    "divisors",
    "drift_law",
    "drift_radius",
    "drift_tolerance",
    "guaranteed_margin",
    "exactness_frontier",
    "dropped_primes",
    "euler_phi",
    "totient_sum",
    "grids",
    "instance_ceiling",
    "margin_under_drift",
    "ladder",
    "margin_ladder",
    "mass",
    "maximal_divisors",
    "period_of",
    "minimum_cost_design",
    "minimum_probes",
    "price_exact",
    "price_form_census",
    "price_form_sweep",
    "price_one_per_maximal",
    "prime_factorization",
    "prime_product",
    "primes_of",
    "radical",
    "seats",
    "solve_lattice_schedule",
    "universe",
    "universe_from_factorization",
    "widths",
    "witness_margin",
]
