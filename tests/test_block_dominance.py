"""The lemma the law rests on, checked period by period.

The law says a period is mispriced exactly when one of its divisors has a
drop-one block that saves.  The direction "a saving block misprices" needs
nothing: the block is a cover, so the closed form that prices its parts is
overcharging.  The other direction needs a statement about covers that are
**not** blocks:

    no antichain with least common multiple ``L`` saves more than the block of
    ``L``,

because a cover made of several grids can always be folded into the single
grid their least common multiple names, and that grid's own block is the
cheapest way to buy it whole.  With the lemma, a cover built from parts never
beats paying the parts, so no block saving means the closed form is optimal.

Read as a ranking -- the block saves *most* among the antichains it shares a
least common multiple with -- the lemma is false, and ``1049`` is the prime
that says so: on the primorial of the first 176 primes the pair antichain wins
(see ``test_strong_lemma_frontier``).  What the law uses is the weaker **sign**
reading, that no antichain saves where no block does, and that one is proved in
``pricing.sign_theorem_census`` from the prime factorisation alone: ``sum 1/p
<= 1`` forces every antichain to weigh at most ``L * sum 1/p``, which leaves the
saving at ``<= 1/2`` and integrality finishes. The census below is therefore
the check on the implementation, not the evidence for the claim.

What the tests pin is that it is checked rather than assumed:
``block_dominance_census`` recomputes the ceiling for every least common
multiple in range, and a counterexample arrives as a failing test instead of a
rumour.
"""

import unittest

from exact_lattice_quantiles.pricing import (antichain_saving_ceiling,
                                             block_catalogue,
                                             block_dominance_census, block_of,
                                             block_saving)

# Quoted in the README: the ceiling of a period where a second saving set
# exists that is not a block -- and loses to it.
SECOND_BEST = 1050


class TheBlockCarriesEverySaving(unittest.TestCase):
    def test_the_ceiling_names_the_block_itself(self):
        for period, saving, block in block_catalogue(2000):
            ceiling, witness = antichain_saving_ceiling(period)
            self.assertEqual(ceiling, saving, f"ceiling of {period}")
            self.assertEqual(tuple(sorted(witness)), block,
                             f"witness of {period}")

    def test_the_census_finds_no_violation(self):
        census = block_dominance_census(2000)
        self.assertEqual(census["violations"], [])
        self.assertEqual(census["periods"], 1999)
        self.assertEqual(len(census["violations"]), 0)

    def test_the_census_counts_the_periods_where_a_saving_is_in_play(self):
        census = block_dominance_census(2000)
        in_play = sum(1 for n in range(2, 2001) if block_saving(n) > 0)
        self.assertEqual(census["in_play"], in_play)
        self.assertEqual(in_play, 82)

    def test_a_second_saving_set_loses_to_the_block(self):
        ceiling, witness = antichain_saving_ceiling(SECOND_BEST)
        self.assertEqual(ceiling, 92)
        self.assertEqual(tuple(sorted(witness)), block_of(SECOND_BEST)[0])
        # {42, 150, 350, 525} also saves -- 8 points, with 42 no 1050 / p --
        # so the lemma is about the maximum and not about the shape.
        self.assertEqual(sum(d // 2 for d in (42, 150, 350, 525)) - 1050 // 2, 8)

    def test_the_sign_reading_holds_where_the_ranking_reading_does(self):
        census = block_dominance_census(2000)
        self.assertEqual(census["sign_violations"], [])
        self.assertEqual(census["sign_verdict"],
                         "no antichain saves where no block does")

    def test_the_reciprocal_sum_is_not_the_criterion(self):
        """``1/2 + 1/3 + 1/5`` is above one and still need not save a point.

        ``30`` prices its block at ``7 + 5 + 3 = 15`` against ``30 // 2`` for the
        one grid; ``42`` and ``66`` lose more.  ``60`` and ``90`` are the same
        shape with a square factor and do save one point each -- which is why
        the boundary is not "four prime factors" either.
        """
        for period in (12, 30, 36, 42, 66, 70, 78):
            self.assertEqual(antichain_saving_ceiling(period), (0, ()),
                             f"{period} should not save")
        for period in (60, 90):
            self.assertEqual(antichain_saving_ceiling(period)[0], 1)


if __name__ == "__main__":
    unittest.main()
