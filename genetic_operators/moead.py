#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MOEA/D Genetic Operators
Yachay Tech University - Phage Therapy Group

Contains:
  - CrossoverOperatorMOEAD : SBX-style crossover adapted for discrete protein sequences
  - SelectionOperatorMOEAD : Tchebycheff-based neighbourhood selection/survival
"""

from typing import List, Tuple
import copy
import random
import numpy as np

from genetic_operators.crossoever_operators import uniform_crossover
from genetic_operators.individual import Individual
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.decomposition.tchebicheff import Tchebicheff
from scipy.spatial.distance import cdist

# ═══════════════════════════════════════════════════════════════════════════════
#  SelectionOperatorMOEAD
# ═══════════════════════════════════════════════════════════════════════════════

class SelectionOperatorMOEAD:
    """
    MOEA/D neighbourhood selection and survival operator.

    Encapsulates:
    * Reference-direction weight vectors (Das–Dennis)
    * Neighbourhood structure
    * Tchebycheff decomposition
    * ``parent_selection``  – picks two parents from a sub-problem's neighbourhood
    * ``selection``         – classical MOEA/D survival update

    Parameters
    ----------
    n_obj        : int   – number of objectives
    n_neighbors  : int   – neighbourhood size (T)
    n_partitions : int   – Das–Dennis partitions (controls number of weight vectors)
    """

    def __init__(self, n_obj: int, n_neighbors: int, n_partitions: int = 4):
        self.n_obj       = n_obj
        self.ref_dirs    = get_reference_directions(
            "das-dennis", n_dim=n_obj, n_partitions=n_partitions
        )
        self.decomp      = Tchebicheff()
        dists            = cdist(self.ref_dirs, self.ref_dirs)
        self.neighbors   = np.argsort(dists, axis=1)[:, :n_neighbors]

    # ------------------------------------------------------------------
    #  Parent selection  (for crossover)
    # ------------------------------------------------------------------

    def parent_selection(self, population: List[Individual], sub_problem_idx: int, n_parent_pairs: int = 1
    ) -> List[Tuple[Individual, Individual]]:
        """
        Select parent pairs uniformly at random from a sub-problem neighbourhood.

        Parameters
        ----------
        population       : current population (len == number of sub-problems)
        sub_problem_idx  : index of the sub-problem being evolved
        n_parent_pairs   : number of parent pairs sampled from the neighbourhood
                           (default: 1). It cannot exceed the available
                           neighbourhood size.

        Returns
        -------
        List[(parent_a, parent_b)] – both parents are references into
        ``population``.
        """
        if len(population) < 2:
            raise ValueError("At least two individuals are required for parent selection.")
        if n_parent_pairs < 1:
            raise ValueError("n_parent_pairs must be at least 1.")

        neighbor_ids = [
            j for j in self.neighbors[sub_problem_idx % len(self.neighbors)]
            if j < len(population)
        ]

        if len(neighbor_ids) < 2:
            neighbor_ids = list(range(len(population)))

        if n_parent_pairs > len(neighbor_ids):
            raise ValueError(
                f"n_parent_pairs ({n_parent_pairs}) cannot exceed "
                f"neighbourhood size ({len(neighbor_ids)})."
            )

        return [
            tuple(population[j] for j in random.sample(neighbor_ids, 2))
            for _ in range(n_parent_pairs)
        ]

    # ------------------------------------------------------------------
    #  Survival / replacement  (after mutation)
    # ------------------------------------------------------------------

    def selection(self, population: List[Individual], offspring: List[Individual]) -> List[Individual]:
        """
        Classical MOEA/D survival update with Tchebycheff decomposition.

        Each offspring[i] is compared against the individuals in the
        neighbourhood of sub-problem ``i % len(ref_dirs)``.  An offspring
        replaces a neighbour if its Tchebycheff score is no worse.

        When the same offspring wins multiple neighbourhood slots, subsequent
        replacements receive a ``copy.deepcopy`` of the offspring so that
        ``move_population`` can detect duplicates by path and copy the file
        instead of moving it.

        Parameters
        ----------
        population : List[Individual]  – current population
        offspring  : List[Individual]  – one offspring per sub-problem

        Returns
        -------
        List[Individual]  – updated population (same length)
        """
        F_pop  = np.array([ind.F for ind in population])
        z_star = np.min(F_pop, axis=0)   # ideal point (minimisation space)

        new_population: List[Individual] = list(population)
        replaced_by: dict = {}           # offspring_idx → first slot won

        for i, child in enumerate(offspring):
            f_child      = child.F
            neighbor_ids = self.neighbors[i % len(self.neighbors)]

            for j in neighbor_ids:
                if j >= len(new_population):
                    continue

                w_j = self.ref_dirs[j]

                score_current = self.decomp.do(
                    new_population[j].F[None, :],
                    weights=w_j[None, :],
                    ideal_point=z_star,
                )[0, 0]

                score_child = self.decomp.do(
                    f_child[None, :],
                    weights=w_j[None, :],
                    ideal_point=z_star,
                )[0, 0]

                if score_child <= score_current:
                    if i not in replaced_by:
                        # First slot this offspring wins – use the object directly
                        replaced_by[i] = j
                        new_population[j] = child
                    else:
                        # Same offspring wins again – deep-copy, keep same pdb path
                        # so move_population detects the duplicate and copies the file
                        clone     = copy.deepcopy(child)
                        clone.pdb = child.pdb
                        new_population[j] = clone

        return new_population