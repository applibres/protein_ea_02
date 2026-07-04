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
    nr           : int   – maximum number of neighbour replacements per offspring
    n_subproblems: int   – active MOEA/D subproblems (should match population size)
    """

    def __init__(
        self,
        n_obj: int,
        n_neighbors: int,
        n_partitions: int = 4,
        nr: int = 2,
        n_subproblems: int = None,
    ):
        self.n_obj = n_obj
        all_ref_dirs = get_reference_directions(
            "das-dennis", n_dim=n_obj, n_partitions=n_partitions
        )

        if n_subproblems is None:
            n_subproblems = len(all_ref_dirs)
        n_subproblems = int(n_subproblems)
        if n_subproblems < 1:
            raise ValueError("n_subproblems must be >= 1.")
        if n_subproblems > len(all_ref_dirs):
            raise ValueError(
                f"Requested n_subproblems={n_subproblems}, but Das-Dennis generated "
                f"only {len(all_ref_dirs)} directions. Increase n_partitions."
            )

        self.ref_dirs = self._select_active_ref_dirs(all_ref_dirs, n_subproblems)
        self.n_subproblems = len(self.ref_dirs)
        self.decomp = Tchebicheff()

        neighborhood_size = min(max(1, int(n_neighbors)), self.n_subproblems)
        if self.n_subproblems >= 2:
            neighborhood_size = max(2, neighborhood_size)

        dists = cdist(self.ref_dirs, self.ref_dirs)
        self.neighbors = np.argsort(dists, axis=1)[:, :neighborhood_size]
        self.nr = max(1, int(nr))

    @staticmethod
    def _select_active_ref_dirs(ref_dirs: np.ndarray, n_subproblems: int) -> np.ndarray:
        """
        Select a well-spread deterministic subset of reference directions.

        Uses greedy farthest-point sampling so active subproblems preserve
        coverage of the original Das–Dennis simplex directions.
        """
        n_total = len(ref_dirs)
        if n_subproblems == n_total:
            return ref_dirs.copy()

        selected_indices = [0]
        selected_mask = np.zeros(n_total, dtype=bool)
        selected_mask[0] = True
        min_dist = cdist(ref_dirs, ref_dirs[[0]])[:, 0]

        while len(selected_indices) < n_subproblems:
            candidate_scores = min_dist.copy()
            candidate_scores[selected_mask] = -np.inf
            next_idx = int(np.argmax(candidate_scores))
            selected_indices.append(next_idx)
            selected_mask[next_idx] = True

            dist_to_new = cdist(ref_dirs, ref_dirs[[next_idx]])[:, 0]
            min_dist = np.minimum(min_dist, dist_to_new)

        return ref_dirs[np.array(selected_indices, dtype=int)]

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
        if len(population) != self.n_subproblems:
            raise ValueError(
                f"Population size ({len(population)}) must match number of active "
                f"subproblems ({self.n_subproblems})."
            )
        if len(population) < 2:
            raise ValueError("At least two individuals are required for parent selection.")
        if n_parent_pairs < 1:
            raise ValueError("n_parent_pairs must be at least 1.")

        neighbor_ids = self.neighbors[sub_problem_idx % self.n_subproblems].tolist()
        if len(neighbor_ids) < 2:
            raise RuntimeError(
                "Neighbourhood size must be >= 2 when population has at least two individuals."
            )

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
        if len(population) != self.n_subproblems:
            raise ValueError(
                f"Population size ({len(population)}) must match number of active "
                f"subproblems ({self.n_subproblems})."
            )
        if len(offspring) != self.n_subproblems:
            raise ValueError(
                f"Offspring size ({len(offspring)}) must match number of active "
                f"subproblems ({self.n_subproblems})."
            )

        F_pop  = np.array([ind.F for ind in population])
        z_star = np.min(F_pop, axis=0)   # ideal point (minimisation space)

        new_population: List[Individual] = list(population)
        replaced_by: dict = {}           # offspring_idx → first slot won

        for i, child in enumerate(offspring):
            f_child      = child.F
            # Classical MOEA/D: update the ideal point incrementally
            # with every new offspring before neighborhood replacement.
            z_star       = np.minimum(z_star, f_child)
            neighbor_ids = self.neighbors[i % self.n_subproblems]
            replacements = 0

            for j in neighbor_ids:
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
                    replacements += 1
                    if replacements >= self.nr:
                        break

        return new_population
