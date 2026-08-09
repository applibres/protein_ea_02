#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NSGA-III Genetic Operators
Yachay Tech University - Phage Therapy Group

Contains:
  - SelectionOperatorNSGA3 : reference-direction survival + global parent selection
"""

from typing import List, Tuple
import random

import numpy as np
from scipy.spatial.distance import cdist

from genetic_operators.individual import Individual
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.core.population import Population
from pymoo.core.problem import Problem
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from pymoo.util.ref_dirs import get_reference_directions


class SelectionOperatorNSGA3:
    """
    NSGA-III selection/survival operator.

    Encapsulates:
    * Reference-direction vectors (Das-Dennis)
    * Global parent selection (no neighbourhood)
    * NSGA-III reference-direction survival

    Parameters
    ----------
    n_obj        : int  - number of objectives
    n_partitions : int  - Das-Dennis partitions (controls # reference directions)
    pop_size     : int  - number of individuals to keep after survival
    """

    def __init__(self, n_obj: int, n_partitions: int = 4, pop_size: int = 100):
        self.n_obj = int(n_obj)
        self.pop_size = int(pop_size)

        if self.n_obj < 1:
            raise ValueError("n_obj must be >= 1.")
        if n_partitions < 1:
            raise ValueError("n_partitions must be >= 1.")
        if self.pop_size < 2:
            raise ValueError("pop_size must be >= 2 for NSGA-III mating/survival.")

        all_ref_dirs = get_reference_directions(
            "das-dennis", n_dim=self.n_obj, n_partitions=int(n_partitions)
        )

        if self.pop_size > len(all_ref_dirs):
            raise ValueError(
                f"Requested pop_size={self.pop_size}, but Das-Dennis generated only "
                f"{len(all_ref_dirs)} directions. Increase n_partitions."
            )

        self.ref_dirs = self._select_active_ref_dirs(all_ref_dirs, self.pop_size)

        # NSGA-III survival helper from pymoo.
        self._algorithm = NSGA3(pop_size=self.pop_size, ref_dirs=self.ref_dirs)
        self._survival = self._algorithm.survival

        # Minimal problem object required by pymoo survival API.
        self._problem = self._build_problem(self.n_obj)

    @staticmethod
    def _build_problem(n_obj: int):
        # Compatibility with pymoo versions where constraints arg name changed.
        try:
            return Problem(
                n_var=1,
                n_obj=n_obj,
                n_ieq_constr=0,
                n_eq_constr=0,
                xl=0,
                xu=1,
            )
        except TypeError:
            return Problem(n_var=1, n_obj=n_obj, n_constr=0, xl=0, xu=1)

    @staticmethod
    def _select_active_ref_dirs(ref_dirs: np.ndarray, n_subproblems: int) -> np.ndarray:
        """
        Select a deterministic well-spread subset of reference directions.

        Uses greedy farthest-point sampling so active directions preserve
        coverage of the original Das-Dennis simplex directions.
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

    def parent_selection(
        self,
        population: List[Individual],
        sub_problem_idx: int,
        n_parent_pairs: int = 1,
    ) -> List[Tuple[Individual, Individual]]:
        """
        Select parent pairs from the whole population using binary tournament
        on non-dominated rank.

        The ``sub_problem_idx`` argument is accepted for interface compatibility
        with MOEA/D parent selection, but it is not used in NSGA-III.
        """
        _ = sub_problem_idx

        if len(population) != self.pop_size:
            raise ValueError(
                f"Population size ({len(population)}) must match configured "
                f"pop_size ({self.pop_size})."
            )
        if len(population) < 2:
            raise ValueError("At least two individuals are required for parent selection.")
        if n_parent_pairs < 1:
            raise ValueError("n_parent_pairs must be at least 1.")

        F_population = np.array([ind.F for ind in population], dtype=float)
        if F_population.ndim != 2 or F_population.shape[1] != self.n_obj:
            raise ValueError(
                f"Each individual must provide F with {self.n_obj} objectives; "
                f"received shape {F_population.shape}."
            )

        # Binary tournament on non-dominated rank (global mating pool).
        fronts = NonDominatedSorting().do(F_population)
        ranks = np.full(len(population), fill_value=np.iinfo(np.int32).max, dtype=np.int32)
        for rank, front in enumerate(fronts):
            ranks[np.asarray(front, dtype=int)] = rank

        population_ids = list(range(len(population)))

        def tournament_pick() -> int:
            a, b = random.sample(population_ids, 2)
            rank_a = int(ranks[a])
            rank_b = int(ranks[b])
            if rank_a < rank_b:
                return a
            if rank_b < rank_a:
                return b
            return random.choice((a, b))

        parent_pairs: List[Tuple[Individual, Individual]] = []
        for _ in range(n_parent_pairs):
            idx_a = tournament_pick()
            idx_b = tournament_pick()
            if idx_a == idx_b:
                alternatives = [idx for idx in population_ids if idx != idx_a]
                idx_b = random.choice(alternatives)
            parent_pairs.append((population[idx_a], population[idx_b]))

        return parent_pairs

    def selection(self, population: List[Individual], offspring: List[Individual]) -> List[Individual]:
        """
        NSGA-III survival update.

        Combines current population and offspring and keeps exactly ``pop_size``
        individuals according to pymoo's NSGA-III survival operator.
        """
        if len(population) != self.pop_size:
            raise ValueError(
                f"Population size ({len(population)}) must match configured "
                f"pop_size ({self.pop_size})."
            )
        if len(offspring) != self.pop_size:
            raise ValueError(
                f"Offspring size ({len(offspring)}) must match configured "
                f"pop_size ({self.pop_size})."
            )

        combined = list(population) + list(offspring)
        F_combined = np.array([ind.F for ind in combined], dtype=float)

        if F_combined.ndim != 2 or F_combined.shape[1] != self.n_obj:
            raise ValueError(
                f"Each individual must provide F with {self.n_obj} objectives; "
                f"received shape {F_combined.shape}."
            )

        pop_pymoo = Population.new("F", F_combined)
        pop_pymoo.set("idx", np.arange(len(combined), dtype=int))

        selected = self._do_survival(pop_pymoo)
        selected_idx = selected.get("idx")

        if selected_idx is None:
            raise RuntimeError("NSGA-III survival did not preserve candidate indices.")
        if len(selected_idx) != self.pop_size:
            raise RuntimeError(
                f"NSGA-III survival returned {len(selected_idx)} individuals; "
                f"expected {self.pop_size}."
            )

        return [combined[int(i)] for i in selected_idx]

    def _do_survival(self, pop_pymoo: Population) -> Population:
        """Version-tolerant adapter for pymoo survival signatures."""
        try:
            return self._survival.do(
                self._problem,
                pop_pymoo,
                n_survive=self.pop_size,
                algorithm=self._algorithm,
            )
        except TypeError:
            try:
                return self._survival.do(
                    self._problem,
                    pop_pymoo,
                    n_survive=self.pop_size,
                )
            except TypeError:
                return self._survival.do(self._problem, pop_pymoo, self.pop_size)
