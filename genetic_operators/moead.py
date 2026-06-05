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
#  CrossoverOperatorMOEAD
# ═══════════════════════════════════════════════════════════════════════════════

class CrossoverOperatorMOEAD:
    """
    Uniform crossover for discrete protein sequences in MOEA/D.

    Given two parent individuals (parent_a, parent_b), it produces ONE child
    by sampling each position from parent_a or parent_b with probability
    ``crossover_prob`` (default 0.5).  Positions where both parents agree are
    always inherited unchanged.

    The child PDB is written to ``output_file`` by the caller (prot_problem).
    This class only computes *which* positions differ between the two parents
    and *which* amino acids the child should receive at those positions, so
    that ``prot_problem.crossover_population`` can apply the structural changes
    via Rosetta.

    Parameters
    ----------
    crossover_prob : float
        Probability of inheriting an allele from parent_b (vs parent_a).
        Default is 0.5 (uniform crossover).
    """

    def __init__(self, crossover_prob: float = 0.5):
        self.crossover_prob = crossover_prob

    # ------------------------------------------------------------------
    def get_crossover_args(
        self,
        parent_a: Individual,
        parent_b: Individual,
        output_file: str,
    ) -> Tuple[str, str, List[int], List[str]]:
        """
        Compute the arguments needed by ``prot_problem.crossover_population``.

        Parameters
        ----------
        parent_a : Individual
            Structural donor (its PDB is used as the base structure).
        parent_b : Individual
            Sequence donor for the swapped positions.
        output_file : str
            Destination PDB path for the child.

        Returns
        -------
        tuple
            ``(pdb_file, output_file, positions_to_mutate, aminoacids_to_place)``
            ready to be passed as a single element of the ``args`` list to
            ``prot_problem.crossover_population``.

            *positions_to_mutate* uses **1-based interface indices** (the index
            in the interface residue list, NOT Rosetta absolute numbering).
            ``prot_problem._crossover_worker`` already maps these through
            ``get_individual_seq`` so the caller must pass the absolute Rosetta
            positions; see note below.

        Notes
        -----
        ``prot_problem.crossover_population`` expects absolute Rosetta residue
        positions.  The interface positions stored in ``Individual.sequence()``
        are in the same order as ``aa_pos_list``, so the caller
        (``pymoo_sga_protein``) must resolve the indices to absolute positions
        using ``self.my_protein_problem.aa_pos_list`` before building the args.
        This method returns **sequence-level indices** (0-based) so the caller
        can do that mapping.

        Returns a 4-tuple:
            pdb_file          – parent_a PDB (structural base)
            output_file       – child output path
            seq_indices       – 0-based positions in the interface sequence
                                where parent_b's allele was chosen
            child_aas         – single-letter AA codes for those positions
        """
        seq_a = parent_a.sequence()   # list[str]
        seq_b = parent_b.sequence()   # list[str]

        seq_indices: List[int] = []
        child_aas:   List[str] = []

        for idx, (aa_a, aa_b) in enumerate(zip(seq_a, seq_b)):
            if aa_a == aa_b:
                continue  # identical – no structural change needed
            # Uniform crossover: pick from parent_b with probability crossover_prob
            if random.random() < self.crossover_prob:
                seq_indices.append(idx)
                child_aas.append(aa_b)

        return parent_a.pdb, output_file, seq_indices, child_aas

    # ------------------------------------------------------------------
    def make_child_individual(
        self,
        parent_a: Individual,
        parent_b: Individual,
        child_aa: List[str],
        output_file: str,
    ) -> Individual:
        """
        Build an ``Individual`` shell for the child (before evaluation).

        Parameters
        ----------
        parent_a : Individual
        parent_b : Individual
        child_aa  : list[str]
            Full interface AA sequence of the child (length == len(parent_a.sequence())).
        output_file : str

        Returns
        -------
        Individual
        """
        child = Individual(list(child_aa))
        child.pdb    = output_file
        child.father = f"{parent_a.id}+{parent_b.id}"
        child.nmut   = sum(
            1 for a, b in zip(child_aa, parent_a.sequence()) if a != b
        )
        return child


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

    def parent_selection(
        self,
        population: List[Individual],
        sub_problem_idx: int,
    ) -> Tuple[Individual, Individual]:
        """
        Select two parents from the neighbourhood of sub-problem ``sub_problem_idx``.

        The first parent is ``population[sub_problem_idx]`` itself; the second
        is a random individual from its neighbourhood (excluding itself).

        Parameters
        ----------
        population       : current population (len == number of sub-problems)
        sub_problem_idx  : index of the sub-problem being evolved

        Returns
        -------
        (parent_a, parent_b) – both are references into ``population``
        """
        neighbor_ids = self.neighbors[sub_problem_idx % len(self.neighbors)]

        # Exclude sub_problem_idx itself for parent_b
        candidates = [j for j in neighbor_ids if j != sub_problem_idx and j < len(population)]

        if not candidates:
            # Fallback: random individual from the whole population
            candidates = [j for j in range(len(population)) if j != sub_problem_idx]

        parent_a = population[sub_problem_idx]
        parent_b = population[random.choice(candidates)]
        return parent_a, parent_b

    # ------------------------------------------------------------------
    #  Survival / replacement  (after mutation)
    # ------------------------------------------------------------------

    def selection(
        self,
        population: List[Individual],
        offspring: List[Individual],
    ) -> List[Individual]:
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