import random

import ollama
from typing import Dict, List, Optional, Sequence, Tuple

from genetic_operators.individual import Individual


IndividualRecord = Dict[str, object]

_ollama_client = None
_ollama_client_host = None


def _get_ollama_client(host: str) -> ollama.Client:
    global _ollama_client, _ollama_client_host
    if _ollama_client is None or _ollama_client_host != host:
        _ollama_client = ollama.Client(host=host)
        _ollama_client_host = host
    return _ollama_client


def _normalize_sequence_value(sequence_value: object, label: str) -> str:
    if isinstance(sequence_value, str):
        sequence = "".join(sequence_value.strip().split())
    elif isinstance(sequence_value, (list, tuple)):
        try:
            sequence = "".join(str(token).strip() for token in sequence_value)
        except Exception as exc:
            raise ValueError(f"{label} must be a sequence-like value.") from exc
    else:
        raise ValueError(f"{label} must be a string or list/tuple of residues.")

    if not sequence:
        raise ValueError(f"{label} must include a non-empty sequence.")
    return sequence


def _get_sequence(individual: IndividualRecord, label: str) -> str:
    sequence = individual.get("sequence")
    return _normalize_sequence_value(sequence, f"{label}['sequence']")


def _get_fitness(individual: IndividualRecord, label: str) -> Tuple[float, ...]:
    fitness = individual.get("fitness")
    if isinstance(fitness, (int, float)):
        return (float(fitness),)
    if not isinstance(fitness, (tuple, list)):
        raise ValueError(f"{label} must include 'fitness' as a tuple/list of numbers.")
    try:
        return tuple(float(value) for value in fitness)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} fitness must contain only numeric values.") from exc


def _format_primary_fitness(fitness: Sequence[float]) -> str:
    return f"dG_separated/dSASAx100={float(fitness[0])}"


def _format_population_context(individuals: Sequence[IndividualRecord]) -> str:
    rows = []
    for idx, individual in enumerate(individuals):
        sequence = _get_sequence(individual, f"individuals[{idx}]")
        fitness = _get_fitness(individual, f"individuals[{idx}]")
        if not fitness:
            raise ValueError(f"individuals[{idx}] must include at least one fitness value.")
        rows.append(
            f"Individual {idx}: sequence={sequence}; "
            f"fitness=({_format_primary_fitness(fitness)})"
        )
    return "\n".join(rows)


def llm(
    parent1: IndividualRecord,
    parent2: IndividualRecord,
    individuals: List[IndividualRecord],
    objective_description: str,
    sequence_initial: Optional[str] = None,
    model_name: str = "gemma4:e4b",
    temperature: float = 0.3,
    ollama_host: str = "http://localhost:11435",
) -> str:
    parent1_seq = _get_sequence(parent1, "parent1")
    parent2_seq = _get_sequence(parent2, "parent2")
    parent1_fitness = _get_fitness(parent1, "parent1")
    parent2_fitness = _get_fitness(parent2, "parent2")

    if len(parent1_seq) != len(parent2_seq):
        raise ValueError("parent1 and parent2 sequences must have the same length.")
    if not parent1_fitness:
        raise ValueError("parent1 must include at least one fitness value.")
    if not parent2_fitness:
        raise ValueError("parent2 must include at least one fitness value.")

    normalized_individuals: List[IndividualRecord] = []
    for idx, individual in enumerate(individuals):
        sequence = _get_sequence(individual, f"individuals[{idx}]")
        if len(sequence) != len(parent1_seq):
            raise ValueError(
                f"individuals[{idx}] sequence length ({len(sequence)}) must match "
                f"parent sequence length ({len(parent1_seq)})."
            )
        normalized_individuals.append({**individual, "sequence": sequence})

    client = _get_ollama_client(ollama_host)

    reference_sequence = (
        _normalize_sequence_value(sequence_initial, "sequence_initial")
        if sequence_initial is not None
        else None
    )
    if reference_sequence is not None and len(reference_sequence) != len(parent1_seq):
        raise ValueError(
            f"sequence_initial length ({len(reference_sequence)}) must match "
            f"parent sequence length ({len(parent1_seq)})."
        )

    reference = (
        f"Reference interface sequence:\n{reference_sequence}\n\n"
        if reference_sequence
        else ""
    )
    system_prompt = (
        "You are an expert protein engineer specializing in protein-protein interface "
        "optimization. Your task is to perform a biologically informed crossover "
        "between two parent protein sequences to produce a child with improved "
        "binding affinity."
    )

    user_prompt = (
        f"Perform one crossover event using INTERFACE sequences only.\n"
        f"Objective: {objective_description}\n"
        f"The first and only fitness shown is dG_separated/dSASAx100; lower values are better "
        f"and indicate a more favorable protein-protein interaction. Minimize this value.\n\n"
        f"{reference}"
        f"Parent 1 interface sequence:\n{parent1_seq}\n"
        f"Fitness: ({_format_primary_fitness(parent1_fitness)})\n\n"
        f"Parent 2 interface sequence:\n{parent2_seq}\n"
        f"Fitness: ({_format_primary_fitness(parent2_fitness)})\n\n"
        f"Using your expertise in protein-protein interface optimization, choose the best "
        f"crossover between the two parents. Prefer "
        f"recombining residues from the parents at interface positions, but you may choose "
        f"a different standard amino acid at an interface position if it is strongly "
        f"biochemically justified.\n"
        f"Do not return either parent unchanged.\n"
        f"Keep exactly the same interface-sequence length: {len(parent1_seq)}.\n"
        f"Output ONLY the child interface sequence, with no explanation."
    )

    response = client.chat(
        model=model_name,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
        options={"temperature": temperature},
    )

    child_sequence = response["message"]["content"].strip()
    if len(child_sequence) != len(parent1_seq):
        raise ValueError(
            f"LLM returned a sequence with length {len(child_sequence)}, "
            f"expected {len(parent1_seq)}."
        )
    return child_sequence


def _record_from_individual(individual: Individual, sequence: str) -> IndividualRecord:
    fitness = individual.fitness
    if fitness is None:
        fitness = individual.F
    if fitness is None:
        raise ValueError(f"Individual {individual.id} does not have fitness values.")

    return {
        "sequence": sequence,
        "fitness": (float(fitness[0]),),
    }


def _clean_child_sequence(child_sequence: str) -> str:
    return "".join(child_sequence.strip().upper().split())

def _primary_fitness(ind: Individual) -> float:
    fitness = ind.fitness if ind.fitness is not None else ind.F
    if fitness is None:
        raise ValueError(f"Individual {ind.id} does not have fitness values for LLM crossover.")
    return float(fitness[0])


def llm_crossover(
    parent_a: Individual, parent_b: Individual, population: List[Individual], sequence_initial: Optional[str] = None, model_name: str = "gemma4:e4b", temperature: float = 0.3,
    ollama_host: str = "http://ollama:11434"
) -> Tuple[List[int], List[str]]:
    
    """Return interface-level crossover edits proposed by the LLM."""
    seq_a = parent_a.sequence()
    seq_b = parent_b.sequence()

    objective_description = (
            "Improve the protein-protein interface by minimizing "
            "dG_separated/dSASAx100, the first fitness value passed to the LLM."
        )
    
    individuals = [
            {
                "sequence": ind.sequence(),
                "fitness": (_primary_fitness(ind),),
            }
            for ind in population
        ]

    if len(seq_a) != len(seq_b):
        raise ValueError("parent_a and parent_b interface sequences must have the same length.")

    parent1 = _record_from_individual(parent_a, seq_a)
    parent2 = _record_from_individual(parent_b, seq_b)

    child_sequence = _clean_child_sequence(
        llm(parent1=parent1, parent2=parent2, individuals=individuals, objective_description=objective_description,
            sequence_initial=sequence_initial, model_name=model_name, temperature=temperature, ollama_host=ollama_host
        )
    )

    if len(child_sequence) != len(seq_a):
        raise ValueError(
            f"LLM child interface sequence length ({len(child_sequence)}) must match "
            f"parent interface sequence length ({len(seq_a)})."
        )

    allowed_aas = set("ACDEFGHIKLMNPQRSTVWY")
    invalid_aas = sorted(set(child_sequence) - allowed_aas)
    if invalid_aas:
        raise ValueError(f"LLM child sequence contains invalid amino acids: {invalid_aas}")

    if child_sequence == seq_a or child_sequence == seq_b:
        raise ValueError("LLM child sequence must not be identical to either parent.")

    seq_indices: List[int] = []
    child_aas: List[str] = []

    for idx, aa_child in enumerate(child_sequence):
        if seq_a[idx] == aa_child:
            continue
        seq_indices.append(idx)
        child_aas.append(aa_child)

    return seq_indices, child_aas


def uniform_crossover(parent_a: Individual, parent_b: Individual, crossover_prob: float = 0.5) -> Tuple[List[int], List[str]]:
    seq_a = parent_a.sequence()   # list[str]
    seq_b = parent_b.sequence()
    seq_indices: List[int] = []
    child_aas: List[str] = []

    for idx, (aa_a, aa_b) in enumerate(zip(seq_a, seq_b)):
        if aa_a == aa_b:
            continue  # identical - no structural change needed
        # Uniform crossover: pick from parent_b with probability crossover_prob
        if random.random() < crossover_prob:
            seq_indices.append(idx)
            child_aas.append(aa_b)

    return seq_indices, child_aas
