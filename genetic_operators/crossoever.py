import ollama
from typing import List, Optional


def optimizar_proteina(
    sequence_initial: str,
    parent1_seq: str,
    parent1_scores: List[float],
    parent2_seq: str,
    parent2_scores: List[float],
    objective_names: List[str],
    objective_description: str,
    mutable_positions: Optional[List[int]] = None,
    model_name: str = "llama3.1:8b",
    temperature: float = 0.3,
    ollama_host: str = "http://ollama:11434"
) -> str:

    if len(parent1_scores) != len(objective_names):
        raise ValueError("parent1_scores y objective_names deben tener el mismo tamaño")

    if len(parent2_scores) != len(objective_names):
        raise ValueError("parent2_scores y objective_names deben tener el mismo tamaño")

    client = ollama.Client(host=ollama_host)

    scores1 = ", ".join(f"{n}={v}" for n, v in zip(objective_names, parent1_scores))
    scores2 = ", ".join(f"{n}={v}" for n, v in zip(objective_names, parent2_scores))

    positions = f"Mutable positions: {mutable_positions}\n" if mutable_positions else ""

    system_prompt = (
        "You are an expert in protein engineering, directed evolution and "
        "multi-objective optimization. Generate a novel protein variant "
        "with improved fitness while preserving biochemical plausibility."
    )

    user_prompt = (
        f"Perform one round of directed evolution.\n"
        f"Objective: {objective_description}\n\n"
        f"Reference sequence:\n{sequence_initial}\n\n"
        f"{positions}"
        f"Parent 1:\n{parent1_seq}\nScores: {scores1}\n\n"
        f"Parent 2:\n{parent2_seq}\nScores: {scores2}\n\n"
        f"Generate ONE child sequence.\n"
        f"Use crossover, mutation or rational design.\n"
        f"Consider all objectives simultaneously.\n"
        f"Do not return either parent.\n"
        f"Keep the same sequence length.\n"
        f"Output ONLY the sequence."
    )

    response = client.chat(
        model=model_name,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ],
        options={"temperature": temperature}
    )

    return response["message"]["content"].strip()


#if __name__ == "__main__":
#
#    child = optimizar_proteina(
#        sequence_initial="MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE",
#        parent1_seq="QHVR",
#        parent1_scores=[0.81, -120.5, 0.42],
#        parent2_seq="RLIV",
#        parent2_scores=[0.73, -115.2, 0.61],
#        objective_names=["binding", "energy", "solubility"],
#        objective_description="maximize binding and solubility while minimizing energy",
#        mutable_positions=[39, 40, 41, 54]
#    )
#
#    print(child)