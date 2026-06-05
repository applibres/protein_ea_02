import pandas as pd
import json
import csv
import os
import re


def csvToTree(csv_data, output):
    df = pd.read_csv(csv_data)
    df['id']       = df['id'].astype(object)
    df['fitness']  = df['fitness'].astype(object)
    df['sequence'] = df['sequence'].astype(object)
    df['nmut']     = df['nmut'].astype(object)

    frst = df.iloc[0]
    df   = df[df['id'] != 'Original']

    tree  = {"name": str(frst['id']), "fitness": frst['fitness'],
             "sequence": frst['sequence'], "nmut": frst['nmut']}
    nodes = {tree['name']: tree}

    for _, row in df.iterrows():
        name   = row['id']
        father = row['father']
        if nodes.get(name):
            continue
        node_father = nodes.get(father)
        if node_father:
            node = {"name": name, "fitness": row['fitness'],
                    "sequence": row['sequence'], "nmut": row['nmut']}
            node_father.setdefault('children', []).append(node)
            nodes[name] = node

    with open(output, 'w') as f:
        json.dump(tree, f, indent=2)


def hamming_distance(seq1, seq2):
    return sum(a != b for a, b in zip(seq1, seq2))


def save_population_to_csv(population, generation, savefile_path):
    """Save a list of pymoo Individuals to CSV (append mode)."""
    file_exists = os.path.isfile(savefile_path)

    with open(savefile_path, mode='a', newline='') as file:
        writer = csv.writer(file)
        if not file_exists:
            writer.writerow(['generation', 'id', 'father', 'nmut',
                             'pdb_file', 'fitness', 'sequence'])

        for ind in population:
            writer.writerow([
                generation,
                ind.id,
                ind.father,
                ind.nmut,
                ind.pdb,
                ','.join(map(str, ind.fitness)),   # ind.F  (np.ndarray)
                ind.sequence(),              # ind.sequence()
            ])


def save_HallofFame(hof, savefile_path):
    """Save Hall-of-Fame individuals to CSV."""
    file_exists = os.path.isfile(savefile_path)
    seen = set()

    with open(savefile_path, mode='a', newline='') as file:
        writer = csv.writer(file)
        if not file_exists:
            writer.writerow(['generation', 'id', 'father', 'nmut',
                             'pdb_file', 'fitness', 'sequence'])

        for ind in hof:
            if ind.id in seen:
                continue
            seen.add(ind.id)
            writer.writerow([
                ind.id.split('-')[0],
                ind.id,
                ind.father,
                ind.nmut,
                ind.pdb,
                ','.join(map(str, ind.fitness)),   # ind.F  (np.ndarray)
                ind.sequence(),              # ind.sequence()
            ])


def get_sequence(generation, indiv, src_path):
    file = os.path.join(src_path, f"g{generation}/pop_g{generation}_AA.txt")
    with open(file, "r") as f:
        lines = f.readlines()

    line     = lines[indiv].strip()
    sequence = re.findall(r"[A-Za-z]", line.split(']')[0])
    return sequence


def add_sc(sc_file):
    df = pd.read_csv(sc_file, sep=r"\s+", skiprows=1)
    df.drop(columns=["SCORE:"], inplace=True)
    return df


def read_scfiles(path, output):
    original_sequence = get_sequence(0, 0, path)
    dirs = [
        os.path.join(path, d) for d in os.listdir(path)
        if os.path.isdir(os.path.join(path, d)) and d != 'tmp'
    ]
    df = pd.DataFrame()
    for dir in dirs:
        generation = int(re.search(r"g(\d+)", dir).group(1))
        for file in os.listdir(dir):
            if file.endswith(".sc"):
                n_indiv  = int(re.search(r"_(\d+)", file).group(1))
                sequence = get_sequence(generation, n_indiv, path)
                df_      = add_sc(os.path.join(dir, file))
                df_['generation'] = generation
                df_['sequence']   = "".join(sequence)
                df_['nmut']       = hamming_distance(original_sequence, sequence)
                df = pd.concat([df_, df], axis=0, ignore_index=True)

    df.drop(columns=["description"], inplace=True)
    df.to_csv(output, index=False)


def save_evolutions_statistics(logbook, file, nobj):
    """Save per-generation stats from a plain list-of-dicts logbook."""
    with open(file, "w") as stat_file:
        header = ["gen"] + [
            f"{stat}_obj{i}"
            for stat in ["avg", "min", "max", "std"]
            for i in range(nobj)
        ]
        stat_file.write(",".join(header) + "\n")

        for entry in logbook:
            row = [str(entry["gen"])]
            for stat in ["avg", "min", "max", "std"]:
                row += [str(v) for v in entry[stat]]
            stat_file.write(",".join(row) + "\n")


def save_population_aa(gendir, gen, pop, offspring_output_pdbfiles):
    """Save AA representation of the population to a text file."""
    path = os.path.join(gendir, f"pop_g{gen}_AA.txt")
    with open(path, "w") as output_file:
        for ind, pdb in zip(pop, offspring_output_pdbfiles):
            output_file.write(f"{ind.sequence()} {ind.fitness} {pdb}\n")


def num2str(num, n_digits=2):
    return f"{num:0>{n_digits}}"