import pandas as pd
import json
import csv
import os
import re

def csvToTree(csv_data, output):
    df = pd.read_csv(csv_data)
    df['id'] = df['id'].astype(object)
    df['fitness'] = df['fitness'].astype(object)
    df['sequence'] = df['sequence'].astype(object)
    frst = df.iloc[0]
    df = df[df['id'] != 'Original']
    tree = {"name":str(frst['id']), "fitness":frst['fitness'].split(','), "sequence":frst['sequence']}
    nodes = {tree['name']:tree}

    for _, row in df.iterrows():
        name = row['id']
        father = row['father']
        fitness = row['fitness'].split(',')
        sequence = row['sequence']

        if not nodes.get(name):
            node_father = nodes.get(father)
            if node_father:
                node = {"name":name, "fitness":fitness, "sequence":sequence}
                if not node_father.get("children"): node_father['children'] = []
                node_father['children'].append(node)
                nodes[name] = node
    
    with open(output, 'w') as f:
        json.dump(tree, f, indent=2)

def hamming_distance(seq1, seq2):
    return sum(aa1 != aa2 for aa1, aa2 in zip(seq1, seq2))

# Custom function to save the population to a CSV file
def save_population_to_csv(population, generation, savefile_path):
    file_exists = os.path.isfile(savefile_path)

    with open(savefile_path, mode='a', newline='') as file:
        writer = csv.writer(file)
        if not file_exists:
            writer.writerow(['generation', 'id', 'father', 'pdb_file', 'fitness', 'sequence'])

        for ind in population:
            writer.writerow([
                generation,
                ind.id,
                ind.father,
                ind.pdb,
                ','.join(map(str, ind.fitness.values)),
                ''.join(map(str, ind))
            ])

def save_HallofFame(hof, savefile_path):
    file_exists = os.path.isfile(savefile_path)
    seen = set()

    with open(savefile_path, mode='a', newline='') as file:
        writer = csv.writer(file)
        if not file_exists:
            writer.writerow(['generation','id', 'father', 'pdb_file', 'fitness', 'sequence'])
        
        for ind in hof:
            if ind.id in seen: continue
            seen.add(ind.id)
            writer.writerow([
                ind.id.split('-')[0],
                ind.id,
                ind.father,
                ind.pdb,
                ','.join(map(str, ind.fitness.values)),
                ''.join(map(str, ind))
            ])

def add_sc(sc_file):
    df = pd.read_csv(sc_file, sep=r"\s+", skiprows=1)
    df.drop(columns=["SCORE:"], inplace=True)
    return df

def read_scfiles(path, output):
    dirs = [os.path.join(path, d) for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)) and d != 'tmp']
    df = pd.DataFrame()
    for dir in dirs:
        for file in os.listdir(dir):
            generation = int(re.search(r"g(\d+)", dir).group(1))
            if file.endswith(".sc"):
                df_ = add_sc(os.path.join(dir, file))
                df_['generation'] = generation
                df = pd.concat([df_,df], axis=0, ignore_index=True)
    
    df.drop(columns=["description"], inplace=True)
    df.to_csv(output, index=False)