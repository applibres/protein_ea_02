import pandas as pd
import json

def csvToTree(csv_data, output):
    df = pd.read_csv(csv_data)
    df['id'] = df['id'].astype(object)
    df['fitness'] = df['fitness'].astype(float)
    frst = df.iloc[0]
    df = df[df['id'] != 'Original']
    tree = {"name":str(frst['id']), "fitness":frst['fitness']}
    nodes = {tree['name']:tree}

    for _, row in df.iterrows():
        name = row['id']
        father = row['father']
        fitness = row['fitness']

        if not nodes.get(name):
            node_father = nodes.get(father)
            if node_father:
                node = {"name":name, "fitness":fitness}
                if not node_father.get("children"): node_father['children'] = []
                node_father['children'].append(node)
                nodes[name] = node
    
    with open(output, 'w') as f:
        json.dump(tree, f, indent=2)

