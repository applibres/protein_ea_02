import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def hierarchy_pos(G, root=None, width=1.0, vert_gap=0.2, vert_loc=0, xcenter=0.5):
    if not nx.is_tree(G):
        return nx.spring_layout(G)
    if root is None:
        root = [n for n, d in G.in_degree() if d == 0][0]
    def _hierarchy_pos(G, root, width=1.0, vert_gap=0.2, vert_loc=0,
                       xcenter=0.5, pos=None, parent=None):
        if pos is None:
            pos = {root: (xcenter, vert_loc)}
        else:
            pos[root] = (xcenter, vert_loc)
        children = list(G.successors(root))
        if len(children) != 0:
            dx = width / len(children)
            nextx = xcenter - width / 2 - dx / 2
            for child in children:
                nextx += dx
                pos = _hierarchy_pos(G, child, width=dx, vert_gap=vert_gap,
                                     vert_loc=vert_loc - vert_gap, xcenter=nextx, pos=pos,
                                     parent=root)
        return pos
    return _hierarchy_pos(G, root, width, vert_gap, vert_loc, xcenter)


def graphTree(csv_data, output):
    df = pd.read_csv(csv_data)
    
    # Create graph
    G = nx.DiGraph()
    
    # Add nodes and edges
    for _, row in df.iterrows():
        node_id = row['id']
        father_id = row['father']
        G.add_node(node_id, label=node_id)
        if father_id not in ["Original", "", None]:
            G.add_edge(father_id, node_id)
    
    pos = hierarchy_pos(G)
    
    # Draw graph
    plt.figure(figsize=(14, 10))
    nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=2000, arrows=True, font_size=10)
    plt.title("Individual and their muations Tree")
    plt.savefig(f"{output}/generations_tree.png", format="png", dpi=300, bbox_inches='tight')

