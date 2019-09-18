# TODO: fix a bug somewhere, my "DORs" aren't connected.
# But i'm pretty sure my algo is correct

import matplotlib.pyplot as plt
import networkx as nx
from ecoli_adjlist_generator import *

def find_dors(graph):
    # the bottom row (genes being regulated) in a DOR have an identical set of predecessors and cannot have self loop
    
    # get list of predecessors for every node, excluding nodes with a self loop.
    # the keys of the dict are the node, the value is a list of predecessors
    nodes = set(map(lambda e: e[0], graph.edges)) # nodes with incoming edges (ie. with predecessors)
    regulon_regulator_map = {}
    for node in nodes:
        potential_regulon = tuple(sorted(graph.predecessors(node)))
        
        # skip if there is a self loop or only 1 gene in regulon
        if node in potential_regulon or len(potential_regulon)==1:
            continue 
        if potential_regulon not in regulon_regulator_map:
            regulon_regulator_map[potential_regulon] = []
        regulon_regulator_map[potential_regulon].append(node) # node is the regulator
    
    DOR_list = []
    for regulons, regulators in regulon_regulator_map.items():
        # skip if there is only 1 regulator (then it's a SIM, not a DOR)
        if len(regulators) == 1:
            continue

        # check that the regulators do not regulate each other
        # self loops are fine
        # ie. every regulator is not contained in any other regulator's predecessor/neighbour list.
        for reg in regulators:
            for successor in graph.successors(reg):
                if successor != reg and successor in regulators:
                    # regulators are regulating each other, so not a DOR
                    continue
        DOR_list.append((tuple(regulons),tuple(regulators)))

    return DOR_list

if __name__ == '__main__':
    adj_list = open_adj_list_json(ADJLIST_JSON)
    graph = default_adj_list_to_graph(adj_list)
    all_dors = find_dors(graph)
    print(all_dors)
    print(len(all_dors))
    if len(all_dors) > 0:
        dor = all_dors[3]
        all_nodes = dor[0] + dor[1]
        dor_subgraph = graph.subgraph(all_nodes)
        nx.draw(dor_subgraph,with_labels=True,edge_color='gray')
        plt.show()
