import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def find_ffl(graph):
    adj_matrix = nx.to_numpy_matrix(graph)

    # Preprocess adjacency_matrix
    rows, cols = adj_matrix.shape

    for i in range(rows):
        if adj_matrix[i,i] == 1:
            adj_matrix[i,:] = 0
            adj_matrix[:,i] = 0

    adj_matrix = adj_matrix - adj_matrix.T
    adj_matrix[adj_matrix == -1] = 0

    counter = 0

    for i in range(rows):
        for j in range(cols):
            if adj_matrix[i,j] == 1:
                for k in range(rows):
                    if adj_matrix[k,j] == 1 and adj_matrix[i,k] == 1:
                        counter += 1

    return counter

def find_SIMS(graph, transcription_factors, group=True):
    # Get all nodes from graph
    nodes = graph.nodes()

    # Get index of transcription factors
    # transcription_factor_indices = []
    # for tf in transcription_factors:
    #     transcription_factor_indices.append(nodes.index(tf))

    # Get candidate nodes
    candidate_nodes = []
    for node, degree in graph.in_degree():
        if degree == 1:
            candidate_nodes.append(node)

    SIMS = {}
    # Sort candidates by transcription transcription_factors
    for tf in transcription_factors:
        tf_targets = []
        for candidate in candidate_nodes:
            if graph.has_edge(tf, candidate):
                tf_targets.append(candidate)
        SIMS[tf] = tf_targets

    if not group:
        return SIMS

    else:
        SIM_modules = []
        for tf in list(SIMS.keys()):
            module = SIMS[tf] + [tf,]
            SIM_modules.append(module)
        return SIM_modules

# Returns INDEX of terminal nodes
def find_terminal_nodes(graph, indegree=1):
    # Convert graph to adjacency matrix
    adj_matrix = nx.to_numpy_matrix(graph)
    rows, cols = adj_matrix.shape


    # Locate nodes with an outdegree of 0
    candidate_nodes = []
    for row in range(rows):
        #print(adj_matrix[row])
        if np.sum(adj_matrix[row]) == 0:
            candidate_nodes.append(row)

    terminal_nodes = []
    # Locate nodes with an indegree of ${indegree}
    for node in candidate_nodes:
        if np.sum(adj_matrix[:,node]) == indegree:
            terminal_nodes.append(node)

    return terminal_nodes

if __name__ == "__main__":
    graph = nx.DiGraph()
    graph.add_edge("TF1","A")
    graph.add_edge("TF1","B")
    graph.add_edge("TF1","C")
    graph.add_edge("TF2","A")
    graph.add_edge("TF2","D")
    SIM_modules = find_SIMS(graph, ["TF1","TF2"])
    #TFs = list(SIM_modules.keys())
    #sub_graph = SIM_modules[TFs[0]] + [TFs[0],]
    #nx.draw(nx.subgraph(graph,SIM_modules[1]), with_labels=True)
    nx.draw(graph, with_labels=True)
    plt.show()
