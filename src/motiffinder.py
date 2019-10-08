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
    """Returns all the SIMs in the graph, which includes "SIMs" with 1 gene.

    Args:
        graph (NetworkX.DiGraph)
        transcription_factors : list of names of TFs in the graph (ie. nodes with out 
            degree >= 1)
        group (bool) : whether to return a list of SIMs, where each SIM is a flattened 
            list of the genes followed by the TF (True), or a dict indexed by TF (False)

    Returns:
        The SIMs in this graph (including "SIMs" with 1 gene) as either a list of dict
        (see the group arg)
    """

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

def find_DORs(graph, group = False):
    """Returns a list of all the complete (fully connected) DORs in the graph.

    Basis of algorithm:
    1. the bottom row (genes being regulated) in a DOR have an identical set 
       of predecessors and cannot have self loop
    2. get list of predecessors for every node, excluding nodes with a self loop.

    Args:
        graph (NetworkX.DiGraph)
        group (bool) : True to return each DOR as a flat tuple of regulators 
            followed by regulons False to return each DOR as a nested tuple of 
            (regulators, regulons)

    Returns:
        a list of fully connected DORs in the graph (see group argument for 
        details on the format of each DOR), or an empty list if there are none
    """

    nodes = set(map(lambda e: e[0], graph.edges)) # nodes with incoming edges
    regulator_regulon_map = {} # regulators are predecessors
    for node in nodes:
        regulators = tuple(sorted(graph.predecessors(node)))
        
        # and skip nodes with only 1 regulator (then it's a potential SIM) or no regulators
        if len(regulators) <= 1:
            continue

        # don't allow self loops
        if node in regulators:
            if regulators in regulator_regulon_map:
                del regulator_regulon_map[regulators]
            continue

        if regulators not in regulator_regulon_map:
            regulator_regulon_map[regulators] = []
        regulator_regulon_map[regulators].append(node)
   
    DORs = []

    for regulators, regulon in regulator_regulon_map.items():
        regulon = tuple(regulon)
        exclude = False
        
        if len(regulon) <= 1:
            continue
        # check regulators are not regulating each other
        regulators_subgraph = graph.subgraph(regulators)
        
        for edge in regulators_subgraph.edges:
            if edge[0] != edge[1]:
                exclude = True
                break
        if exclude:
            continue
        # check regulons are not regulating regulators also
        for gene in regulon:
            for reg in regulators:
                if graph.has_edge(gene,reg):
                    exclude = True
                    break
            if exclude:
                break
        if exclude:
            continue
        
        if group:
            DORs.append(regulators + regulon)
        else:
            DORs.append((regulators,regulon))

    return DORs

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
