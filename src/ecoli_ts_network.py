"""ECOLI TRANSCRIPTION NETWORK
This module contains functions to create a graph representing an E Coli transcription
network. The network is from this paper: https://www.pnas.org/content/114/38/10286 
(Dataset S1, which is an xlsx file. The data from the "gene" sheet is used).

The representation of the transcription network is fairly simplistic: for every TF-gene
pair, only the positive/negative nature of the regulation and level of evidence
(1, 2 or inf) is stored. For genes regulated by more than 1 TF, it is also assumed
that there are no interactions between the TFs (ie. no logic).

FUNCTIONS AVAILABLE:
* open_graph (most useful)
* graph_high_confidence (most useful)
* open_adjlist
* csv_to_adjlist
* save_adjlist

GRAPH FORMAT:
Each node represents both the gene and protein product, which can be a TF.
Each edge represents a TF-gene interaction. The starting node is the TF and
the destination node is the gene being regulated. Every interaction also has
the following attributes:
    1. is positive regulation (bool)
    2. evidence level (1, 2 or inf)

The graph can be created in 2 formats: an adjacency list (dict), or NetworkX.DiGraph.
The DiGraph format is probably better as it's easier to extend to work with other
networks and allows us to use the NetworkX graph functions.

Format of NetworkX.DiGraph:
    The additional attributes of each interaction are stored as edge attributes:
    1. is_positive (bool): is postive regulation
    2. evidence (float - 1, 2 or inf): evidence level

Format of adjacency list (dict):
    Key: name of TF
    Value: tuple of interactions between this TF and a gene. Each interaction is
           a tuple of (gene name (str), is positive regulation (bool), evidence
           level (1, 2 or inf)).

NON-PYTHON DEPENDENCIES:
1. ecoli_ts_network.csv : csv_to_adjlist
    This is from the gene sheet of the xlsx file from the paper (see top of doc).
2. ecoli_ts_network.json : open_adjlist and open_graph
    This file can be recreated using csv_to_adjlist and using save_adjlist on
    the resulting adjacency list.

NOTE:
There are some entries in the CSV/excel file with multiple TFs (eg. inhrA;inhrB).
I've assumed that this means both TFs individually regulate the gene. This just my
own guess - I couldn't find an explanation for it in the paper.
"""

import pandas as pd
import json
import networkx as nx
import motiffinder as mf

CSV_FILE = "ecoli_ts_network.csv"
ADJLIST_JSON = "ecoli_ts_network.json"

def csv_to_adjlist(csv_filename = CSV_FILE):
    """Creates an adjacency list of the E Coli transcription network from the CSV file
    from the paper

    NOTE:
    There are some entries in the CSV file with multiple TFs (eg. inhrA;inhrB).
    I've assumed that this means both TFs individually regulate the gene. This just my
    own guess - I couldn't find an explanation for it in the paper.

    Args:
        csv_filename (str) : name of CSV file to read from

    Returns:
        A dict which is the adjacency list
    """

    START = 'TF'
    END = 'gene'
    DIRECTION = 'effect'
    EVIDENCE = 'ev_level'

    data = pd.read_csv(csv_filename)

    adj_list = {}

    for i in range(len(data)):
        starts = data.at[i, START].split(';') # multiple genes separated by ; for start
        end = data.at[i, END]
        is_plus = data.at[i, DIRECTION] == '+'
        evidence = data.at[i, EVIDENCE]

        for start in starts:
            if start not in adj_list:
                adj_list[start] = []
            adj_list[start].append((end, is_plus, evidence))

    return adj_list

def save_adjlist(adjlist, filename = ADJLIST_JSON):
    """Saves the given adjacency list to a JSON file.

    Args:
        adjlist (dict): adjacency list
        filename (str, optional) : name of JSON file to save to. Default is the value of
        ADJLIST_JSON (ecoli_ts_network.json currently).
    """

    json.dump(adjlist, open(filename,'w+'))


def open_adjlist(filename = ADJLIST_JSON):
    with open(filename) as f:
        adjlist = json.load(f)
    return adjlist

def adjlist_to_graph(adjlist, predicate = lambda tf, gene : True):
    """Creates a NetworkX.DiGraph from the given the adjacency list (dict).

    Args:
        adjlist (dict): adjacency list
        predicate (function, optional): predicate for whether to include a particular
            TF-gene interaction in the graph. Should be a function that takes which 
            takes in:
            1. transcription factor (str)
            2. gene (gene name, is_positive_regulation and confidence_level)
            and returns True for pairs to include in the graph.
            (default: function that always returns True).
    
    Returns:
        A NetworkX.DiGraph
    """

    graph = nx.DiGraph()
    for start, neighbours in adjlist.items():
        for edge in neighbours:
            if predicate(start, edge):
                graph.add_edge(start,edge[0], is_positive=edge[1], evidence=edge[2])
    return graph

def adjlist_to_graph_high_confidence(adjlist):
    """Creates a NetworkX.DiGraph which includes only high confidence TF-gene pairs (2 or inf)."""
    return adjlist_to_graph(adjlist,  lambda tf, gene : gene[2] >= 2 )

def open_graph(filename = ADJLIST_JSON):
    """Returns a NetworkX.DiGraph of the E Coli transcription network.

    Args:
        filename (str) : name of JSON file containing the adjacency list of the graph. 
            See the module documentation for the format of the adjacency list and how 
            to create it.

    Returns:
        A NetworkX.DiGraph
    """

    return adjlist_to_graph(open_adjlist(filename))

def high_evidence_graph(graph):
    """Returns a new NetworkX.DiGraph containing edges with evidence level of 2 or inf
    
    Args:
        graph (NetworkX.DiGraph)

    Returns:
        A new NetworkX.DiGraph
    """

    new_graph = nx.DiGraph()
    for edge in graph.edges:
        if graph.get_edge_data(*edge)['evidence'] >= 2:
            new_graph.add_edge(*edge)
    return new_graph

def get_transcription_factors(graph):
    """Return a list of transcription factors in the graph.
    """

    return list(filter(lambda node : tuple(graph.successors(node)), graph.nodes))

# Testing code
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # adjlist = csv_to_adjlist()
    # save_adjlist(adjlist)
    graph = open_graph()
    hi_conf = graph_high_confidence(graph)
    print(graph.size(), hi_conf.size())
    print(len(graph), len(hi_conf))

    # adj_list = open_adj_list_json(ADJLIST_JSON)
    # graph = default_adj_list_to_graph(adj_list)
    # adj_csv = pd.read_csv("ecoli_tn_genes.csv")
    # transcription_factors = adj_csv[adj_csv.ev_level > 1].tf.unique()
    # SIMS = mf.find_SIMS(graph, transcription_factors)
    # SIMS = list(filter(lambda x: len(x) > 1, SIMS))
    # SIMS = list(sorted(SIMS, key = lambda x: len(x)))
    # nx.draw(nx.subgraph(graph, SIMS[-1]), with_labels=True)
    # plt.show()
    # print(graph.edges['cra','pdhR']) # verify infinity is present
    # print(graph.edges['sgrR','sroA']) # verify 2.0 is present
    # subgraph = graph.subgraph(['cra','pdhR','sgR','sroA','thiP'])
    # nx.draw(subgraph,with_labels=True,edge_color='gray')
    # plt.show()
