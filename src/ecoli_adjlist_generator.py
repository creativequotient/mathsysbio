# TODO: clean up and document

import pandas as pd
import json
import networkx as nx
import motiffinder as mf
CSV_FILE = "ecoli_tn_genes.csv"
ADJLIST_JSON = "ecoli_ts_network.json"

def csv_to_adj_list(csv_file_name):
    START = 'TF'
    END = 'gene'
    DIRECTION = 'effect'
    CONFIDENCE = 'ev_level'

    data = pd.read_csv(csv_file_name)

    adj_list = {}

    for i in range(len(data)):
        starts = data.at[i, START].split(';') # multiple genes separated by ; for start
        end = data.at[i, END]
        is_plus = data.at[i, DIRECTION] == '+'
        confidence = data.at[i, CONFIDENCE]

        for start in starts:
            if start not in adj_list:
                adj_list[start] = []
            adj_list[start].append((end,is_plus,confidence))

    print(adj_list['nhaR'])
    print(len(adj_list.keys()))

    return adj_list

def adj_list_to_json(adj_list, filename):
    # save as JSON
    json.dump({'README': 'Adjacency list for E Coli transcriptional regulatory network'\
               'from https://www.pnas.org/content/114/38/10286 (Dataset S1).'\
               'All data is included, format of (neighbour) nodes is a tuple of'\
               '(gene : str, is_positive_regulation : bool, evidence_level)',
               'DATA': adj_list },
              open(filename,'w+')
              )

def adj_list_to_graph(adj_list, edge_predicate):
    graph = nx.DiGraph()
    for start, neighbours in adj_list.items():
        for edge in neighbours:
            if edge_predicate(edge):
                graph.add_edge(start,edge[0], is_positive=edge[1], confidence=edge[2])
    return graph

def open_adj_list_json(filename):
    with open(filename) as f:
        adj_list = json.load(f)['DATA']
    return adj_list

def default_adj_list_to_graph(adjlist):
    return adj_list_to_graph(adj_list,  lambda e : e[2] == 2 or e[2] == float('inf'))

import matplotlib.pyplot as plt

if __name__ == '__main__':
    adj_list = open_adj_list_json(ADJLIST_JSON)
    graph = default_adj_list_to_graph(adj_list)
    adj_csv = pd.read_csv("ecoli_tn_genes.csv")
    transcription_factors = adj_csv[adj_csv.ev_level > 1].TF.unique()
    SIMS = mf.find_SIMS(graph, transcription_factors)
    SIMS = list(filter(lambda x: len(x) > 1, SIMS))
    SIMS = list(sorted(SIMS, key = lambda x: len(x)))
    nx.draw(nx.subgraph(graph, SIMS[-1]), with_labels=True)
    plt.show()
    # print(graph.edges['cra','pdhR']) # verify infinity is present
    # print(graph.edges['sgrR','sroA']) # verify 2.0 is present
    # subgraph = graph.subgraph(['cra','pdhR','sgR','sroA','thiP'])
    # nx.draw(subgraph,with_labels=True,edge_color='gray')
    # plt.show()
