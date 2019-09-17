import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import multiprocessing as mp
import functools

def get_random_nodes(indegrees):
    sum_degrees = np.sum(indegrees)
    probability_dist = indegrees / sum_degrees
    output = np.random.choice(np.arange(len(indegrees)), (2), True, probability_dist)
    return output

def generate_scale_free_graph(nodes, edges):

    # Generate BA-graph with n = nodes, m = 1
    graph = nx.generators.random_graphs.barabasi_albert_graph(nodes, 1)
    edge_list = graph.edges()

    # Compute in-degrees
    indegrees = np.zeros(nodes)

    for i, j in edge_list:
        indegrees[i] += 1
        indegrees[j] += 1

    additional_edges = edges - nodes + 1

    for _ in range(additional_edges):
        while True:
            i, j = get_random_nodes(indegrees)
            edge_list = graph.edges()
            if (i, j) in edge_list or (j, i) in edge_list:
                continue
            else:
                indegrees[i] += 1
                indegrees[j] += 1
                graph.add_edge(i, j)
                break

    return graph

def generate_scale_free_graphs(edges, nodes, n, data_path):
    try:
        path = os.path.join(data_path, "{0}_{1}".format(nodes, edges))
        os.mkdir(path)
    except Exception as e:
        pass

    frequency = []

    for i in range(n):
        graph = generate_scale_free_graph(nodes, edges)
        output_path = os.path.join(path, "{0}_{1}_{2}.edgelist".format(nodes, edges, i))
        nx.write_edgelist(graph, output_path)
        frequency.append(len(list(graph.selfloop_edges())))

    stats_path = os.path.join(data_path, "stats", "{0}_{1}_frequency.csv".format(nodes, edges))

    freq_arr = np.array(frequency, dtype=np.int32)

    np.savetxt(stats_path, freq_arr, delimiter=",")


if __name__ == "__main__":
    nodes = 420
    edges = 520
    p = mp.Pool(processes=8)
    edges = [520 + i * 50 for i in range(1000)]
    partial = functools.partial(generate_scale_free_graphs, nodes=nodes, n=1000, data_path="data")
    p.map(partial, edges)
