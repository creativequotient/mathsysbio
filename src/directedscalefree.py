import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import multiprocessing as mp
import functools
import random
from motiffinder import find_ffl
from undirectedscalefree import generate_scale_free_graph

'''
Get indegree and outdegree distribution
Argument(s):
    graph: networkx Graph or DiGraph
Return(s):
    dictionary of indegrees and outdegrees
'''
def get_degree_distributions(graph):
    nodes = graph.number_of_nodes()
    indegrees = [0 for i in range(nodes)]
    outdegrees = [0 for i in range(nodes)]
    for outgoing_node, incoming_node in graph.edges():
        outdegrees[int(outgoing_node)] += 1
        indegrees[int(incoming_node)] += 1
    return {"indegrees":indegrees,
            "outdegrees":outdegrees}


'''
Get random node(s) based on degree distribution
Argument(s):
    degrees: list of degrees with element at index i corresponding to degree of node i
    n: number of nodes to select
    replace: select nodes with/without replacement
Return(s):
    node(s) randomly selected based on degree distribution
'''
def get_random_nodes(degrees, n, replace=True):
    sum_degrees = np.sum(degrees)
    probability_dist = degrees / sum_degrees
    output = np.random.choice(np.arange(len(degrees)), (n), replace, probability_dist)
    return output


'''
Generate random gene network consisting of n nodes and (n - 1) edges following
scale free outdegree distribution and poisson indegree distribution
Argument(s):
  nodes: number of nodes within the the network
Return(s):
  networkx Graph object of n nodes and (n-1) directed edges
'''
def minimal_gene_network(nodes):

    # Initial setup of 2 node graph with directed edge from node 0 to node 1
    graph = nx.DiGraph()
    graph.add_edge(0, 1)
    outdegree_distribution = [5,5]

    for node in range(2, nodes):
        #print("Creating node {}".format(node))
        #outdegree_distribution += [random.randint(1,node**3)]
        outdegree_distribution += [1]
        while True:
            outgoing_node = get_random_nodes(outdegree_distribution[:-1], 1)[0]
            incoming_node = node

            #incoming_node, outgoing_node = outgoing_node, incoming_node

            #incoming_node = random.randint(0, node - 1)
            # if random.randint(0,2) == 0:
            #     outgoing_node = get_random_nodes(outdegree_distribution[:-1], 1)[0]
            #     incoming_node = node
            # else:
            #     outgoing_node = node
            #     incoming_node = random.randint(0, node - 2)

            if (outgoing_node, incoming_node) in graph.edges():
                continue
            else:
                graph.add_edge(outgoing_node, incoming_node)
                outdegree_distribution[outgoing_node] += 1
                break

    return graph


def simulate_gene_network2(nodes, edges):
    # graph = minimal_gene_network(nodes)
    ba = nx.generators.barabasi_albert_graph(nodes, 1)
    graph = nx.DiGraph()
    for i, j in ba.edges():
        if random.randint(0,1) == 0:
            graph.add_edge(i,j)
        else:
            graph.add_edge(j,i)

    outdegree_distribution = get_degree_distributions(graph)["outdegrees"]

    for edge in range(edges - nodes + 1):
        while True:
            outgoing_node = get_random_nodes(outdegree_distribution, 1)[0]
            incoming_node = random.randint(0, nodes - 1)
            if (outgoing_node, incoming_node) in graph.edges():
                continue
            else:
                graph.add_edge(outgoing_node, incoming_node)
                outdegree_distribution[outgoing_node] += 1
                break

    return graph

def simulate_gene_network(nodes, edges):
    # graph = minimal_gene_network(nodes)
    ba = generate_scale_free_graph(nodes,edges)
    graph = nx.DiGraph()
    for i, j in ba.edges():
        if random.randint(0,1) == 0:
            graph.add_edge(i,j)
        else:
            graph.add_edge(j,i)
    return graph

def generate_directed_scale_free_graphs(edges, nodes, n, data_path):
    try:
        path = os.path.join(data_path, "{0}_{1}".format(nodes, edges))
        os.mkdir(path)
    except Exception as e:
        pass

    frequency = []

    for i in range(n):
        graph = simulate_gene_network(nodes, edges)
        output_path = os.path.join(path, "{0}_{1}_{2}.edgelist".format(nodes, edges, i))
        nx.write_edgelist(graph, output_path)
        frequency.append(find_ffl(graph))

    stats_path = os.path.join("..", "stats", "ba_graphs", "feedforwardloops", "{0}_{1}_frequency.csv".format(nodes, edges))

    freq_arr = np.array(frequency, dtype=np.int32)

    np.savetxt(stats_path, freq_arr, delimiter=",")


if __name__ == "__main__":
    nodes = 420
    edges = 420
    p = mp.Pool(processes=16)
    edges = [edges + i * 10 for i in range(500)]
    partial = functools.partial(generate_directed_scale_free_graphs, nodes=nodes, n=2000, data_path="../graphs/ba_graphs")
    p.map(partial, edges)

    # graph = simulate_gene_network(200, 10000)
    # graph_in_degree = [x[1] for x in graph.in_degree()]
    # graph_out_degree = [x[1] + 1 for x in graph.out_degree()]
    # max_out_degree = max(graph_out_degree)
    # print(graph_out_degree)
    # out_degrees = np.ones(max_out_degree + 1)
    # for i in graph_out_degree:
    #     out_degrees[i] += 1
    # print(out_degrees)
    # print(graph.edges())
    # degrees = get_degree_distributions(graph)
    # indegrees = degrees["indegrees"]
    # outdegrees = degrees["outdegrees"]
    # fig, (ax1, ax2) = plt.subplots(2)
    # ax1.scatter(np.arange(len(graph_in_degree)), graph_in_degree)
    # ax1.hist(graph_in_degree)
    # ax2.scatter(np.arange(2, len(out_degrees) + 1), out_degrees[1:])
    # nx.draw(graph)
    # ax2.loglog()
    # plt.show()
    # nx.draw(graph)
    # plt.show()
