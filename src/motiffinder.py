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
