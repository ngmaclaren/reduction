import numpy as np
import networkx as nx
from networkx.generators.community import LFR_benchmark_graph
# import matplotlib.pyplot as plt

seed = 12345

N = 1000 # number of nodes
p = 0.05 # probability of an edge, Erdős-Rényi graph
m = 2 # number of edges to add, Barabási-Albert graph
phk = 0.1 # target clustering coefficient for Holme-Kim graph

# Erdős-Rényi graph
er = nx.erdos_renyi_graph(N, p, seed = seed)

# Barabási-Albert graph
ba = nx.barabasi_albert_graph(N, m, seed = seed)

# Holme-Kim graph
hk = nx.powerlaw_cluster_graph(N, m, phk, seed = seed)

# Lancichinetti-Fortunato-Radicchi graph, using close to the NetworkX example settings
lfr_orig = LFR_benchmark_graph(N, 3, 1.5, .1, average_degree = 4, min_community = 20, seed = seed)
lfr_orig.remove_edges_from(nx.selfloop_edges(lfr_orig))
lcc = max(nx.connected_components(lfr_orig), key = len)
lfr = lfr_orig.subgraph(lcc).copy()
# nx.draw(lfr, node_size = 10, with_labels = False)
# plt.show()

# Write graphs as edgelists for compatability
nx.write_edgelist(er, '../data/er.txt', data = False)
nx.write_edgelist(ba, '../data/ba.txt', data = False)
nx.write_edgelist(hk, '../data/hk.txt', data = False)
nx.write_edgelist(lfr, '../data/lfr.txt', data = False)
