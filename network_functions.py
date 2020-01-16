import networkx as nx
import pandas as pd
import numpy as np
import math

SMALL_VAL = 0.001

# node weight = perturbed_SI x (perturbed_SI/control_SI)
def get_activated_response_node_weight(perturbed_SI, control_SI):
	return perturbed_SI * (perturbed_SI/control_SI)

# node weight = control_SI x (control_SI/perturbed_SI)
def get_repressed_response_node_weight(perturbed_SI, control_SI):
	return control_SI * (control_SI/perturbed_SI)

# node weight = |log2(perturbed_SI) - log2(control_SI)|
def get_abs_FC_node_weight(perturbed_SI, control_SI):
	return abs(math.log2(perturbed_SI) - math.log2(control_SI))

# node weight = perturbed_SI/control_SI
def get_activated_FC_node_weight(perturbed_SI, control_SI):
        return perturbed_SI/control_SI

# node weight = control_SI/perturbed_SI
def get_repressed_FC_node_weight(perturbed_SI, control_SI):
	return control_SI/perturbed_SI

def get_edge_weight(n1, n2):
	temp = math.sqrt(n1 * n2)
	return 1/temp

def get_path_cost(G, path):
	cost = 0.0
	for i in range(len(path) - 1):
		cost += G[path[i]][path[i+1]]['weight']
	return cost

def get_edges_in_path(path, G):
	edges_in_path = set()
	for i in range(len(path) - 1):
		edges_in_path.add((path[i], path[i+1], G[path[i]][path[i+1]]['weight']))
	return edges_in_path

# SI is a pandas dataframe, indexed by gene, having only 1 column
# node weight = SI itself
def get_highest_activity_network(SI, nw):
	han = nx.DiGraph()
	for edge in nw.edges():
		if edge[0] in SI.index and edge[1] in SI.index: # got gene expression values
			if SI.loc[edge[0]][0] > 0 and SI.loc[edge[1]][0] > 0: # SI is positive for both nodes
				n1 = SI.loc[edge[0]][0]
				n2 = SI.loc[edge[1]][0]
				han.add_edge(edge[0], edge[1], weight = get_edge_weight(n1, n2))
	return han

# SI is given as a pandas dataframe, indexed by gene
# column 0 -> perturbed gene expression values
# column 1 -> control gene expression values
# node weight = perturbed_SI x (perturbed_SI/control_SI)
def get_activated_response_network(SI, nw):
	response_nw = nx.DiGraph()
	for edge in nw.edges():
		if edge[0] in SI.index and edge[1] in SI.index: # got gene expression values
			if SI.loc[edge[0]][0] > 0 and SI.loc[edge[0]][1] > 0: # SI is positive for node 1
				if SI.loc[edge[1]][0] > 0 and SI.loc[edge[1]][1] > 0: # SI is positive for node 2
					n1 = get_activated_response_node_weight(SI.loc[edge[0]][0], SI.loc[edge[0]][1])
					n2 = get_activated_response_node_weight(SI.loc[edge[1]][0], SI.loc[edge[1]][1])
					response_nw.add_edge(edge[0], edge[1], weight = get_edge_weight(n1, n2))
	return response_nw

# SI is given as a pandas dataframe, indexed by gene
# column 0 -> perturbed gene expression values
# column 1 -> control gene expression values
# node weight = control_SI x (control_SI/perturbed_SI)
def get_repressed_response_network(SI, nw):
	response_nw = nx.DiGraph()
	for edge in nw.edges():
		if edge[0] in SI.index and edge[1] in SI.index: # got gene expression values
			if SI.loc[edge[0]][0] > 0 and SI.loc[edge[0]][1] > 0: # SI is positive for node 1
				if SI.loc[edge[1]][0] > 0 and SI.loc[edge[1]][1] > 0: # SI is positive for node 2
					n1 = get_repressed_response_node_weight(SI.loc[edge[0]][0], SI.loc[edge[0]][1])
					n2 = get_repressed_response_node_weight(SI.loc[edge[1]][0], SI.loc[edge[1]][1])
					response_nw.add_edge(edge[0], edge[1], weight = get_edge_weight(n1, n2))
	return response_nw

# SI is given as a pandas dataframe, indexed by gene
# column 0 -> perturbed gene expression values
# column 1 -> control gene expression values
# node weight = perturbed_SI/control_SI
def get_activated_FC_network(SI, nw):
	fc_nw = nx.DiGraph()
	for edge in nw.edges():
		if edge[0] in SI.index and edge[1] in SI.index: # got gene expression values
			if SI.loc[edge[0]][0] > 0 and SI.loc[edge[0]][1] > 0: # SI is positive for node 1
				if SI.loc[edge[1]][0] > 0 and SI.loc[edge[1]][1] > 0: # SI is positive for node 2
					n1 = get_activated_FC_node_weight(SI.loc[edge[0]][0], SI.loc[edge[0]][1])
					n2 = get_activated_FC_node_weight(SI.loc[edge[1]][0], SI.loc[edge[1]][1])
					fc_nw.add_edge(edge[0], edge[1], weight = get_edge_weight(n1, n2))
	return fc_nw

# SI is given as a pandas dataframe, indexed by gene
# column 0 -> perturbed gene expression values
# column 1 -> control gene expression values
# node weight = control_SI/perturbed_SI
def get_repressed_FC_network(SI, nw):
	fc_nw = nx.DiGraph()
	for edge in nw.edges():
		if edge[0] in SI.index and edge[1] in SI.index: # got gene expression values
			if SI.loc[edge[0]][0] > 0 and SI.loc[edge[0]][1] > 0: # SI is positive for node 1
				if SI.loc[edge[1]][0] > 0 and SI.loc[edge[1]][1] > 0: # SI is positive for node 2
					n1 = get_repressed_FC_node_weight(SI.loc[edge[0]][0], SI.loc[edge[0]][1])
					n2 = get_repressed_FC_node_weight(SI.loc[edge[1]][0], SI.loc[edge[1]][1])
					fc_nw.add_edge(edge[0], edge[1], weight = get_edge_weight(n1, n2))
	return fc_nw

# SI is given as a pandas dataframe, indexed by gene
# column 0 -> perturbed gene expression values
# column 1 -> control gene expression values
# node weight = |log2(perturbed_SI) - log2(control_SI)|
# Since a subtraction is done here, and during edge weight calculation this will become the denominator,
# we need to avoid a value of 0. So if the node weight is about to evaluate to 0, we add SMALL_VAL to 
# one of the expressions
def get_abs_FC_network(SI, nw):
        abs_FC_nw = nx.DiGraph()
        for edge in nw.edges():
                if edge[0] in SI.index and edge[1] in SI.index: # got gene expression values
                        if SI.loc[edge[0]][0] > 0 and SI.loc[edge[0]][1] > 0: # SI is positive for node 1
                                if SI.loc[edge[1]][0] > 0 and SI.loc[edge[1]][1] > 0: # SI is positive for node 2
                                        n1_perturbed_SI = SI.loc[edge[0]][0]
                                        n1_control_SI = SI.loc[edge[0]][1]
                                        if n1_perturbed_SI == n1_control_SI:
                                                n1_control_SI += SMALL_VAL
                                        n1 = get_abs_FC_node_weight(n1_perturbed_SI, n1_control_SI)
                                        n2_perturbed_SI = SI.loc[edge[1]][0]
                                        n2_control_SI = SI.loc[edge[1]][1]
                                        if n2_perturbed_SI == n2_control_SI:
                                                n2_control_SI += SMALL_VAL
                                        n2 = get_abs_FC_node_weight(n2_perturbed_SI, n2_control_SI)
                                        if n1 == n2:
                                                n2 += SMALL_VAL
                                        abs_FC_nw.add_edge(edge[0], edge[1], weight = get_edge_weight(n1, n2))
        return abs_FC_nw

def get_sp_costs(G):
	sp_costs = nx.all_pairs_dijkstra_path_length(G)

	# Convert to a pandas dataframe
	sp_costs_df = pd.DataFrame.from_dict({i+'#'+j: sp_costs[i][j] for i in sp_costs.keys() for j in sp_costs[i].keys() if i!=j}, orient='index')

	return sp_costs_df


def get_all_sp_paths_costs(G):
	spaths_dict = nx.all_pairs_dijkstra_path(G)

	# Convert to a pandas dataframe
	# Index will be the full path, with nodes separated by #
	# The column value will be the cost of that path
	spaths_costs_df = pd.DataFrame.from_dict({'#'.join(spaths_dict[i][j]): get_path_cost(G, spaths_dict[i][j]) for i in spaths_dict.keys() for j in spaths_dict[i].keys() if i!=j}, orient='index')

	return spaths_costs_df

# Return all-pairs-shortest-path costs normalized by the number of hops in the path
def get_normalized_sp_costs(G):
	shortest_paths = nx.all_pairs_dijkstra_path(G)
	print("Done calculating shortest paths")

	# Normalize by number of hops
	norm_sp_costs = defaultdict(lambda: defaultdict(float))
	for src, dest_paths in shortest_paths.items():
		for dest, path in dest_paths.items():
			if len(path) > 1: # paths to self are ignored
				norm_sp_costs[src][dest] = get_path_cost(G, path)/(len(path)-1)
			else:
				norm_sp_costs[src][dest] = 0.0

	# Convert to a pandas dataframe
	norm_sp_costs_df = pd.DataFrame.from_dict({i+'#'+j: norm_sp_costs[i][j] for i in norm_sp_costs.keys() for j in norm_sp_costs[i].keys()}, orient='index')

	return norm_sp_costs_df

# Paths has the form ['a#b#c', 'b#c'], i.e it is a list of strings where each string is a path
# We need to return a dataframe indexed by the same list, with the cost of each path in the given graph
def get_costs_of_given_paths(G, paths):
        same_path_costs = pd.DataFrame(index = tuple(paths), columns = [0])
        for path in paths:
                same_path_costs.loc[path] = get_path_cost(G, path.split('#'))
        return same_path_costs
