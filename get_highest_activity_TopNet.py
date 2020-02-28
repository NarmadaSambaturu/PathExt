# Code to create the highest activity TopNet (HA TopNet) given inputs PPI and gene expression data,

# NOTE: This code assumes the microarray data has
#	Column 0 -> Gene symbols, 
# 	Columns 1 through m -> normalised SI
#	Header for the sample of interest is given as a user input

import pandas as pd
import random
import networkx as nx
import numpy as np
import sys
import math

import microarray_functions as mic_fun
import network_functions as net_fun
import percentile_functions as perc_fun

def get_edges_in_path(path, G):
	edges_in_path = set()
	for i in range(len(path) - 1):
		edges_in_path.add((path[i], path[i+1], G[path[i]][path[i+1]]['weight']))
	return edges_in_path

# SI is given as a pandas dataframe, indexed by gene
# column 0 -> disease gene expression values
def combine_data_get_sp_paths_costs_ha(G_unweighted, SI, ha_nw_fname):
	# Map gene expression values onto unweighted network
	G_ha = net_fun.get_highest_activity_network(SI, G_unweighted)
	print("Got highest activity base network with ", len(G_ha.nodes()), " nodes and ", len(G_ha.edges()), " edges")

	# Drop SI values for genes which don't map to response network
	genes_to_drop = set(SI.index) - set(G_ha.nodes())
	SI = SI.drop(genes_to_drop)

	# Write the highest activity network to file
	nx.write_weighted_edgelist(G_ha, ha_nw_fname, delimiter = '\t')
	# Get all-pairs-shortest-path costs
	# Return value is a pandas dataframe, indexed by the string src#dest
	Pij = net_fun.get_all_sp_paths_costs(G_ha)
	print("Got shortest path costs for ", Pij.shape[0], " node-pairs")

	return Pij, G_ha


if len(sys.argv) != 8:
	print("argv[1] = microarray data file (tab-delimited, with header)")
	print("argv[2] = name of sample to study")
	print("argv[3] = unweighted (directed) network file")
	print("argv[4] = percentile threshold")
	print("argv[5] = path length threshold")
	print("argv[6] = output file for highest activity base network")
	print("argv[7] = output file for HA TopNet")
	sys.exit(1)

# Set inputs
data_fname = sys.argv[1]
sample_of_interest = sys.argv[2] # This is the sample we want to study
unweighted_nw_fname = sys.argv[3]
percentile = float(sys.argv[4]) # We'll only keep paths whose cost < this threshold
path_length_thresh = int(sys.argv[5]) # We'll only keep paths with length >= this threshold
ha_nw_fname = sys.argv[6] # This is the base network
ha_topnet_fname = sys.argv[7]

# Read microarray data
# Column 0 -> gene labels. Make this the index
# Columns 1 through m -> various samples
SI = pd.read_csv(data_fname, sep = '\t', na_values=' ')
SI = SI.set_index(SI.columns[0])
SI = SI[[sample_of_interest]]
#SI = mic_fun.restructure_SI(SI, sample_of_interest, control_sample) # This gives SI restructured such that
								     # column 0 -> perturbation to study, 
								     # column 1 -> control, 
								     # columns 2 to m -> other perturbations
print("Got microarray data with ", SI.shape[0], " rows and ", SI.shape[1], " columns")

# Read unweighted network
G_unweighted = nx.read_edgelist(unweighted_nw_fname, delimiter = "\t", nodetype = str, create_using = nx.DiGraph())
print("Read network with ", len(G_unweighted.nodes()), " nodes and ", len(G_unweighted.edges()), " edges")

# Compute Pij
Pij, G_ha = combine_data_get_sp_paths_costs_ha(G_unweighted, SI, ha_nw_fname) # Give an empty list as the paths argument
Pij = perc_fun.get_Pij_percentile_norm_cost(Pij, percentile, path_length_thresh)
print("After taking percentile cutoff, Pij has ", Pij.shape[0], " rows and ", Pij.shape[1], " columns")
print(Pij.head())

# Get edges in the top paths
topnet_edges = set()
for path_str in Pij.index:
	path = path_str.split('#')
	topnet_edges = topnet_edges.union(get_edges_in_path(path, G_ha))
print("Got ", len(topnet_edges), " edges in TopNet")

# Write edges in significant paths
with open(ha_topnet_fname, 'w') as outfile:
	for edge in topnet_edges:
		outfile.write( edge[0] + "\t" )
		outfile.write( edge[1] + "\t" )
		outfile.write( str(edge[2]) + "\n" )
print("Done writing TopNet")
