# Code to create an activated response network given inputs PPI and gene expression data,
# and find top x percentile shortest paths in the actual data.
# We also create a background of randomized data by permuting the input gene expression data.
# For each resulting randomized network, we compute the cost of the paths that are in the
# top x percentile in the actual data.

# This code assumes the omics data is already normalized, and has the following form
# Column 0 -> gene symbols
# Columns 1 through m -> omics values
# Each column can correspond to a different condition (eg: individual patient samples),
# or a summarised value (eg: median expression across a cohort)
# Column names for the control and the perturbation of interest are given as inputs

import pandas as pd
import random
random.seed(1)
import networkx as nx
import numpy as np
import sys
import math

import microarray_functions as mic_fun
import network_functions as net_fun
import percentile_functions as perc_fun

# SI is given as a pandas dataframe, indexed by gene
# column 0 -> disease gene expression values
# column 1 -> healthy gene expression values
def combine_data_get_sp_paths_costs_activated(G_unweighted, SI, response_nw_fname, paths):
	SI_relevant = mic_fun.get_relevant_SI(SI)

	# Map gene expression values onto unweighted network
	G_response = net_fun.get_activated_response_network(SI_relevant, G_unweighted)
	print("Got response network with ", len(G_response.nodes()), " nodes and ", len(G_response.edges()), " edges")

	# Drop SI values for genes which don't map to response network
	genes_to_drop = set(SI.index) - set(G_response.nodes())
	SI = SI.drop(genes_to_drop)

	if len(paths) == 0: # Actual data
		# Write the response network to file
		nx.write_weighted_edgelist(G_response, response_nw_fname, delimiter = '\t')
		# Get all-pairs-shortest-path costs
		# Return value is a pandas dataframe, indexed by the string src#dest
		Pij = net_fun.get_all_sp_paths_costs(G_response)
		print("Got shortest path costs for ", Pij.shape[0], " node-pairs")
	elif len(paths) > 0: # Randomized data
		# Get cost of the same paths that are the shortest in the actual dataset
		Pij = net_fun.get_costs_of_given_paths(G_response, paths)
		print("Got cost of paths which are shortest in the actual data")

	return Pij, SI


if len(sys.argv) != 10:
	print("argv[1] = microarray data file (tab-delimited, with header)")
	print("argv[2] = name of perturbation sample to study")
	print("argv[3] = name of control sample")
	print("argv[4] = unweighted (directed) network file")
	print("argv[5] = percentile threshold")
	print("argv[6] = path length threshold")
	print("argv[7] = number of randomizations")
	print("argv[8] = output file for activated response base network")
	print("argv[9] = output file prefix for Pij")
	sys.exit(1)

# Set inputs
data_fname = sys.argv[1]
perturbation_sample = sys.argv[2] # This is the perturbation we want to study
control_sample = sys.argv[3]
unweighted_nw_fname = sys.argv[4]
percentile = float(sys.argv[5]) # We'll only keep paths whose cost < this threshold
path_length_thresh = int(sys.argv[6]) # We'll only keep paths with length >= this threshold
num_trials = int(sys.argv[7])
response_nw_fname = sys.argv[8]
output_fname_prefix = sys.argv[9]

# Read microarray data
# Column 0 -> gene labels. Make this the index
# Columns 1 through m -> a control, and various perturbed conditions
SI = pd.read_csv(data_fname, sep = '\t', na_values=' ')
SI = SI.set_index(SI.columns[0])
SI = mic_fun.restructure_SI(SI, perturbation_sample, control_sample) # This gives SI restructured such that
								     # column 0 -> perturbation to study, 
								     # column 1 -> control, 
								     # columns 2 to m -> other perturbations
print("Read microarray data with ", SI.shape[0], " rows and ", SI.shape[1], " columns")

# Read unweighted network
G_unweighted = nx.read_edgelist(unweighted_nw_fname, delimiter = "\t", nodetype = str, create_using = nx.DiGraph())
print("Read network with ", len(G_unweighted.nodes()), " nodes and ", len(G_unweighted.edges()), " edges")

# Compute Pij in actual data
Pij, SI = combine_data_get_sp_paths_costs_activated(G_unweighted, SI, response_nw_fname, []) # Give an empty list as the paths argument
print("After dropping genes which don't map to our network, got SI for ", SI.shape[0], " genes and ", SI.shape[1], " samples")
Pij = perc_fun.get_Pij_percentile_norm_cost(Pij, percentile, path_length_thresh)
print("After taking percentile cutoff, Pij has ", Pij.shape[0], " rows and ", Pij.shape[1], " columns")
print(Pij.head())
Pij.to_csv(output_fname_prefix+"_actual.txt", sep = "\t")

# Randomized data
for i in range(num_trials):
	print("######################## Trial ", i, " ########################")

	randomized_SI = mic_fun.get_randomized_single_sample_mult_pert(SI)
	print("After shuffling, got SI values for ", randomized_SI.shape[0], " genes and ", randomized_SI.shape[1], " samples")

	# Give shortest paths in actual dataset as the paths argument
	Pij, randomized_SI = combine_data_get_sp_paths_costs_activated(G_unweighted, randomized_SI, response_nw_fname, list(Pij.index))

	Pij.to_csv(output_fname_prefix+"_"+str(i)+".txt", sep = "\t")
