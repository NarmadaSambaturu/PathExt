import numpy as np

def get_Pij_percentile_norm_cost(Pij, percentile, path_length_thresh):
	# Add new columns 1 -> path length (number of hops)
	# 2 -> normalised cost (cost/path length)
	# Only retain paths whose length >= path_length_thresh
	Pij[1] = [len(x.split('#')) - 1 for x in list(Pij.index)]
	Pij[2] = Pij[0]/Pij[1]
	Pij = Pij.loc[Pij[1] >= path_length_thresh]

	# Sort by normalized cost. Normalization is by dividing path cost by number of hops
	#Pij[2] = Pij[0]/Pij[1]
	Pij = Pij.sort_values(by = Pij.columns.values[-1])
	#print("Done sorting by normalized cost")

	# Get cost corresponding to percentile
	costs = list(Pij[2])
	path_cost_thresh = np.percentile(costs, percentile)
	#print("Path cost threshold = ", path_cost_thresh)

	# Keep only rows with length >= path_length_thresh AND cost < path_cost_thresh
	Pij = Pij.loc[Pij[2] < path_cost_thresh]
	#print("Done retaining only ijs with path length >= ", path_length_thresh, " AND path cost < ",  path_cost_thresh)

	# Drop the column with normalized costs
	Pij.drop([1, 2], axis = 1, inplace = True)
	#print("Done dropping the column with normalized costs")

	return Pij
