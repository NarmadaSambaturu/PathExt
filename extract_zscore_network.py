import networkx as nx
import csv
import sys
import pandas as pd

def get_edges_in_path(path, G):
	edges_in_path = set()
	for i in range(len(path) - 1):
		edges_in_path.add((path[i], path[i+1], G[path[i]][path[i+1]]['weight']))
	return edges_in_path

if len(sys.argv) != 5:
	print("argv[1] = response network (weighted)")
	print("argv[2] = pij file with zscore p-values and BH correction")
	print("argv[3] = q-score cutoff")
	print("argv[4] = output file")
	sys.exit(1)

nw_fname = sys.argv[1]
zscores_fname = sys.argv[2]
pval_cutoff = float(sys.argv[3])
out_fname = sys.argv[4]

# Read response network
G_response = nx.read_weighted_edgelist(nw_fname, delimiter = "\t", create_using = nx.DiGraph())
print("Got response network with ", len(G_response.nodes()), " nodes and ", len(G_response.edges()), " edges")

# Read zscores file
# Only retain information where BH-adjusted p-value < threshold
significant_paths = set() # node1#node2#...#nodek
with open(zscores_fname) as f:
	reader = csv.reader(f, delimiter = "\t")
	next(reader) # skip header
	for row in reader:
		if row[3] == 'NA':
			continue
		path_str = row[0]
		BH_pval = float(row[3])
		zscore = float(row[1])
		if BH_pval <= pval_cutoff: #and zscore < 0: # actual Pij < randomised Pij
			significant_paths.add(path_str)
print(len(significant_paths), " node pairs have significant zscores")

################ Make a network with the edges involved in the paths Pij

# Get edges in the significant paths
topnet_edges = set()
for path_str in significant_paths:
	path = path_str.split('#')
	topnet_edges = topnet_edges.union(get_edges_in_path(path, G_response))
print("Got ", len(topnet_edges), " edges in TopNet")

# Write edges in significant paths
with open(out_fname, 'w') as outfile:
	for edge in topnet_edges:
		outfile.write( edge[0] + "\t" )
		outfile.write( edge[1] + "\t" )
		outfile.write( str(edge[2]) + "\n" )
print("Done writing TopNet")
