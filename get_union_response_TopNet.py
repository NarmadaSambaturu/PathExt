import networkx as nx
import sys

if len(sys.argv) != 4:
	print("argv[1] = activated response TopNet (weighted)")
	print("argv[2] = repressed respone TopNet (weighted)")
	print("argv[3] = output file for union response TopNet (unweighted)")
	sys.exit(1)

# Set inputs
act_topnet_fname = sys.argv[1]
rep_topnet_fname = sys.argv[2]
response_topnet_fname = sys.argv[3]

# Read activated and repressed response TopNets
act_topnet = nx.read_weighted_edgelist(act_topnet_fname, delimiter = '\t', create_using = nx.DiGraph())
print("Got Activated Response TopNet with ", len(act_topnet.nodes()), " nodes and ", len(act_topnet.edges()), " edges")

rep_topnet = nx.read_weighted_edgelist(rep_topnet_fname, delimiter = '\t', create_using = nx.DiGraph())
print("Got Repressed Response TopNet with ", len(rep_topnet.nodes()), " nodes and ", len(rep_topnet.edges()), " edges")

# Get union network
response_topnet = nx.DiGraph()
response_topnet.add_edges_from(set(act_topnet.edges()).union(set(rep_topnet.edges())))
print("Got union Response TopNet with ", len(response_topnet.nodes()), " nodes and ", len(response_topnet.edges()), " edges")

# Write union network to file
nx.write_edgelist(response_topnet, response_topnet_fname, delimiter = '\t', data = False)
print("Done writing Response TopNet to file")
