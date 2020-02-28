# PathExt: a general framework for path-based mining of omics-integrated biological networks

We provide a computational tool - PathExt, which identifies differentially active paths when a control is available, and most active paths otherwise, in an omics-integrated biological network. The sub-network comprising such paths, referred to as the TopNet, captures the most relevant genes and processes underlying the specific biological context. The TopNet forms a well-connected graph, reflecting the tight orchestration in biological systems. Two key advantages of PathExt are (i) it can extract characteristic genes and pathways even when only a single sample is available, and (ii) it can be used to study a system even in the absence of an appropriate control.

The inputs to PathExt are (a) a directed gene network and (b) gene-centric omics data for the conditions of interest. The omics data can represent a variety of quantities pertaining to the node, such as gene expression level, differential expression, protein, metabolite level, etc., in one or more conditions. The output of PathExt is a sub-network, that we refer to as the TopNet, consisting of the most significant differential or active paths, and is interpreted based on the application context.

***************************************************************************************
Executing a command without any arguments will generate a list of required inputs.

# Comparing two conditions

#### Activated Response TopNet - captures up-regulated genes and processes <br />

$ bash get_Activated_Response_TopNet.sh <br />
argv[1] = microarray data file <br />
argv[2] = name of perturbation sample to study <br />
argv[3] = name of control sample <br />
argv[4] = unweighted network file <br />
argv[5] = percentile threshold <br />
argv[6] = path length threshold <br />
argv[7] = q-score cutoff <br />
argv[8] = number of randomizations <br />
argv[9] = output directory <br />
argv[10] = file name for base response network (we'll put it in the output directory) <br />
argv[11] = file name for TopNet (we'll put it in the output directory) <br />

Inputs:
1. microarray data file <br />
Tab-delimited file of normalized signal intensities. Each column can correspond to a different condition (eg: individual patient samples), or a summarised value (eg: median expression across a cohort). 

2. name of perturbation sample to study <br />
3. name of control sample <br />
Column names for the perturbation of interest (eg: a particular patient ID) and control are given as inputs.

4. unweighted network file <br />
Knowledge-based network

5. percentile threshold <br />
Paths whose cost is below this percentile threshold will be used to compute the TopNet.

6. path length threshold <br />
Minimum number of edges a path should have for it to be considered for TopNet creation. A good default value is 2.

7. q-score cutoff <br />
FDR-corrected cutoff to be used to consider a path significantly different from a randomized background.

Output files generated: <br />
1. Activated Response base network <br />
Base network after integrating omics data with knowledge based network. Node weight used is N_i = SI x FC where N_i is the weight of node i, and SI is the normalized signal intensity, or expression level, of a particular gene. FC = SI_perturbed/SI_control is the fold change in expression values. Edge cost = 1/sqrt(N_i x N_j). Tab-delimited file.

2. Activated Response TopNet <br />
TopNet comprising of edges from highly up-regulated and statistically significant paths. Tab-delimited file.

3. Pij_zscores_fdr.txt <br />
Tab-delimited file with details of paths considered in TopNet creation (after applying path length and percentile thresholds). This file shows the z-score of the actual path cost vs randomized path costs, as well as the FDR-corrected q-score for each such path.

Example: <br />
$ bash get_Activated_Response_TopNet.sh test_data/GSE71200_SI.txt GSM1829740 GSM1829696 test_data/small_Mtb_network.txt 0.5 2 0.05 100 test_data/results/ Activated_Response_base_network.txt Activated_Response_TopNet.txt

The example input omics data, test_data/GSE71200_SI.txt, is taken from GEO accession number GSE71200 [1]. It contains transcriptomic data for mycobacterium tuberculosis (M.tb) exposed for 16 hours to various concentrations of antibiotics.
GSM1829740 corresponds to M.tb exposed to twice the minimum inhibitory concentration of isoniazid. <br />
The example input knowledge-based network, test_data/small_Mtb_network.txt, is a sub-network of the network published in [2].



#### Repressed Response TopNet - captures down-regulated genes and processes <br />

$ bash get_Repressed_Response_TopNet.sh <br />
argv[1] = microarray data file <br />
argv[2] = name of perturbation sample to study <br />
argv[3] = name of control sample <br />
argv[4] = unweighted network file <br />
argv[5] = percentile threshold <br />
argv[6] = path length threshold <br />
argv[7] = q-score cutoff <br />
argv[8] = number of randomizations <br />
argv[9] = output directory <br />
argv[10] = file name for base response network (we'll put it in the output directory) <br />
argv[11] = file name for TopNet (we'll put it in the output directory) <br />

Example: <br />
$ bash get_Repressed_Response_TopNet.sh test_data/GSE71200_SI.txt GSM1829740 GSM1829696 test_data/small_Mtb_network.txt 0.5 2 0.05 100 test_data/results/ Repressed_Response_base_network.txt Repressed_Response_TopNet.txt



#### Response TopNet - union of Activated and Repressed Response TopNets, provides a holistic view of the active, altered genes and processes. <br />

$ python get_union_response_TopNet.py <br />
argv[1] = activated response TopNet (weighted) <br />
argv[2] = repressed respone TopNet (weighted) <br />
argv[3] = output file for union response TopNet (unweighted) <br />

Example: <br />
$ python get_union_response_TopNet.py test_data/results/Activated_Response_TopNet.txt test_data/results/Repressed_Response_TopNet.txt test_data/results/Response_TopNet.txt




***************************************************************************************

# Studying a single condition, or working in the absence of an appropriate control

#### Highest Activity TopNet (HA TopNet) - captures highly active processes in the condition of interest

$ python get_highest_activity_TopNet.py <br />
argv[1] = microarray data file (tab-delimited, with header) <br />
argv[2] = name of sample to study <br />
argv[3] = unweighted (directed) network file <br />
argv[4] = percentile threshold <br />
argv[5] = path length threshold <br />
argv[6] = output file for highest activity base network <br />
argv[7] = output file for HA TopNet <br />

Example: <br />
$ python get_highest_activity_TopNet.py test_data/GSE71200_SI.txt GSM1829740 test_data/small_Mtb_network.txt 0.5 2 test_data/results/HA_base_network.txt test_data/results/HA_TopNet.txt



***************************************************************************************

If you find this code useful, please cite <br />
Sambaturu, Narmada, et al. "PathExt: a general framework for path-based mining of omics-integrated biological networks." bioRxiv (2020).


### References

[1] Ma S, Minch KJ, Rustad TR, Hobbs S et al. Integrated Modeling of Gene Regulatory and Metabolic Networks in Mycobacterium tuberculosis. PLoS Comput Biol 2015 Nov;11(11):e1004543. PMID: 26618656 <br />
[2] Mishra, S., Shukla, P., Bhaskar, A., Anand, K., Baloni, P., Jha, R. K., Mohan, A., Rajmani, R. S., Nagaraja, V., Chandra, N., et al. (2017). Efficacy of β-lactam/β-lactamase inhibitor combination is linked to whib4-mediated changes in redox physiology of mycobacterium tuberculosis. Elife, 6, e25624.
