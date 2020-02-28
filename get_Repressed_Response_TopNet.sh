#!/bin/sh

trim_last_slash () {
        wdir=$1
        if [ "${wdir: -1}" == "/" ]
        then
                wdir="${wdir%?}"
        fi
        echo "${wdir}"
}

if [ "$#" -ne 11 ]
then
	echo "argv[1] = microarray data file"
	echo "argv[2] = name of perturbation sample to study"
	echo "argv[3] = name of control sample"
	echo "argv[4] = unweighted network file"
	echo "argv[5] = percentile threshold"
	echo "argv[6] = path length threshold"
	echo "argv[7] = q-score cutoff"
	echo "argv[8] = number of randomizations"
	echo "argv[9] = output directory"
	echo "argv[10] = file name for base response network (we'll put it in the output directory)"
	echo "argv[11] = file name for TopNet (we'll put it in the output directory)"
        exit
fi

# Set inputs
data_fname=${1}
perturbation_sample=${2}
control_sample=${3}
unweighted_nw_fname=${4}
percentile=${5}
path_length_thresh=${6}
qscore_thresh=${7}
num_trials=${8}
out_dir=$(trim_last_slash ${9})
response_nw_fname=${10}
topnet_fname=${11}

# Create a temporary directory in the output directory
# We will put all the Pij files here
# Once the z-score has been calculated we will remove the Pij files as well as the temp directory
mkdir ${out_dir}/temp

# Map the microarray data onto the unweighted network to get a repressed response network.
# Compute top 'percentile' shortest paths in this (actual) network. Write these paths and costs.
# Randomize the microarray data 'num_trials' times, resulting in 'num_trials' randomized response networks.
# Compute the cost of the same paths as in the top 'percentile' shortest paths in the actual response network.
# Write these paths and costs into output files.
python repressed_response_Pijs.py ${data_fname} ${perturbation_sample} ${control_sample} ${unweighted_nw_fname} ${percentile} ${path_length_thresh} ${num_trials} ${out_dir}/${response_nw_fname} ${out_dir}/temp/Pij

# Calculate z-score and corresponding p-value for each path
python fdr_rand_pijs_boxcox.py ${out_dir}/temp ${out_dir}/Pij_zscores.txt

# Delete all temporary files
rm -rf ${out_dir}/temp

# Carry out FDR (benjamini-hochberg)
python benjamini_hochberg_boxcox.py ${out_dir}/Pij_zscores.txt ${qscore_thresh} ${out_dir}/Pij_zscores_fdr.txt

# Extract top-net based on paths whose q-score <= qscore_thresh
python extract_fdr_network.py ${out_dir}/${response_nw_fname} ${out_dir}/Pij_zscores_fdr.txt ${qscore_thresh} ${out_dir}/${topnet_fname}
