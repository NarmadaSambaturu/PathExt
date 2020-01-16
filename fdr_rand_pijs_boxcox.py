import sys
import scipy.stats
import scipy.stats.mstats
import numpy as np
import pandas as pd
import os

def concat_randomised_pij_values(pij, rand_pij_files):
	i = 0
	for f in rand_pij_files:
		print(f)
		temp = pd.read_csv(f, sep = '\t')
		temp = temp.set_index(temp.columns.values[0])
		temp.columns = [str(i)]
		pij = pd.concat([pij, temp], axis=1)
		i += 1
	return pij

def read_pij_files(dir_name):
	# Get all pij fnames
	actual_data_fname = ""
	rand_pij_files = []
	for f in os.listdir(dir_name):
		if os.path.isfile(os.path.join(dir_name, f)) and f.endswith(".txt"):
			if "actual" in f:
				actual_data_fname = os.path.join(dir_name, f)
			else:
				rand_pij_files.append(os.path.join(dir_name, f))
	print("Found ", len(rand_pij_files), " files")
	print("Actual data fname = ", actual_data_fname)

	# Read actual data
	pij = pd.read_table(actual_data_fname)
	pij = pij.set_index(pij.columns.values[0]) # This allows us to access values by i#j
	pij.columns = ['actual']
	print("Actual pij has ", pij.shape[0], " rows and ", pij.shape[1], " columns")
	print(pij.head())

	# Read randomisations and concat to pij dataframe
	pij = concat_randomised_pij_values(pij, rand_pij_files)
	print("After reading all randomized pijs, we have ", pij.shape[0], " rows and ", pij.shape[1], " columns")
	print(pij.head())
	return pij

# We get a value, and a list as inputs
# Apply Box-Cox transformation to the list to make sure it is normally distributed
# Carry out normaltest to verify that it is normally distributed
# Get z-score of that value in that list
# Get p-value for that z-score, with the assumption that the distribution is normally distributed
# Return tuple (normaltest before Box-Cox, normaltest after Box-Cox, z-score, two sided p-value of z-score)
def get_zscore_pval_boxcox(val, list_vals):
	list_vals.append(val)
	normaltest_before_boxcox = scipy.stats.mstats.normaltest(list_vals)
	list_vals_boxcox = scipy.stats.boxcox(list_vals)[0]
	normaltest_after_boxcox = scipy.stats.mstats.normaltest(list_vals_boxcox)
	new_val = list_vals_boxcox[-1]
	list_vals = list_vals[:-1]
	zscore = (new_val - np.mean(list_vals_boxcox))/np.std(list_vals_boxcox)
	two_sided_zscore_pval = scipy.stats.norm.sf(abs(zscore))*2
	return (normaltest_before_boxcox[1], normaltest_after_boxcox[1], zscore, two_sided_zscore_pval)


if len(sys.argv) != 3:
	print("argv[1] = folder with all pij files")
	print("argv[2] = output file")
	sys.exit(1)

# Set inputs
dir_name = sys.argv[1]
out_fname = sys.argv[2]

################# Read pij files #################
pij = read_pij_files(dir_name)

################# Compute z-scores #################
zscores = {} # key -> node1#node2, value -> tuple (normaltest before Box-Cox, normaltest after Box-Cox, z-score, two sided p-value of z-score)
rownum = 0
for row in pij.itertuples():
        l = list(row[1:])
        nodes = row[0]. split('#')
        if nodes[0] == nodes[1]:
                continue
        zscores[row[0]] = get_zscore_pval_boxcox(l[0], l[1:])
        if rownum % 10000 == 0:
                print("Done working with ", rownum, "rows")
        rownum += 1
print("Done computing zscores")

# Print output
with open(out_fname, 'w') as outfile:
	outfile.write( "ij" + "\t" )
	outfile.write( "normaltest_before_boxcox" + "\t" )
	outfile.write( "normaltest_after_boxcox" + "\t" )
	outfile.write( "zscore_after_boxcox" + "\t" )
	outfile.write( "two_sided_pval_for_zscore" + "\n" )
	for ij, tup in zscores.items():
		outfile.write( ij + "\t" )
		outfile.write( "\t".join([str(x) for x in tup]) + "\n" )
