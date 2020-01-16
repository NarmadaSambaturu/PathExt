import pandas as pd
import statsmodels.stats.multitest as bh
import sys

if len(sys.argv) != 4:
	print("argv[1] = file with i#j and z-score p-values")
	print("argv[2] = alpha value (family-wise error rate)")
	print("argv[3] = output file")
	sys.exit(1)

# Set inputs
zscore_pvals_fname = sys.argv[1]
alpha_val = float(sys.argv[2])
out_fname = sys.argv[3]

# Read z-scores and corresponding p-values
zscore_pvals = pd.read_table(zscore_pvals_fname)
zscore_pvals.columns = ('ij', 'normaltest_before_boxcox', 'normaltest_after_boxcox', 'zscore_after_boxcox', 'two_sided_pval_for_zscore')
zscore_pvals = zscore_pvals.drop(['normaltest_before_boxcox', 'normaltest_after_boxcox'], axis = 1)
zscore_pvals = zscore_pvals.set_index(zscore_pvals.columns.values[0])
zscore_pvals = zscore_pvals.dropna()
print(zscore_pvals.head())
print("Done reading input file")
print("After dropping nans, shape of dataframe is ", zscore_pvals.shape)

# Carry out Benjamini-Hochberg p-value correction
bh_output = bh.multipletests(zscore_pvals['two_sided_pval_for_zscore'], alpha = alpha_val, method = 'fdr_bh')
reject = bh_output[0]
bh_pval = bh_output[1]
print("Done with BH test")

# Print output
with open(out_fname, 'w') as f:
	f.write("i#j" + "\t")
	f.write("zscore_after_boxcox" + "\t")
	f.write("zscore_pval" + "\t")
	f.write("bh_corrected_pval" + "\t")
	f.write("reject" + "\n")
	index = 0
	for row in zscore_pvals.itertuples():
		f.write(row[0] + "\t")
		f.write(str(row[1]) + "\t")
		f.write(str(row[2]) + "\t")
		f.write(str(bh_pval[index]) + "\t")
		f.write(str(reject[index]) + "\n")
		index += 1
