import pandas as pd
import math
import numpy as np
import random


############## Methods for the case with a single perturbed and a single control sample,
############## but with multiple perturbations

# Create a new dataframe with the same index as the input, and with the specified column name
# For each row, pick a random value from that row
def get_randomized_colnum(SI, colnum):
	randomized_vals = pd.DataFrame(index = tuple(SI.index), columns = [colnum])
	for row in SI.itertuples():
		list_row = list(row)[1:]
		list_row = [x for x in list_row if not math.isnan(x)] # remove NaNs
		randomized_vals.loc[row[0]] = random.choice(list_row)
	return randomized_vals

def get_randomized_single_sample_mult_pert(SI):
	randomized_perturbed = get_randomized_colnum(SI, 0)
	randomized_control = get_randomized_colnum(SI, 1)
	return pd.concat([randomized_perturbed, randomized_control], axis=1)

# We are only interested in column 0 (perturbed) and column 1 (control)
# So we return a pandas dataframe with only the relevant columns
def get_relevant_SI(SI):
	perturbed_data = SI.iloc[:,0]
	control_data = SI.iloc[:,1]
	SI_relevant = pd.concat([perturbed_data, control_data], axis=1)
	return SI_relevant

def restructure_SI(SI, perturbation_sample, control_sample):
	perturbed_data = SI.loc[:, perturbation_sample]
	control_data = SI.loc[:, control_sample]
	SI = SI.drop([perturbation_sample, control_sample], axis = 1)
	SI_restructured = pd.concat([perturbed_data, control_data, SI], axis = 1)
	return SI_restructured




############## Methods for the case with multiple perturbed and control samples

# Return a dataframe with each row shuffled around independently
def shuffle_disease_healthy(SI):
        shuffled_SI = pd.DataFrame(index = tuple(SI.index), columns = tuple(SI.columns))
        for row in SI.itertuples():
                list_row = list(row)[1:]
                shuffled_SI.loc[row[0]] = np.random.permutation(list_row)
        return shuffled_SI

# Here we have k1 perturbed samples and k2 control samples
# k1 + k2 = m (total number of samples)
# We return the median of the first k1 columns as the representative perturbed sample
# and the median of the next k2 columns as the representative control sample
def get_median_SI(SI, num_perturbed_samples, num_control_samples):
        perturbed_median = SI.iloc[:,0:num_perturbed_samples].median(1)
        control_median = SI.iloc[:,num_perturbed_samples:(num_control_samples+num_perturbed_samples)].median(1)
        SI_median = pd.concat([perturbed_median, control_median], axis=1)
        return SI_median
