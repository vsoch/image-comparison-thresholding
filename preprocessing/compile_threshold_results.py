#!/usr/bin/python

# -------------------------------------------------------------
# Handling-missingness:  compile threshold permutation results

# This script will compile results from piecewise analyses, likely needs big memory node
import os
import pandas
import sys
import pickle
import numpy as np
import nibabel as nib
from glob import glob

top_directory = os.path.abspath("../")
inputs_directory = "%s/data/permutation_example" %(top_directory)
output_directory = "%s/data/results" %(top_directory)
contrast_file = "%s/doc/hcp_contrasts_id_filter.csv" %(top_directory)
contrasts = pandas.read_csv(contrast_file,sep="\t").id.tolist()

# Lets work from the output directory
os.chdir(output_directory)
notdone = []

# 500 was done for experiment
nruns = len("%s/*" %(glob(inputs_directory)))
for r in range(0,nruns):
    print "Processing %s" %(r)
    inputs = np.sort(glob("%s/%s/comparisons/*.pkl" %(inputs_directory,r))).tolist()
    input_ids = ["%s_%s" %(r,x.replace(".pkl","").replace("%s/%s/comparisons/" %(inputs_directory,r),"")) for x in inputs]
    # Find contrast ids that are missing
    contrast_ids = ["%s" %(x.replace(".pkl","").replace("%s/%s/comparisons/" %(inputs_directory,r),"")) for x in inputs]
    missing = ["%s_%s" %(r,c) for c in contrasts if c not in contrast_ids]
    if len(missing) > 0:
        print "ERROR: Missing output for run %s, re-run when runs available." %(r)
        notdone.append(r)
    else:
        # Read in the first input to get threshold ids, etc.
        tmp = pickle.load(open(inputs[0],"rb"))
        # Our dataframe index must have image id, thresh, and direction (to be unique)
        row_index = []
        thresholds = np.unique(tmp["thresh"])
        directions = ["posneg","pos"]
        for i in input_ids:
            for t in thresholds:
                for d in directions:
                    row_index.append("%s_%s_%s" %(i,t,d))
        # We will save all results to these lists
        pearsons_pi = pandas.DataFrame(index=row_index) 
        pearsons_pd = pandas.DataFrame(index=row_index)     # pearsonr scores
        spearmans_pi = pandas.DataFrame(index=row_index) 
        spearmans_pd = pandas.DataFrame(index=row_index)    # spearmanr scores
        pd_sizes = pandas.DataFrame(index=row_index) 
        pi_sizes = pandas.DataFrame(index=row_index)        # size of final masks
        nanlog_pd = pandas.DataFrame(index=row_index) 
        nanlog_pi = pandas.DataFrame(index=row_index)       # "success","nan_fewer_3_values","nan_no_overlap"
        # And we will keep a column of thresholds, and one of directions
        pearsons_pi["direction"] = [x.split("_")[4] for x in row_index]
        pearsons_pd["direction"] = [x.split("_")[4] for x in row_index]
        spearmans_pi["direction"] = [x.split("_")[4] for x in row_index]
        spearmans_pd["direction"] = [x.split("_")[4] for x in row_index]
        pi_sizes["direction"] = [x.split("_")[4] for x in row_index]
        pd_sizes["direction"] = [x.split("_")[4] for x in row_index]
        nanlog_pd["direction"] = [x.split("_")[4] for x in row_index]
        nanlog_pi["direction"] = [x.split("_")[4] for x in row_index]
        pearsons_pi["thresh"] = [x.split("_")[3] for x in row_index]
        pearsons_pd["thresh"] = [x.split("_")[3] for x in row_index]
        spearmans_pi["thresh"] = [x.split("_")[3] for x in row_index]
        spearmans_pd["thresh"] = [x.split("_")[3] for x in row_index]
        nanlog_pd["thresh"] = [x.split("_")[3] for x in row_index]
        nanlog_pi["thresh"] = [x.split("_")[3] for x in row_index]
        pi_sizes["thresh"] = [x.split("_")[3] for x in row_index]
        pd_sizes["thresh"] = [x.split("_")[3] for x in row_index]
        for i in inputs:
            tmp = pickle.load(open(i,"rb"))
            # Each image is an entire column, with threshold and image_ids in rows
            input_id = "%s_%s" %(r,tmp["id"]) # This will be the column name
            uids = ["%s_%s_%s" %(r,tmp["idB"][x],tmp["size_ids"][x]) for x in range(0,len(tmp["idB"]))]
            pearsons_pi.loc[uids,input_id] = tmp["svi_pearson"]
            pearsons_pd.loc[uids,input_id] = tmp["cca_pearson"]
            spearmans_pi.loc[uids,input_id] = tmp["svi_spearman"]
            spearmans_pd.loc[uids,input_id] = tmp["cca_spearman"]
            # Save all mask sizes differences, again will be nan for image vs itself.
            pd_sizes.loc[uids,input_id] = tmp["sizes"]["cca"].tolist()
            pi_sizes.loc[uids,input_id] = tmp["sizes"]["svi"].tolist()
            # Nanlog - did we append a nan and why?
            nanlog_pd.loc[uids,input_id] = tmp["nanlog_cca"]
            nanlog_pi.loc[uids,input_id] = tmp["nanlog_svi"]
        pearsons_pd.to_csv("%s_pearson_cca.tsv" %(r),sep="\t")
        pearsons_pi.to_csv("%s_pearson_svi.tsv" %(r),sep="\t")
        spearmans_pd.to_csv("%s_spearman_cca.tsv" %(r),sep="\t")
        spearmans_pi.to_csv("%s_spearman_svi.tsv" %(r),sep="\t")
        pd_sizes.to_csv("%s_sizes_cca.tsv" %(r),sep="\t")
        pi_sizes.to_csv("%s_sizes_svi.tsv" %(r),sep="\t")
        nanlog_pd.to_csv("%s_nanlog_cca.tsv" %(r),sep="\t")
        nanlog_pi.to_csv("%s_nanlog_svi.tsv" %(r),sep="\t")
