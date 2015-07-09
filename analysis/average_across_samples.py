#!/usr/bin/python

# ---------------------------------------------------------------------------------
# Handling-missingness

# This script will compile results across samples, and take an average
import os
import pandas
import sys
import pickle
import numpy as np
import nibabel as nib
from glob import glob

top_directory = os.path.abspath("../")
inputs_directory = "%s/permutation_example" %(top_directory)
output_directory = "%s/results" %(top_directory)
contrast_file = "%s/doc/hcp_contrasts_id_filter.csv" %(top_directory)
contrasts = pandas.read_csv(contrast_file,sep="\t").id.tolist()
os.chdir(output_directory)

# Parse data files

pearsons_cca = glob("*pearson_cca.tsv")
pearsons_svi = glob("*pearson_svi.tsv")
spearmans_cca = glob("*spearman_cca.tsv")
spearmans_svi = glob("*spearman_svi.tsv")

inputs = {"pearsons_cca":pearsons_cca,
          "pearsons_svi":pearsons_svi,
          "spearmans_svi":spearmans_svi,
          "spearmans_cca":spearmans_cca}

for input_key,files in inputs.iteritems():
    print "Processing %s" %input_key
    # Read into list of data frames
    dfs = numpy.zeros((1316,47,len(files)))
    for ff in range(0,len(files)):
        f = files[ff]
        number =  f.split("_")[0]
        print number
        tmp = pandas.read_csv(f,sep="\t",index_col=0)
        # To merge data frames we need same rows and columns
        tmp.index = [x.replace("%s_TASK" %number,"TASK") for x in tmp.index]
        tmp.columns = [x.replace("%s_TASK" %number,"TASK") for x in tmp.columns]
        # This info is represented in index names
        tmp = tmp.drop("thresh",axis=1)
        tmp = tmp.drop("direction",axis=1)
        dfs[:,:,ff] = tmp
    # Now take mean across the 500 dataframes
    df = numpy.zeros((1316,47))
    for i in range(0,1316):
       for j in range(0,47):
           df[i,j] = numpy.mean(dfs[i,j,:])
    df = pandas.DataFrame(df)
    df.index = tmp.index
    df.columns = tmp.columns
    df.to_csv("%s.tsv" %input_key,sep="\t")


# Parse size files
sizes_cca = glob("*sizes_cca.tsv")
sizes_svi = glob("*sizes_svi.tsv")
inputs = {"sizes_cca":sizes_cca,"sizes_svi":sizes_svi}

for input_key,files in inputs.iteritems():
    print "Processing %s" %input_key
    # Read into list of data frames
    dfs = numpy.zeros((1316,47,len(files)))
    for ff in range(0,len(files)):
        f = files[ff]
        number =  f.split("_")[0]
        print number
        tmp = pandas.read_csv(f,sep="\t",index_col=0)
        # To merge data frames we need same rows and columns
        tmp.index = [x.replace("%s_TASK" %number,"TASK") for x in tmp.index]
        tmp.columns = [x.replace("%s_TASK" %number,"TASK") for x in tmp.columns]
        # This info is represented in index names
        tmp = tmp.drop("thresh",axis=1)
        tmp = tmp.drop("direction",axis=1)
        dfs[:,:,ff] = tmp
    # Now take mean across the 500 dataframes
    df = numpy.zeros((1316,47))
    for i in range(0,1316):
       for j in range(0,47):
           sizes = numpy.mean(dfs[i,j,:])
           if (sizes==0).any():
               df[i,j] = 0.0
           else:
               df[i,j] = numpy.mean(dfs[i,j,:])
    df = pandas.DataFrame(df)
    df.index = tmp.index
    df.columns = tmp.columns
    df.to_csv("%s.tsv" %input_key,sep="\t")
