#!/usr/bin/python

# ---------------------------------------------------------------------------------
# handling-missingness: test threshold instance

# Here I want to do a simple experiment to show how pearson and spearman correlations change with 
# two strategies for handling missing data.

# complete case analysis [intersection, voxels in both mask]
# single-value imputation [union, voxels in either mask]

import os
import time
import pandas
import sys
import pickle
import numpy as np
import nibabel as nib
from glob import glob
import image_transformations as IT
from scipy.stats import pearsonr, spearmanr
from nilearn.masking import apply_mask

groupA = sys.argv[1]
thresholds = [float(x) for x in sys.argv[2].split(",")]
standard = sys.argv[3]
output_pkl = sys.argv[4]
map_id = sys.argv[5]
# groupA_path,thresholds,standard,output_pkl,map_id

# standard brain mask that will be used for "brain mask" and to eliminate out of brain voxels for all masks
brain_mask = nib.load(standard)

# Read in main image
mrA = nib.load(groupA)
input_directory = os.path.dirname(groupA)
contrast_task = os.path.basename(groupA).replace("_groupA_tstat1.nii.gz","")

# Get degrees of freedom (it's 46)
try:
    dofA = len(open("%s/../copesA.txt" %(input_directory),"r").readlines())-2
except:
    dofA = 46

# Convert to Z score maps
mrA = IT.t_to_z(mrA,dofA)
  
# 3/12/2015: confirmed that thresholding at 0.0 is equivalent to original image
  
# Complete case analysis, single value imputation, multiple-imputation
cca_pearson = []
svi_pearson = []
cca_spearman = []
svi_spearman = []

# Edge cases:
# No surviving voxels: we append nan
# Thresholded image is 0: all correlations nan
# Fewer than 3 remaining voxels to compare: we append nan
# Suggested by Russ 3/13/2015, and I agree

# We also will save sizes of each, and a log to keep track of nan behavior
sizes = pandas.DataFrame(columns=["cca","svi"])
nanlog_cca = []
nanlog_svi = []

# Get group B images
groupB_maps = glob("%s/*groupB_tstat1.nii.gz" %(input_directory))

print "Calculating mask varieties [CCA,SVI] vs thresholded..."
idx = 0
size_ids = []
ids = []
thresh_id = []

# Human readable labels for absolute value
#  absolute_value == True means we take positive and negative
#  absolute_value == False means positive only

# Get dof for B maps
try:
    dofB = len(open("%s/../copesB.txt" %(input_directory),"r").readlines())-2
except: 
    dofB = 46

for groupB in groupB_maps:
    mrB = nib.load(groupB)
    mrB = IT.t_to_z(mrB,dofB)    
    contrast_task_B = os.path.basename(groupB).replace("_groupB_tstat1.nii.gz","")
    for thresh in thresholds:
	    print "Processing threshold %s" %(thresh)
	    for absolute_value in ["pos","posneg"]:
		if absolute_value == "posneg": 
		    # Group A is always unthresholded, B is thresholded
		    mrthresh = IT.threshold_abs(mrB,thresholds=[thresh])[thresh]
		else: 
		    mrthresh = IT.threshold_pos(mrB,thresholds=[thresh])[thresh]
                ids.append(contrast_task_B)
                thresh_id.append(thresh)
		# 3/12/2015: confirmed that first returns +/- values, second returns only positive  
		# If the image is empty thresholded, we append NaN
		if len(np.unique(mrthresh.get_data()))==1:
		    cca_pearson.append(np.nan)
		    svi_pearson.append(np.nan)
		    cca_spearman.append(np.nan)
		    svi_spearman.append(np.nan)
		    sizes.loc[idx] = [0,0]
		    nanlog_cca.append("nan_mrthresh_empty")
		    nanlog_svi.append("nan_mrthresh_empty")
		    size_ids.append("%s_%s" %(thresh,absolute_value))
		else:
		    # Generate a union (svi) and intersection (cca) mask
		    ccamask = IT.get_pairwise_deletion_mask(mrA,mrthresh,brain_mask)
		    svimask = IT.get_pairwise_inclusion_mask(mrA,mrthresh,brain_mask)
		    # COMPLETE CASE ANALYSIS (OLD PAIRWISE DELETION, intersection)
		    # Calculate correlation if there is overlap
		    if len(np.unique(ccamask.get_data())) == 2:
		       datacca = apply_mask([mrA,mrthresh],ccamask) 
		       # We need at least 3 values
		       if np.shape(datacca)[1] > 2: 
		           cca_pearson.append(pearsonr(datacca[0],datacca[1])[0])
		           cca_spearman.append(spearmanr(datacca[0],datacca[1])[0])
		           nanlog_cca.append("success")
		       else: 
		           cca_pearson.append(np.nan)
		           cca_spearman.append(np.nan)
		           nanlog_cca.append("nan_fewer_3_values")
		    # Otherwise (no overlap) it is nan
		    else: 
		        cca_pearson.append(np.nan)
		        cca_spearman.append(np.nan)
		        nanlog_cca.append("nan_no_overlap")    
		    # SINGLE VALUE IMPUTATION (old pairwise inclusion, union)
		    # Calculate correlation if there is overlap
		    if len(np.unique(svimask.get_data())) == 2:
		       datasvi = apply_mask([mrA,mrthresh],svimask) 
		       # We need at least 3 values
		       if np.shape(datasvi)[1] > 2: 
		           svi_pearson.append(pearsonr(datasvi[0],datasvi[1])[0])
		           svi_spearman.append(spearmanr(datasvi[0],datasvi[1])[0])
		           nanlog_svi.append("success")
		       else: 
		           svi_pearson.append(np.nan)
		           svi_spearman.append(np.nan)
		           nanlog_svi.append("nan_fewer_3_values")
		    # Otherwise (no overlap) it is nan
		    else: 
		        svi_pearson.append(np.nan)
		        svi_spearman.append(np.nan)
		        nanlog_svi.append("nan_no_overlap")
		    # Save sizes of all masks
		    sizes.loc[idx] = [len(datacca[0]),len(datasvi[0])]
		    size_ids.append("%s_%s" %(thresh,absolute_value))
	        idx+=1

# Save all data to output dictionary
output = {"id":map_id,
          "idB":ids,
          "thresh":thresh_id,
          "cca_pearson":cca_pearson,
          "svi_pearson":svi_pearson,
          "cca_spearman":cca_spearman,
          "svi_spearman":svi_spearman,
          "sizes":sizes,"size_ids":size_ids,
          "nanlog_cca":nanlog_cca,
          "nanlog_svi":nanlog_svi,
          "mrA_dof":dofA,
          "mrB_dof":dofB}

pickle.dump(output,open(output_pkl,"wb"))
