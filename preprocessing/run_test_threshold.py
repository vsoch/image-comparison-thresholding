#!/usr/bin/python

# ---------------------------------------------------------------------------------
# handling-missingness: run test threshold instance

# This batch script will prepare and submit jobs for running on a SLURM cluster
# This experiment will prepare a flat database of similarity metrics assessed over a set of maps
# for different thresholds, to be used for higher level analysis

import os
import pickle
import time
import pandas
import numpy
import nibabel
from glob import glob
from clusterhcp.database import get_hcp_paths
from clusterhcp.stats import select_unrelated_samples

# This is the top directory of the HCP data
top_directory = "/path/to/HCP"
basedir = os.path.abspath("../")
outdirectory = "%s/data/permutations" %(basedir)

# This is the size of the groups that we want to generate
size = 46

# Number of permutation runs
nruns = 500

# Read in contrast list
input_file = "%s/doc/hcp_contrasts_id_filter.csv" %(basedir)

# Read in input file
contrasts = pandas.read_csv(input_file,sep="\t")

### STEP 0: Try loading each file
for con in contrasts.iterrows():     
    task = con[1]["task"]
    contrast = con[1]["contrasts"]
    map_id = con[1]["id"]
    print "Testing %s" %(map_id)
    paths = get_hcp_paths(top_directory, tasks=task, contrasts=contrast)
    for p in paths:
        try:
            tmp = nibabel.load(p)
        except:
            print "ERROR: problem with image %s" %(p)


### STEP 1: First we will define groups A and B for 500 runs ############################
nruns=500
# Generate the subject groups first, and then run them over iterations
for i in range(0,nruns):
    # Top level of output directory is for the iteration
    output_directory = "%s/%s" %(outdirectory,i)
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
        os.mkdir("%s/maps" %(output_directory))
        os.mkdir("%s/comparisons" %(output_directory))
    maps_directory = "%s/maps" %(output_directory)
    # Get a group of subjects for one of the tasks/contrasts
    task = 'tfMRI_WM'
    contrast = '0BK_BODY'
    paths = get_hcp_paths(top_directory, tasks=task, contrasts=contrast) # default is to return 465 subjects with all tasks
    # get unrelated samples, default is two groups
    A,B = select_unrelated_samples(paths,size=size)
    subjectsA = [int(s.split("/")[7]) for s in A]
    subjectsB = [int(s.split("/")[7]) for s in B]
    outputA = "\n".join([str(x) for x in subjectsA])
    outputB = "\n".join([str(x) for x in subjectsB])
    fileyA = open("%s/copesA.txt" %(output_directory),"wb")
    fileyB = open("%s/copesB.txt" %(output_directory),"wb") 
    fileyA.writelines(outputA)
    fileyB.writelines(outputB)
    fileyA.close()
    fileyB.close()
    

### STEP 2: Process groups, A and B ############################
for con in contrasts.iterrows(): 
    print "Processing contrast %s" %(con[1].id)  
    task = con[1]["task"]
    contrast = con[1]["contrasts"]
    map_id = con[1]["id"]
    paths = get_hcp_paths(top_directory, tasks=task, contrasts=contrast)
    for i in range(0,nruns):
        # Top level of output directory is for the iteration
        output_directory = "%s/%s" %(outdirectory,i)
        maps_directory = "%s/maps" %(output_directory)
        # Load the groups
        subjectsA = open("%s/copesA.txt" %(output_directory),"r").readlines()
        subjectsB = open("%s/copesB.txt" %(output_directory),"r").readlines()
        subjectsA = [x.replace("\n","") for x in subjectsA]
        subjectsB = [x.replace("\n","") for x in subjectsB]
        A = [s for s in paths if s.split("/")[7] in subjectsA]
        B = [s for s in paths if s.split("/")[7] in subjectsB]
        groupA = ",".join(A)
        groupB = ",".join(B)
        filey = ".job/%s_%s.job" %(i,map_id)
        filey = open(filey,"w")
        filey.writelines("#!/bin/bash\n")
        filey.writelines("#SBATCH --job-name=%s_%s\n" %(i,map_id))
        filey.writelines("#SBATCH --output=.out/%s_%s.out\n" %(i,map_id))
        filey.writelines("#SBATCH --error=.out/%s_%s.err\n" %(i,map_id))
        filey.writelines("#SBATCH --time=1:00\n")
        filey.writelines("module load fsl\n")
        filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_maps_tacc.py %s %s %s %s\n" %(groupA,groupB,maps_directory,map_id)) 
        filey.close()
        os.system("sbatch -p russpold .job/%s_%s.job" %(i,map_id))


### STEP 3: Find any missing maps ###########################################################
for con in contrasts.iterrows(): 
    print "Processing contrast %s" %(con[1].id)  
    task = con[1]["task"]
    contrast = con[1]["contrasts"]
    map_id = con[1]["id"]
    paths = get_hcp_paths(top_directory, tasks=task, contrasts=contrast)
    for i in range(0,nruns):
        # Top level of output directory is for the iteration
        output_directory = "%s/%s" %(outdirectory,i)
        maps_directory = "%s/maps" %(output_directory)
        outfileA = "%s/%s_groupA_tstat1.nii.gz" %(maps_directory,map_id)
        outfileB = "%s/%s_groupB_tstat1.nii.gz" %(maps_directory,map_id)
        if not os.path.exists(outfileA):
            print "Missing %s" %(outfileA)        
            subjectsA = open("%s/copesA.txt" %(output_directory),"r").readlines()
            subjectsA = [x.replace("\n","") for x in subjectsA]
            A = [s for s in paths if s.split("/")[7] in subjectsA]
            groupA = ",".join(A)
            filey = ".job/%s_A_%s.job" %(i,map_id)
            filey = open(filey,"w")
            filey.writelines("#!/bin/bash\n")
            filey.writelines("#SBATCH --job-name=%s_A_%s\n" %(i,map_id))
            filey.writelines("#SBATCH --output=.out/%s_A_%s.out\n" %(i,map_id))
            filey.writelines("#SBATCH --error=.out/%s_A_%s.err\n" %(i,map_id))
            filey.writelines("#SBATCH --time=1:00\n")
            filey.writelines("module load fsl\n")
            filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_map_single_tacc.py %s %s %s %s\n" %(groupA,"A",maps_directory,map_id)) 
            filey.close()
            os.system("sbatch -p russpold .job/%s_A_%s.job" %(i,map_id))
        if not os.path.exists(outfileB):
            print "Missing %s" %(outfileB)        
            subjectsB = open("%s/copesB.txt" %(output_directory),"r").readlines()
            subjectsB = [x.replace("\n","") for x in subjectsB]
            B = [s for s in paths if s.split("/")[7] in subjectsB]
            groupB = ",".join(B)
            filey = ".job/%s_B_%s.job" %(i,map_id)
            filey = open(filey,"w")
            filey.writelines("#!/bin/bash\n")
            filey.writelines("#SBATCH --job-name=%s_B_%s\n" %(i,map_id))
            filey.writelines("#SBATCH --output=.out/%s_B_%s.out\n" %(i,map_id))
            filey.writelines("#SBATCH --error=.out/%s_B_%s.err\n" %(i,map_id))
            filey.writelines("#SBATCH --time=1:00\n")
            filey.writelines("module load fsl\n")
            filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/make_group_map_single_tacc.py %s %s %s %s\n" %(groupB,"B",maps_directory,map_id))           
            filey.close()
            os.system("sbatch -p russpold .job/%s_B_%s.job" %(i,map_id))


### STEP 4: Which runs are done? ############################

done = []
count=0
for i in range(0,nruns):
   num = glob("/share/PI/russpold/work/IMAGE_COMPARISON/experiment3/permutations/%s/maps/*_tstat1.nii.gz" %(i))
   if len(num)==94: done.append(i)
   else:
       count = count + (94-len(num))

## Delete intermediate files to free up room

for r in range(0,500):
  print r
  os.system("rm /share/PI/russpold/work/IMAGE_COMPARISON/experiment3/permutations/%s/maps/*4D.nii.gz" %(r))


### STEP 5: Calculate similarities for maps that have all images generated ############################

# We will run over a set of thresholds
thresholds = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0]
thresholds = ",".join([str(x) for x in thresholds])

standard = "%s/standard/MNI152_T1_2mm_brain_mask.nii.gz" %(basedir)

for i in [231, 237, 253]: 
    print i   
    output_directory = "%s/%s" %(outdirectory,i)
    maps_directory = "%s/maps" %(output_directory)
    # There will be one of group A for each of group B, unrelated
    image_pairs = pandas.DataFrame()
    groupA_files = numpy.sort(glob("%s/*_groupA_tstat1.nii.gz" %(maps_directory)))
    groupB_files = numpy.sort(glob("%s/*_groupB_tstat1.nii.gz" %(maps_directory)))
    if len(groupA_files) == len(groupB_files):
	    image_pairs["groupA"] = groupA_files
	    image_pairs["groupB"] = groupB_files
	    for run in image_pairs.iterrows():
		mapnameA = run[1].groupA.replace("_groupA_tstat1.nii.gz","")
		mapnameB = run[1].groupB.replace("_groupB_tstat1.nii.gz","")
		if mapnameA == mapnameB:
		    groupA_path = run[1]["groupA"]
		    groupB_path = run[1]["groupB"]
		    contrast_task = groupA_path.split("/")[-1].replace("_groupA_tstat1.nii.gz","")
                    output_pkl = "%s/comparisons/%s.pkl" %(output_directory,contrast_task)
		    if not os.path.exists(output_pkl):
                        filey = ".job/%s_%s.job" %(i,contrast_task)
                        print filey
                        filey = open(filey,"w")
                        filey.writelines("#!/bin/bash\n")
                        filey.writelines("#SBATCH --job-name=%s_%s\n" %(i,contrast_task))
                        filey.writelines("#SBATCH --output=.out/%s_%s.out\n" %(i,contrast_task))
                        filey.writelines("#SBATCH --error=.out/%s_%s.err\n" %(i,contrast_task))
                        filey.writelines("#SBATCH --time=30:00\n")
                        filey.writelines("module load fsl\n")
		        filey.writelines("python /home/vsochat/SCRIPT/python/brainmeta/image_comparison/experiments/experiment3/test_thresholding_tacc.py %s %s %s %s %s\n" %(groupA_path,thresholds,standard,output_pkl,contrast_task))     
                        filey.close()   
                        os.system("sbatch -p russpold .job/%s_%s.job" %(i,contrast_task))
		else:
		    print "Error, mismatch for run %s %s" %(run[1].groupA,run[1].groupB)
    else:
        print "Error: iteration %s is missing maps, fix!" %(i)
