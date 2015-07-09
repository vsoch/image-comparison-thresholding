#!/usr/bin/python

# This script will summarize the tasks that we have for HCP, and output how many subjects per, etc.
# We need to find a group of subjects that have all tasks.

import os
import sys
import nibabel
import pandas
import numpy as np
from glob import glob

# Directory with subject lists
top_directory = os.path.abspath("../")
copes_directory = "%s/data/copes" %(top_directory)
copes_files = glob("%s/*" %(copes_directory))

# How many copes?
len(copes_files)
# 86

# Contrasts file corresponding to the copes files
contrasts_file = pandas.read_csv("%s/doc/hcp_contrasts_id.tsv" %(top_directory),sep="\t")

# This is raw data for all tasks, resting, all subjects. We can get all subject ids from here
datalog = pandas.read_csv("%s/doc/hcp_datalog.tsv" %(top_directory),sep="\t")

# We will save a data frame to keep a lot of who has data for each task!
# Get subids from the complete task data
subids = np.unique([x.split("/")[7] for x in datalog.volume_folders])
# 501

df = pandas.DataFrame(columns=subids)
# Here we will save list of task contrasts
task_contrasts = []
count=1
# We now want to find a set of subjects that have data for all the tasks
for cope_file in copes_files:
  task_contrast = os.path.split(cope_file)[1].replace("_copes.txt","")
  task_contrasts.append(task_contrast)
  cope_file = pandas.read_csv(cope_file,sep="\t",header=None)
  # You may need to adjust this index depending on the path to your cope files
  subids = [x.split("/")[2] for x in cope_file[0]]
  df.loc[count,subids] = 1
  count = count+1

# Add the task_contrasts
df.index = task_contrasts
df.to_csv("%s/doc/hcp_ss_tasklog.tsv" %(top_directory),sep="\t")

# For which subjects do we have data for all tasks?
subids = df.sum(axis=0)==df.shape[0]
subids = subids.index[subids==True]
subids = subids.tolist()
len(subids)
hcp_subjects = pandas.DataFrame()
hcp_subjects["id"] = subids
hcp_subjects.to_csv("%s/doc/hcp_465_with_all_tasks.tsv" %(top_directory),sep="\t")
# 465
