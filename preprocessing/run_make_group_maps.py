#!/usr/bin/python

# -------------------------------------------------------------
# Handling-missingness:  generate group maps

import pandas

# Variables to run jobs
basedir = os.path.abspath("../")
contrasts = pandas.read_csv("%s/doc/hcp_contrasts_id.tsv" %(basedir),sep="\t")
groups = pandas.read_csv("%s/doc/hcp_unrelated_groups.txt" %(basedir),sep="\t",header=None) # THIS FILE IS NOT PROVIDED
group_list_directory = "%s/data/copes" %(basedir)

# Parse into dict
groupdict = dict()
for g in groups.iterrows():
    groupdict[g[1][0]] = g[1][1].split(",")

# Now for each group, for each contrast, we will read in the image file, select the subjects, and run randomise to generate the group maps
for grp,subs in groupdict.iteritems():
  print "Processing %s" %(grp)
  for con in contrasts.iterrows():
    contrast_name = con[1].id
    contrast_directory = "%s/%s" %(output_directory,contrast_name)
    if not os.path.exists(contrast_directory): os.mkdir(contrast_directory)
    ss = ["%s" %(s) for s in subs]
    ss = ",".join(ss)
    contrast_list = "%s/%s_%s_copes.txt" %(group_list_directory,con[1].task,con[1].contrasts)
    contrast_map = "%s/%s_%s_copes_4D.nii.gz" %(group_maps_directory,con[1].task,con[1].contrasts)
    output_nii = "%s/%s_%s_4D.nii.gz" %(contrast_directory,grp,contrast_name)    
    # Write job to file
    filey = ".job/mk_%s_%s.job" %(grp,contrast_name)
    filey = open(filey,"w")
    filey.writelines("#!/bin/bash\n")
    filey.writelines("#SBATCH --job-name=mk_%s_%s\n" %(grp,contrast_name))
    filey.writelines("#SBATCH --output=.out/mk_%s_%s.out\n" %(grp,contrast_name))
    filey.writelines("#SBATCH --error=.out/mk_%s_%s.err\n" %(grp,contrast_name))
    filey.writelines("#SBATCH --time=2-00:00\n")
    filey.writelines("#SBATCH --mem=64000\n")
    filey.writelines("python make_group_maps.py %s %s %s %s %s %s\n" %(grp,contrast_name,output_nii,contrast_map,contrast_list,ss))
    filey.close()
    os.system("sbatch -p russpold .job/mk_%s_%s.job" %(grp,contrast_name))
