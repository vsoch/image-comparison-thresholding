#!/usr/bin/python

# -------------------------------------------------------------
# Handling-missingness:  generate group maps

import sys
from clusterhcp.stats import run_randomise_local

groupA = sys.argv[1].split(",")
groupB = sys.argv[2].split(",")
output_directory = sys.argv[3]
map_id = sys.argv[4]

# groupA,groupB,maps_directory,map_id

# Generate group maps with randomise, to output folder
run_randomise_local(groupA,output_directory,output_prefix="%s_groupA" %(map_id))
run_randomise_local(groupB,output_directory,output_prefix="%s_groupB" %(map_id))
