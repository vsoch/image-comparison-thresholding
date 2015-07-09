#!/usr/bin/python

# We have negative contrasts, and opposite contrasts, and we should eliminate them to not be redundant.
import os
import pandas
import pickle

# Input file with map paths
basedir = os.path.abspath("../")
input_file = "%s/doc/hcp_contrasts_id.tsv" %(basedir)
inputs = pandas.read_csv(input_file,sep="\t")

tokeep = ["TASK06_CON01","TASK06_CON02","TASK06_CON06","TASK01_CON07","TASK01_CON08","TASK01_CON09",
         "TASK05_CON13","TASK05_CON14","TASK05_CON15","TASK03_CON19","TASK03_CON20","TASK03_CON22",
         "TASK02_CON25","TASK02_CON26","TASK02_CON27","TASK07_CON31","TASK07_CON32","TASK07_CON33",
         "TASK07_CON34","TASK07_CON35","TASK07_CON36","TASK07_CON37","TASK07_CON38","TASK07_CON39",
         "TASK07_CON40","TASK07_CON41","TASK07_CON45","TASK07_CON46","TASK07_CON47","TASK07_CON48",
         "TASK07_CON49","TASK07_CON50","TASK07_CON51","TASK07_CON52","TASK04_CON61","TASK04_CON62",
         "TASK04_CON63","TASK04_CON64","TASK04_CON65","TASK04_CON66","TASK04_CON67","TASK04_CON68",
         "TASK04_CON69","TASK04_CON70","TASK04_CON71","TASK04_CON72","TASK04_CON73"]

inputs_filtered = pandas.DataFrame(columns=inputs.columns)
for keeper in tokeep:
  subset = inputs[inputs.id == keeper]
  inputs_filtered = inputs_filtered.append(subset)

# Save to file
output_file = "%s/doc/hcp_contrasts_id_filter.tsv" %(basedir)
inputs_filtered.to_csv(output_file,sep="\t")
