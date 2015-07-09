# Here we will calculate thresholded (with random field theory) maps to compare to a standard threshold of 0.0 to 2.0

from clusterhcp.stats import threshold, get_clusters
from glob import glob
import pandas
import image_transformations as IT
from nilearn.masking import apply_mask
import nibabel
import os

data_directory = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/experiment3/0_permutations"

input_files = glob("%s/*.nii.gz" %data_directory)
mask = "/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain_mask.nii.gz"
output_folder = "/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/experiment3/threshold/random_field_theory"

# These are t maps, and default threshold is t=2.3 p=0.05
for input_image in input_files:
    output_prefix = input_image.split("/")[-1].split(".")[0]
    dlh, volume = threshold(input_image, output_folder, mask=mask)
    clusters = get_clusters(thresh_map=input_image,
                            dlh=dlh,volume=volume,
                            output_directory=output_folder,
                            output_prefix=output_prefix)

# Now we will count the number of voxels and compare the two methods

grf_files = glob("%s/*.nii.gz" %output_folder)
grf_files.sort()
input_files.sort()
voxel_counts = pandas.DataFrame(columns=["rft","thresh_0.0","thresh_1.0","thresh_2.0"])

brain_mask = nibabel.load(mask)
brain_voxels = len(brain_mask.get_data()[brain_mask.get_data()!=0])

for f in range(0,len(grf_files)):
  f1 = IT.t_to_z(nibabel.load(input_files[f]),46)
  f2 = nibabel.load(grf_files[f])
  print "Comparing %s vs %s" %(os.path.basename(grf_files[f]),os.path.basename(input_files[f]))
  mrthresh = IT.threshold_abs(f1,thresholds=[0.0,1.0,2.0])
  masked_mrthresh = [(thresh,apply_mask([nii],brain_mask)) for thresh,nii in mrthresh.iteritems()]
  # Count voxels in each
  image_name = os.path.basename(input_files[f])
  voxel_counts.loc[image_name] = (len(f2.get_data()[f2.get_data()!=0]),
                                 len(masked_mrthresh[0][1][0][masked_mrthresh[0][1][0]!=0]),
                                 len(masked_mrthresh[1][1][0][masked_mrthresh[1][1][0]!=0]),
                                 len(masked_mrthresh[2][1][0][masked_mrthresh[2][1][0]!=0]))
  
voxel_counts["absolute_difference_0.0"] = voxel_counts["thresh_0.0"] - voxel_counts["rft"]
voxel_counts["absolute_difference_1.0"] = voxel_counts["thresh_1.0"] - voxel_counts["rft"]
voxel_counts["absolute_difference_2.0"] = voxel_counts["thresh_2.0"] - voxel_counts["rft"]
voxel_counts["perc_brainvox_0.0"] = voxel_counts["thresh_0.0"] / brain_voxels
voxel_counts["perc_brainvox_1.0"] = voxel_counts["thresh_1.0"] / brain_voxels
voxel_counts["perc_brainvox_2.0"] = voxel_counts["thresh_2.0"] / brain_voxels
voxel_counts["perc_brainvox_rft"] = voxel_counts["rft"] / brain_voxels
voxel_counts["perc_difference_0.0"] = voxel_counts["absolute_difference_0.0"] / brain_voxels
voxel_counts["perc_difference_1.0"] = voxel_counts["absolute_difference_1.0"] / brain_voxels
voxel_counts["perc_difference_2.0"] = voxel_counts["absolute_difference_2.0"] / brain_voxels

voxel_counts.to_csv("/home/vanessa/Documents/Work/BRAINMETA/IMAGE_COMPARISON/experiment3/voxel_counts_rft_vs_thresh.tsv",sep="\t")

# Estimate accuracy with simple interpolation between 2.0 and 3.0
acc1 = 0.9842199809
acc2 = 0.976882272 # cca pearson, 2.0
acc3 = 0.9394306721 # cca pearson, 3.0
pertenth = (acc2 - acc3) / 10 # accuracy lost for each 10th
(pertenth*3)                  # accuracy lost going from 2.0 to 2.3
acc23 = acc2 - pertenth*3
loss_in_acc = acc1 - acc23
# 0.018573188869999924


