import nibabel
from nilearn.plotting import plot_stat_map
from glob import glob
import pandas
import os

outdir = "../img"
# These maps are not provided
input_maps = glob("../../data/results/0/*groupA_tstat1.nii.gz")
os.chdir("../preprocessing/")
import image_transformations as IT

# Load in dataframe to look up image and contrast names
lookup = pandas.read_csv("../doc/hcp_contrasts_id.tsv",sep="\t")

# Get max and min values
max = 0
min = 0
for input_map in input_maps:
    img = nibabel.load(input_map)
    img = IT.t_to_z(img,46).get_data()
    if img.max() > max:
        max = img.max()
    if img.min() < min:
        min = img.min()


for i in range(0,len(input_maps)):
    print "Processing %s of %s" %(i,len(input_maps))
    input_map = input_maps[i]
    mapname = os.path.basename(input_map).replace("_groupA_tstat1.nii.gz","")
    task = lookup.task[lookup.id==mapname].tolist()[0]
    contrast = lookup.contrasts[lookup.id==mapname].tolist()[0]
    title = "%s_%s" %(task,contrast)
    image = nibabel.load(input_map)
    Z = IT.t_to_z(image,46)
    for thresh in range(0,10):     
        print thresh
        pos = IT.threshold_pos(image,thresholds=[thresh])[thresh]
        posneg = IT.threshold_abs(image,thresholds=[thresh])[thresh]
        fig = plot_stat_map(pos,vmax=max,title="%s Positive Only" %title,cut_coords=(0,0,0))
        fig.savefig("%s/%s_%s_pos.png" %(outdir,thresh,mapname))
        fig.close()
        fig = plot_stat_map(posneg,vmax=max,title="%s Positive and Negative" %title,cut_coords=(0,0,0))
        fig.savefig("%s/%s_%s_posneg.png" %(outdir,thresh,mapname))
        fig.close()
