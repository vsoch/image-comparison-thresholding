# Here we can produce a similarity matrix to compare images A and B

import pandas
import nibabel
from pybraincompare.mr.datasets import get_standard_mask
from nilearn.masking import apply_mask
from scipi.stats import pearsonr

standard = get_standard_mask()

from glob import glob

amaps = glob("*groupA*.nii.gz")
bmaps = glob("*groupB*.nii.gz")

# Compare all A maps vs A maps
asim = pandas.DataFrame(index=amaps,columns=amaps)
bsim = pandas.DataFrame(index=bmaps,columns=bmaps)

for i in range(0,len(amaps)):
    print "Processing %s of %s" %(i,len(amaps))
    for j in range(0,len(amaps)):
        amapi = nibabel.load(amaps[i])
        amapj = nibabel.load(amaps[j])
        bmapi = nibabel.load(bmaps[i])
        bmapj = nibabel.load(bmaps[j])
        vectorai = apply_mask(amapi,standard)
        vectoraj = apply_mask(amapj,standard)
        vectorbi = apply_mask(bmapi,standard)
        vectorbj = apply_mask(bmapj,standard)
        asim.loc[amaps[i],amaps[j]] = pearsonr(vectorai,vectoraj)[0]
        bsim.loc[bmaps[i],bmaps[j]] = pearsonr(vectorbi,vectorbj)[0]

asim.to_csv("Amaps_unthresh_pearson_similarity.tsv",sep="\t")
bsim.to_csv("Bmaps_unthresh_pearson_similarity.tsv",sep="\t")
