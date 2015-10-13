### Effects of thresholding on correlation-based image similarity metrics

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.32074.svg)](http://dx.doi.org/10.5281/zenodo.32074)

This repository includes scripts to generate initial subsampling data, and then carry out analyses. Since generation of the samples data requires considerable computational resources and access to private subject data to indicate relatedness (NOT PROVIDED) we are including the output of these methods (with no proprietary data) on which subsequent analyses can be carried out. You can visualize the confusion matrices for our result [here](http://vsoch.github.io/image-comparison-thresholding). The structure of this repository is as follows:

#### preprocessing
Includes scripts for working with raw data to run permutations to generate image comparison scores across thresholds, similarity metrics, directions (positive vs. negative values), and groups of unrelated subjects. Meta data about subjects is not provided, however a researcher looking to reproduce analyses in entirely can [apply for access](http://www.humanconnectome.org/data/data-use-terms/restricted-access-overview.html).

###### make_contrast_paths.py 	
This produces a log file that lists all of the HCP contrasts, from the raw data.  Generation of this file requires the entire HCP dataset, and so the file is provided under docs ([hcp_contrasts_id.tsv](doc/hcp_contrasts_id_filter.tsv)). The script also generates a *_copes.txt file for each cope, which lists the exact paths for each subject to the cope file (for use to generate group maps during the permutations). We are [providing these cope files as an example](data/copes), and obviously these would need to be regenerated for complete paths on a cluster-specific implementation. You will need to modify the first path in this file to point to the complete HCP single subject datasets.

###### filter_contrasts.py
A set of non-redundant contrasts were selected (eg, faces-shapes would be equivalent to shapes-faces), and this file filters the original set of 86 (hcp_contrasts_id) down to a non-redundant set (47). The filtered file [is provided](doc/hcp_contrasts_id_filter.tsv)

###### summarize_hcp.py
This identifies a subset of HCP subjects that have all contrasts for all tasks (465). The output of this script [is provided](doc/hcp_465_with_all_tasks.tsv)

###### image_transformations.py
Provides functions for thresholding images (including positive and negative values, and positive values only), as well as generate masks that represent different strategies for handling missing values in the images (complete case analysis and single value imputation), as well as a function to convert from T images to Z score maps. This function has since been released as a [standalone package](https://github.com/vsoch/TtoZ) for broader use.

###### make_group_maps.py
Since it was computationally feasible, we generated the group maps (for each of a group A and B) for each of the 500 permutations in advance, that way a comparison within a permutation that relied upon a map for each of group A and B could have assurance that the maps already existed.  This script was run via [run_make_group_maps.py](preprocessing/run_make_group_maps.py) and relied upon the [clusterhcp](https://github.com/vsoch/clusterhcp) package. This is the step that requires knowing relatedness of individuals in HCP, specifically groups of unrelated subjects in a file called "hcp_unrelated_groups.txt." While we provide the methods for running the map generation, this data file is not provided in this repository or in cluster hcp. To define lists of unrelated subjects to generate unrelated maps, you will have to pursue acquiring the restricted data on your own.

###### test_threshold.py
This script is intended to be submit by "run_test_threshold.py" in a cluster environment to generate permuted results of similarity scores using different metrics (Pearson and Spearman) across a range of thresholds (+/- and + 0.0 to 13.0) and two comparison strategies (complete case analysis and single value imputation) using brain maps from the subset of subjects identified as having all contrasts, across all 47 task contrasts. The output of this is, for each permutation, a folder named by the permutation number, with two subdirectories, "maps" and "permutations." The maps directory contains maps for groups A and B for each task contrast. The permutations directory contains a pickle file for each contrast that contains the similarity metrics described above. We have [provided an example](data/permutation_example/0). We have removed two files in the top directory, "copesA.txt" and "copesB.txt" that indicated groups of unrelated subjects. The script [run_test_threshold.py](preprocessing/run_test_threshold.py) will submit instances of this script to generate all permutation data for analyses. You will need to change one of the first lines to indicate the path to the single subject raw HCP data.

###### compile_threshold_results.py
This script iterates over the permuation directories described above to generate a set of compiled data files for each permutation. Since we did not provide the entire set of permutations datasets and maps, we are providing these [compiled results files](data/results).


#### analyses
The final analysis is done in R, which includes combining all of the permutation text files described above into R data objects. 

###### threshold_analysis.R
To do this, the main analysis pipeline is [threshold_analysis.R](analysis/threshold_analysis.R) and helper functions are provided in [helper_functions.R](analysis/helper_functions.R). Images and figures output for the paper are generated in this file, and a subset of images with a "command.txt" script are [also provided](data/img) to demonstrate generation of an animation from static images. 

###### average_across_samples.py
To generate our Table 1 (and the confusion matrices for each thresholding, comparison strategy, and similarity metric) we needed to take an average across samples, accomplished with this script.

###### compare_to_tft.py
For our discussion, we wanted to compare the percentage of voxels obtained at our optimal threshold to maps generated with random field theory thresholding, and we did this with [compare_to_rft.py](analysis/compare_to_rft.py) 

#### results
Complete code to generate the results interface, an interactive confusion matrix for the classification task, is also provided as a [flask interface](results) or [static site](https://github.com/vsoch/image-comparison-thresholding/tree/gh-pages). A live version is also [available](http://vsoch.github.io/image-comparison-thresholding).
