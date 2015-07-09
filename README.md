### Steps to do experiment 3

- First summarize hcp data (summarize_hcp.py)
- Prepare meta data files and logs for analysis (prep_test_thresholding.py)
- Generate permuted group maps, and generate pickle results for spearman and pearson at different thresholds (run_test_threshold_slurm.py,test_threshold_tacc.py)
- Compile results (compile_threshold_results.py)
- Analysis in R (threshold_analysis.R, helper_functions.R)

### Extras
- find_missing.py will identify thresholds and directions to re-run
- filter_contrasts.py will remove "neg" contrasts, or inverse of the same (eg, faces-shapes == shapes-faces)
- compare_to_rft.py: will compare threshold result to random field theory
