This folder includes complete results for all comparison methods (complete case analysis and single value imputation) across all thresholds and directions, for a single permutation. For the raw data, please see the [maps](../maps) folder. To load the result for a particular contrast:

      git clone https://github.com/vsoch/image-comparison-thresholding/
      cd image-comparison-thresholding/data/example_permutation/0/comparisons
      python

      import pickle
      data = pickle.load(open("TASK01_CON07.pkl","rb"))

      data.keys()
      ['nanlog_cca',      # includes ['nan_fewer_3_values', 'nan_mrthresh_empty', 'nan_no_overlap', 'success']
       'svi_spearman',    # single value imputation, spearman
       'svi_pearson',     # single value imputation, pearson
       'sizes',           # number of voxels compared
       'cca_pearson',     # complete case analysis, pearson
       'mrA_dof',         # degrees of freedom, map A
       'nanlog_svi',      
       'mrB_dof',         # degrees of freedom, map B
       'thresh',          # threshold for image B
       'size_ids',        # ids that correspond to "sizes" eg, ['0.0_pos'.. '0.0_posneg'...]
       'cca_spearman',    # complete case analysis, spearman
       'idB',             # id of image B
       'id']              # id if image A

