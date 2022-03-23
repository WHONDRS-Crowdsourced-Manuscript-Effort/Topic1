# Output thresholds to distinguish core and satellite

Abbreviations:

- MF = molecular formula
- CS = core-satellite

---

## Files in this folder:

- 'FTICR_crosstable_rep.merged1_all_em.thres_yyyy-mm-dd.csv'
  Cross table including all three thresholding approaches for the dataset where MF present in at least 1 replicate were kept ('rep.merged1').
- 'FTICR_crosstable_rep.merged2_all_em.thres_yyyy-mm-dd.csv'
  Cross table including all three thresholding approaches for the dataset where MF present in at least 2 replicates were kept ('rep.merged2').
- 'FTICR_peaks_crosstable_rep.merged1_all_em.thres_yyyy-mm-dd.csv'
  Cross table including only the emergent threshold approach for the **peak** dataset where peaks present in at least 1 replicate were kept ('rep.merged1').
- 'FTICR_peaks_crosstable_rep.merged2_all_em.thres_yyyy-mm-dd.csv'
  Cross table including only the emergent threshold approach for the **peak** dataset where peaks present in at least 2 replicates were kept ('rep.merged1').

>If there are several dates, please use the **newest version**.

---

## New columns are:

  - 'occupancy_sed': Occupancy in terms of the number of **sediment** sites that the particular MF is found in.
  - 'occupancy_water': Occupancy in terms of the number of **surface water** sites that the particular MF is found in.
  - 'perc.occup_sed': Occupancy in **sediments** relative to the total number of sediment sites, given in percent. (This variable was used for thresholding CS)
  - 'perc.occup_water': Occupancy in **surface waters** relative to the total number of surface water sites, given in percent. (This variable was used for thresholding CS)
  - 'cs.flag.emergent_sed': CS classification based on the second derivative of the frequency-occupancy distribution of sediment samples.
  - 'cs.flag.emergent_water': CS classification based on the second derivative of the frequency-occupancy distribution of surface water samples.
  - 'cs.flag.pca_sed': Literature thresholds compared and best threshold selected based on variance explained in PCA of sediment samples (only MF dataset).
  - 'cs.flag.pca_water': Literature thresholds compared and best threshold selected based on variance explained in PCA of surface water samples (only MF dataset).
  - 'cs.flag.rf_sed': CS classification based on random forest analysis of sediment samples (only MF dataset).
  - 'cs.flag.rf_water': CS classification based on random forest analysis of surface water samples (only MF dataset).
