# Code for paper "Hybrid local and distributed coding in PMd/M1 provides separation and interaction of bilateral arm signals"

## This repository contains code for analysis and visualization, organized by figures, for the paper "Hybrid local and distributed coding in PMd/M1 provides separation and interaction of bilateral arm signals". Some simple plotting code is omitted. Some code requires iterating over sessions, which has not been included as a batch script.

### Tanner Dixon, October 3, 2020


#### Fig 1: RT's, Reach durations, speed profiles

Code: 
create_rt_speed_mat.m
plot_speed_profiles.m

Data: 
trials


#### Fig 3: Firing rate traces

Code:
plot_unit_fr_targets.m

Data: 
unit_data


#### Fig 4: Arm pref and modulation ECDF's

Code:
calc_limb_dedication.m
plot_arm_pref.m

Data:
unit_data


#### Fig 5: Arm pref vs Mod slopes and modulation accounted for

Code:
calc_limb_dedication.m
plot_arm_pref.m
ap_vs_mod_slope_stats.m
dedicated_mod_stats.m

Data:
unit_data


#### Fig 6: PCA dimensionality and covariance alignment

Code:
prep_pca.m
estimate_dim_pca.m
pca_align_bihem.m
pca_align_bihem_epochs.m
crossval_pca_align.m

Data:
unit_data


#### Fig 7: PCA coeff analysis

Code:
calc_limb_dedication.m
crossval_pca_var_and_coeffs.m
plot_pca_coeff_analysis.m

Data:
unit_data


#### Fig 8: Pref/Non-pref analysis

Code:
calc_limb_dedication.m
pref_nonpref_mod_across_time.m
plot_pref_nonpref_mod_across_time.m
prep_lda.m
prep_lda_epochs.m
pp_lda.m
pp_lda_epochs.m
prep_pca_pref_epochs.m
pca_align_pref_epochs.m

Data:
unit_data
