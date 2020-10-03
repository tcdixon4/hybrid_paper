function [princ_var] = pca_align_pref_epochs(unit_data, p, norm_method)

%
% Computes PCA models using trials of the "non-preferred" arm and projects
% held-out activity for both arms through these models. Then, compares the
% amount of variance preserved through the projection for the hand the
% model was trained on (non-pref) vs the hand it was not (pref). For
% analyzing the "distributed subspace".
%
%
% INPUTS: 
%
% unit_data - unit-separated data struct containing hemisphere, brain area,
%             firing rate and other metrics
%             (struct: 1 x num_units)
%
% p - selected number of PCA components to use. integer
%
% norm_method - selected normalization method for unit firing rates. 
%               'rest' for resting std + 1hz
%               'range' for full firing rate range + 5hz
%
%
% OUTPUTS:
%
% princ_var - variance (raw value) preserved through PCA projections
%             (struct: 1x1, fields 'l_pref' and 'r_pref'
%              nested struct: 1x1, fields 'l_hand' and 'r_hand'
%              nested struct: 1x1, fields 'rest', 'prep', and 'move'
%              nested scalar: variance)
%
% 


%% Setup

[X_ecc, X_cen] = prep_pca_pref_epochs(unit_data, norm_method);


%% train pca models on center config data using only non-preferred arm
% trials and project data for both arms, then calculate variance captured

%% left preferring units
X_train = X_cen.l_pref.r_hand.rest;
X_test_l_hand = X_ecc.l_pref.l_hand.rest;
X_test_r_hand = X_ecc.l_pref.r_hand.rest;
[P,~,~,~,~,~] = ...
    pca(X_train,'NumComponents',p, 'Centered',false);
T_l_hand = X_test_l_hand*P;
T_r_hand = X_test_r_hand*P;
princ_var.l_pref.l_hand.rest = trace(cov(T_l_hand));
princ_var.l_pref.r_hand.rest = trace(cov(T_r_hand));


X_train = X_cen.l_pref.r_hand.prep;
X_test_l_hand = X_ecc.l_pref.l_hand.prep;
X_test_r_hand = X_ecc.l_pref.r_hand.prep;
[P,~,~,~,~,~] = ...
    pca(X_train,'NumComponents',p, 'Centered',false);
T_l_hand = X_test_l_hand*P;
T_r_hand = X_test_r_hand*P;
princ_var.l_pref.l_hand.prep = trace(cov(T_l_hand));
princ_var.l_pref.r_hand.prep = trace(cov(T_r_hand));


X_train = X_cen.l_pref.r_hand.move;
X_test_l_hand = X_ecc.l_pref.l_hand.move;
X_test_r_hand = X_ecc.l_pref.r_hand.move;
[P,~,~,~,~,~] = ...
    pca(X_train,'NumComponents',p, 'Centered',false);
T_l_hand = X_test_l_hand*P;
T_r_hand = X_test_r_hand*P;
princ_var.l_pref.l_hand.move = trace(cov(T_l_hand));
princ_var.l_pref.r_hand.move = trace(cov(T_r_hand));


%% right preferring units
X_train = X_cen.r_pref.r_hand.rest;
X_test_l_hand = X_ecc.r_pref.l_hand.rest;
X_test_r_hand = X_ecc.r_pref.r_hand.rest;
[P,~,~,~,~,~] = ...
    pca(X_train,'NumComponents',p, 'Centered',false);
T_l_hand = X_test_l_hand*P;
T_r_hand = X_test_r_hand*P;
princ_var.r_pref.l_hand.rest = trace(cov(T_l_hand));
princ_var.r_pref.r_hand.rest = trace(cov(T_r_hand));


X_train = X_cen.r_pref.r_hand.prep;
X_test_l_hand = X_ecc.r_pref.l_hand.prep;
X_test_r_hand = X_ecc.r_pref.r_hand.prep;
[P,~,~,~,~,~] = ...
    pca(X_train,'NumComponents',p, 'Centered',false);
T_l_hand = X_test_l_hand*P;
T_r_hand = X_test_r_hand*P;
princ_var.r_pref.l_hand.prep = trace(cov(T_l_hand));
princ_var.r_pref.r_hand.prep = trace(cov(T_r_hand));


X_train = X_cen.r_pref.r_hand.move;
X_test_l_hand = X_ecc.r_pref.l_hand.move;
X_test_r_hand = X_ecc.r_pref.r_hand.move;
[P,~,~,~,~,~] = ...
    pca(X_train,'NumComponents',p, 'Centered',false);
T_l_hand = X_test_l_hand*P;
T_r_hand = X_test_r_hand*P;
princ_var.r_pref.l_hand.move = trace(cov(T_l_hand));
princ_var.r_pref.r_hand.move = trace(cov(T_r_hand));






    
    