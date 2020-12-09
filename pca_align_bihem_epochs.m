function [X, princ_align_self, princ_align_lr]...
    = pca_align_bihem_epochs(unit_data, configs, p)

%
% Computes cross-validated alignment of p-dimensional PCA subspaces within 
% discrete task phases.
%
%
% INPUTS: 
%
% unit_data - unit-separated data struct containing hemisphere, brain area,
%             firing rate and other metrics
%             (struct: 1 x num_units)
%
% configs - selected starting configurations of the hands to use. 1:3
%
% p - number of principle components to use. integer
%
%
% OUTPUTS:
%
% X - for each hand, time-locked and trial-separated firing rate traces for
%     all simultaneously recorded units
%     (struct: 1x1, fields 'l_hand' and 'r_hand'
%      nested cell arrays: num_trials x 1
%      nested matrices: num_samples x num_units)
%
% princ_align_self - alignment across time for trials using the same hand
%                    (struct: 1x1, fields 'rest', 'prep', 'move'
%                     nested struct: 1x1, fields 'l_hand' and 'r_hand'.
%                     nested vectors: 1 x 1000, MCCV alignment values)
% 
% princ_align_self - alignment across time for trials using the same hand
%                    (struct: 1x1, fields 'rest', 'prep', 'move'
%                     nested struct: 1x1, fields 'l_on_r' and 'r_on_l'.
%                     nested vectors: 1 x 1000, MCCV alignment values)
% 


%% Setup

X = prep_pca(unit_data, configs, 'both', 'rest');
X_l_hand_total = cat(3,X.l_hand{:});
X_l_hand.rest = permute(X_l_hand_total(11:26,:,:), [3,2,1]);
X_l_hand.prep = permute(X_l_hand_total(37:52,:,:), [3,2,1]);
X_l_hand.move = permute(X_l_hand_total(63:78,:,:), [3,2,1]);
X_r_hand_total = cat(3,X.r_hand{:});
X_r_hand.rest = permute(X_r_hand_total(11:26,:,:), [3,2,1]);
X_r_hand.prep = permute(X_r_hand_total(37:52,:,:), [3,2,1]);
X_r_hand.move = permute(X_r_hand_total(63:78,:,:), [3,2,1]);


%% Calculate PC alignment between and within hands for each epoch

[princ_align_self.rest.l_hand,~,~] = crossval_pca_align(...
    X_l_hand.rest, X_l_hand.rest, 1, p, 1000);
[princ_align_self.rest.r_hand,~,~] = crossval_pca_align(...
    X_r_hand.rest, X_r_hand.rest, 1, p, 1000);
[princ_align_lr.rest.l_on_r,~,~] = crossval_pca_align(...
    X_l_hand.rest, X_r_hand.rest, 0, p, 1000);
[princ_align_lr.rest.r_on_l,~,~] = crossval_pca_align(...
    X_r_hand.rest, X_l_hand.rest, 0, p, 1000);

[princ_align_self.prep.l_hand,~,~] = crossval_pca_align(...
    X_l_hand.prep, X_l_hand.prep, 1, p, 1000);
[princ_align_self.prep.r_hand,~,~] = crossval_pca_align(...
    X_r_hand.prep, X_r_hand.prep, 1, p, 1000);
[princ_align_lr.prep.l_on_r,~,~] = crossval_pca_align(...
    X_l_hand.prep, X_r_hand.prep, 0, p, 1000);
[princ_align_lr.prep.r_on_l,~,~] = crossval_pca_align(...
    X_r_hand.prep, X_l_hand.prep, 0, p, 1000);

[princ_align_self.move.l_hand,~,~] = crossval_pca_align(...
    X_l_hand.move, X_l_hand.move, 1, p, 1000);
[princ_align_self.move.r_hand,~,~] = crossval_pca_align(...
    X_r_hand.move, X_r_hand.move, 1, p, 1000);
[princ_align_lr.move.l_on_r,~,~] = crossval_pca_align(...
    X_l_hand.move, X_r_hand.move, 0, p, 1000);
[princ_align_lr.move.r_on_l,~,~] = crossval_pca_align(...
    X_r_hand.move, X_l_hand.move, 0, p, 1000);

disp('Completed')

