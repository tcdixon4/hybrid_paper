function [X, princ_align_self, princ_align_lr] = ...
    pca_align_bihem(unit_data, configs, p)

%
% Computes cross-validated alignment of p-dimensional PCA subspaces on a
% fine timescale, taking the alignment between all pairwise comparisons of
% timepoints.
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
%                    (struct: 1x1, fields 'l_hand' and 'r_hand'.
%                     nested matrices: num_samples x num_samples, mean  
%                     alignment values for each pair of timepoints)
% 
% princ_align_lr - alignment across time for trials using different hands
%                  (struct: 1x1, fields 'l_on_r' and 'r_on_l'.
%                   nested matrices: num_samples x num_samples, mean  
%                   alignment values for each pair of timepoints)
% 


%% Setup

% extract time windows: -300:Instruct+520ms, -200:Move+500ms
X = prep_pca(unit_data, configs, 'both', 'rest');
X_l_hand = cat(3,X.l_hand{:});
X_l_hand = X_l_hand([11:52, 53:88],:,:);
X_r_hand = cat(3,X.r_hand{:});
X_r_hand = X_r_hand([11:52, 53:88],:,:);


%% Iterate over all combinations of timepoints, calculating and storing the 
% principal covariance alignment between the two

for a = 1:78
    l_hand_a = squeeze(X_l_hand(a,:,:))';
    r_hand_a = squeeze(X_r_hand(a,:,:))';
    
    for b = 1:78
        l_hand_b = squeeze(X_l_hand(b,:,:))';
        r_hand_b = squeeze(X_r_hand(b,:,:))';
        
        % Train PCA models and calculate cross-validated alignment
        [alignment, ~, ~] = ...
            crossval_pca_align(l_hand_a, l_hand_b, 1, p,100);
        princ_align_self.l_hand(a,b) = mean(alignment);
        
        [alignment, ~, ~] = ...
            crossval_pca_align(r_hand_a, r_hand_b, 1, p,100);
        princ_align_self.r_hand(a,b) = mean(alignment);
        
        [alignment, ~, ~] = ...
            crossval_pca_align(l_hand_a, r_hand_b, 0, p,100);
        princ_align_lr.l_on_r(a,b) = mean(alignment);
        
        [alignment, ~, ~] = ...
            crossval_pca_align(r_hand_a, l_hand_b, 0, p,100);
        princ_align_lr.r_on_l(a,b) = mean(alignment);
        
    end
    
    %% display progress
    percent_complete = floor(a/0.78);
    if ~mod(percent_complete,10)
        time = clock;
        disp(num2str(time(4:5)))
        disp(strcat(num2str(percent_complete), '% complete'))
    end
    
end

    
    