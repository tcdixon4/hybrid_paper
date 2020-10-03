function [coeff, explained, hem] = ...
    crossval_pca_var_and_coeffs(unit_data, phase, norm_method)

%
% Computes variance preserved through PCA projections and returns the
% coefficients of those PCA models, representing the model structure.
%
%
% INPUTS: 
%
% unit_data - unit-separated data struct containing hemisphere, brain area,
%             firing rate and other metrics
%             (struct: 1 x num_units)
%
% phase - selected phase of the task to analyze. 'rest', 'prep' or 'move'
%
% norm_method - selected normalization method for unit firing rates. 
%               'rest' for resting std + 1hz
%               'range' for full firing rate range + 5hz
%
%
% OUTPUTS:
%
% coeff - PCA coefficient matrix
%         (struct: 1x1, fields 'l_hand', 'r_hand'
%          nested matrix: num_units x num_units)
%
% explained - proportion of variance accounted for with each train/test
%             condition combination
%             (struct: 1x1, fields 'l_hand', 'r_hand'
%              nested struct: 1x1, fields 'l_space', 'r_space')
% 
% hem - hemisphere that each unit is in
%       (vector: num_units x 1)
%


%% Prepare the neural data

X.cen = prep_pca(unit_data, 2, 'both', norm_method);
X_holder_lh = prep_pca(unit_data, 3, 'both', norm_method);
X_holder_rh = prep_pca(unit_data, 1, 'both', norm_method);
X.ecc.l_hand = X_holder_lh.l_hand;
X.ecc.r_hand = X_holder_rh.r_hand;
hem = [unit_data.hem];

% clearvars unit_data X_holder_lh X_holder_rh


%% Extract just the data from the phase of interest

switch phase
    case 'rest'
        idx = 11:26;
    case 'prep'
        idx = 37:52;
    case 'move'
        idx = 63:78;
    case 'all'
        idx = [];
end

if ~isempty(idx)
    for trial = 1:length(X.cen.l_hand)
        X.cen.l_hand{trial} = X.cen.l_hand{trial}(idx,:);
    end
    for trial = 1:length(X.cen.r_hand)
        X.cen.r_hand{trial} = X.cen.r_hand{trial}(idx,:);
    end
    for trial = 1:length(X.ecc.l_hand)
        X.ecc.l_hand{trial} = X.ecc.l_hand{trial}(idx,:);
    end
    for trial = 1:length(X.ecc.r_hand)
        X.ecc.r_hand{trial} = X.ecc.r_hand{trial}(idx,:);
    end
end


%% Fit the models on the eccentric configuration

l_hand_cen_mat = vertcat(X.cen.l_hand{:});
r_hand_cen_mat = vertcat(X.cen.r_hand{:});

[l_hand_coeff, ~,~,~,...
    ~, ~] = pca(l_hand_cen_mat, 'Centered',false);
[r_hand_coeff, ~,~,~,...
    ~, ~] = pca(r_hand_cen_mat, 'Centered',false);


%% test the models on the center configuration

l_hand_ecc_mat = vertcat(X.ecc.l_hand{:});
r_hand_ecc_mat = vertcat(X.ecc.r_hand{:});

% Project onto both the native and cross-limb spaces
l_hand_l_space_score = l_hand_ecc_mat*l_hand_coeff;
r_hand_r_space_score = r_hand_ecc_mat*r_hand_coeff;
l_hand_r_space_score = l_hand_ecc_mat*r_hand_coeff;
r_hand_l_space_score = r_hand_ecc_mat*l_hand_coeff;

% calculate the variance explained by each PC
explained.l_hand.l_space = var(l_hand_l_space_score)/...
    sum(var(l_hand_l_space_score));
explained.r_hand.r_space = var(r_hand_r_space_score)/...
    sum(var(r_hand_r_space_score));
explained.l_hand.r_space = var(l_hand_r_space_score)/...
    sum(var(l_hand_r_space_score));
explained.r_hand.l_space = var(r_hand_l_space_score)/...
    sum(var(r_hand_l_space_score));

coeff.l_hand = l_hand_coeff;
coeff.r_hand = r_hand_coeff;





