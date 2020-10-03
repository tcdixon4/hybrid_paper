function [X_ecc, X_cen] = prep_pca_pref_epochs(unit_data, norm_method)

%
% Collects firing rate data from two independent datasets into 
% phase-specific matrices for performing PCA during trials of the 
% "preferred" and "non-preferred" arms of sub-populations divided based on
% arm-preference.
%
%
% INPUTS: 
%
% unit_data - unit-separated data struct containing hemisphere, brain area,
%             firing rate and other metrics
%             (struct: 1 x num_units)
%
% norm_method - selected normalization method for unit firing rates. 
%               'rest' for resting std + 1hz
%               'range' for full firing rate range + 5hz
%
%
% OUTPUTS:
%
% X_ecc - neural data from eccentric starting config
%         (struct: 1x1, fields 'l_pref' and 'r_pref'
%          nested struct: 1x1, fields 'l_hand' and 'r_hand'
%          nested struct: 1x1, fields 'rest', 'prep', and 'move'
%          nested matrices: num_trials*num_samples x num_units)
%
% X_cen - neural data from centered starting config
%         (struct: 1x1, fields 'l_pref' and 'r_pref'
%          nested struct: 1x1, fields 'l_hand' and 'r_hand'
%          nested struct: 1x1, fields 'rest', 'prep', and 'move'
%          nested matrices: num_trials*num_samples x num_units)
% 


%% Compute arm preferences

[arm_pref, ~, ~, ~]...
    = calc_limb_dedication(unit_data, 0);

r_pref_idx.rest = arm_pref.rest>0;
r_pref_idx.prep = arm_pref.prep>0;
r_pref_idx.move = arm_pref.move>0;


%% Set up data matrices for center config

% preallocate matrices for fine-timescale PCA
num_units = length(unit_data);

num_trials_l_hand = ...
    size(vertcat(unit_data(1).ipsi.config(2).target.rest), 1);
num_trials_r_hand = ...
    size(vertcat(unit_data(1).contra.config(2).target.rest), 1);

l_hand.rest = zeros(num_trials_l_hand*16, num_units);
l_hand.prep = zeros(num_trials_l_hand*16, num_units);
l_hand.move = zeros(num_trials_l_hand*16, num_units);
r_hand.rest = zeros(num_trials_r_hand*16, num_units);
r_hand.prep = zeros(num_trials_r_hand*16, num_units);
r_hand.move = zeros(num_trials_r_hand*16, num_units);


%% Log the data for the center config

% loop over all units
for unit = num_units:-1:1
    l_hand_holder = [];
    r_hand_holder = [];
    rest_ipsi = [];
    rest_contra = [];
    
    % if left hemisphere
    if unit_data(unit).hem == 0
        l_hand_holder = [l_hand_holder; ...
            [vertcat(unit_data(unit).ipsi.config(2).target.rest),...
            vertcat(unit_data(unit).ipsi.config(2).target.prep), ...
            vertcat(unit_data(unit).ipsi.config(2).target.move)]];
        r_hand_holder = [r_hand_holder; ...
            [vertcat(unit_data(unit).contra.config(2).target.rest),...
            vertcat(unit_data(unit).contra.config(2).target.prep), ...
            vertcat(unit_data(unit).contra.config(2).target.move)]];
        
    % if right hemisphere
    elseif unit_data(unit).hem == 1
        r_hand_holder = [r_hand_holder; ...
            [vertcat(unit_data(unit).ipsi.config(2).target.rest),...
            vertcat(unit_data(unit).ipsi.config(2).target.prep), ...
            vertcat(unit_data(unit).ipsi.config(2).target.move)]];
        l_hand_holder = [l_hand_holder; ...
            [vertcat(unit_data(unit).contra.config(2).target.rest),...
            vertcat(unit_data(unit).contra.config(2).target.prep), ...
            vertcat(unit_data(unit).contra.config(2).target.move)]];
    end
    
    %% normalize unit firing rates
    
    % calculate rest period mean to subtract out
    for c = 1:3 % use all configs for this part
        rest_ipsi = [rest_ipsi; ...
            vertcat(unit_data(unit).ipsi.config(c).target.rest)];
        rest_contra = [rest_contra;
            vertcat(unit_data(unit).contra.config(c).target.rest)];
    end
    rest_ipsi = rest_ipsi(:, 11:26);
    rest_contra = rest_contra(:, 11:26);
    mu_rest(unit) = mean([rest_ipsi(:); rest_contra(:)]);
    
    switch norm_method
        
        case 'rest'
            % calculate rest period std for Z-scoring
            sigma_rest(unit) = ...
                sqrt(mean([var(rest_ipsi(:)), var(rest_contra(:))])) + 1;
            
            tmp = (l_hand_holder(:,11:26) - mu_rest(unit))/sigma_rest(unit);
            l_hand.rest(:,unit) = tmp(:);
            tmp = (l_hand_holder(:,37:52) - mu_rest(unit))/sigma_rest(unit);
            l_hand.prep(:,unit) = tmp(:);
            tmp = (l_hand_holder(:,63:78) - mu_rest(unit))/sigma_rest(unit);
            l_hand.move(:,unit) = tmp(:);
            tmp = (r_hand_holder(:,11:26) - mu_rest(unit))/sigma_rest(unit);
            r_hand.rest(:,unit) = tmp(:);
            tmp = (r_hand_holder(:,37:52) - mu_rest(unit))/sigma_rest(unit);
            r_hand.prep(:,unit) = tmp(:);
            tmp = (r_hand_holder(:,63:78) - mu_rest(unit))/sigma_rest(unit);
            r_hand.move(:,unit) = tmp(:);
            
        case 'range'
            % calculate overall firing rate range of trial-averaged data
            % for normalization
            fr_cat = nan(1,36*16*3);
            for config = 1:3
                for target = 1:6
                    idx = ((config-1)*6 + target - 1)*48 + 1;
                    % ipsi
                    if isempty(unit_data(unit).ipsi.config(config).target(target).rest)
                        fr_cat(idx:(idx+47)) = nan(1,48);
                    else
                        fr = [unit_data(unit).ipsi.config(config).target(target).rest,...
                            unit_data(unit).ipsi.config(config).target(target).prep,...
                            unit_data(unit).ipsi.config(config).target(target).move];
                        fr_cat(idx:(idx+47)) = mean(fr(:,[11:26,37:52,63:78]));
                    end
                    % contra
                    if isempty(unit_data(unit).contra.config(config).target(target).rest)
                        fr_cat((idx+864):(idx+911)) = nan(1,48);
                    else
                        fr = [unit_data(unit).contra.config(config).target(target).rest,...
                            unit_data(unit).contra.config(config).target(target).prep,...
                            unit_data(unit).contra.config(config).target(target).move];
                        fr_cat((idx+864):(idx+911)) = mean(fr(:,[11:26,37:52,63:78]));
                    end
                end
            end
            
            fr_range = range(fr_cat)+5;
            tmp = (l_hand_holder(:,11:26) - mu_rest(unit))/fr_range;
            l_hand.rest(:,unit) = tmp(:);
            tmp = (l_hand_holder(:,37:52) - mu_rest(unit))/fr_range;
            l_hand.prep(:,unit) = tmp(:);
            tmp = (l_hand_holder(:,63:78) - mu_rest(unit))/fr_range;
            l_hand.move(:,unit) = tmp(:);
            tmp = (r_hand_holder(:,11:26) - mu_rest(unit))/fr_range;
            r_hand.rest(:,unit) = tmp(:);
            tmp = (r_hand_holder(:,37:52) - mu_rest(unit))/fr_range;
            r_hand.prep(:,unit) = tmp(:);
            tmp = (r_hand_holder(:,63:78) - mu_rest(unit))/fr_range;
            r_hand.move(:,unit) = tmp(:);
            
        otherwise
            error('Error. Invalid normalization method name.')
            
    end
    
end

X_cen.l_pref.l_hand.rest = l_hand.rest(:,~r_pref_idx.rest);
X_cen.l_pref.l_hand.prep = l_hand.prep(:,~r_pref_idx.prep);
X_cen.l_pref.l_hand.move = l_hand.move(:,~r_pref_idx.move);
X_cen.l_pref.r_hand.rest = r_hand.rest(:,~r_pref_idx.rest);
X_cen.l_pref.r_hand.prep = r_hand.prep(:,~r_pref_idx.prep);
X_cen.l_pref.r_hand.move = r_hand.move(:,~r_pref_idx.move);

X_cen.r_pref.l_hand.rest = l_hand.rest(:,r_pref_idx.rest);
X_cen.r_pref.l_hand.prep = l_hand.prep(:,r_pref_idx.prep);
X_cen.r_pref.l_hand.move = l_hand.move(:,r_pref_idx.move);
X_cen.r_pref.r_hand.rest = r_hand.rest(:,r_pref_idx.rest);
X_cen.r_pref.r_hand.prep = r_hand.prep(:,r_pref_idx.prep);
X_cen.r_pref.r_hand.move = r_hand.move(:,r_pref_idx.move);


%% Set up data matrices for eccentric config

% preallocate matrices for fine-timescale PCA
num_units = length(unit_data);

num_trials_l_hand = ...
    size(vertcat(unit_data(1).ipsi.config(3).target.rest), 1);
num_trials_r_hand = ...
    size(vertcat(unit_data(1).contra.config(1).target.rest), 1);

l_hand.rest = zeros(num_trials_l_hand*16, num_units);
l_hand.prep = zeros(num_trials_l_hand*16, num_units);
l_hand.move = zeros(num_trials_l_hand*16, num_units);
r_hand.rest = zeros(num_trials_r_hand*16, num_units);
r_hand.prep = zeros(num_trials_r_hand*16, num_units);
r_hand.move = zeros(num_trials_r_hand*16, num_units);


%% Log the data for the eccentric config

% loop over all units
for unit = num_units:-1:1
    l_hand_holder = [];
    r_hand_holder = [];
    rest_ipsi = [];
    rest_contra = [];
    
    % if left hemisphere
    if unit_data(unit).hem == 0
        l_hand_holder = [l_hand_holder; ...
            [vertcat(unit_data(unit).ipsi.config(3).target.rest),...
            vertcat(unit_data(unit).ipsi.config(3).target.prep), ...
            vertcat(unit_data(unit).ipsi.config(3).target.move)]];
        r_hand_holder = [r_hand_holder; ...
            [vertcat(unit_data(unit).contra.config(1).target.rest),...
            vertcat(unit_data(unit).contra.config(1).target.prep), ...
            vertcat(unit_data(unit).contra.config(1).target.move)]];
        
    % if right hemisphere
    elseif unit_data(unit).hem == 1
        r_hand_holder = [r_hand_holder; ...
            [vertcat(unit_data(unit).ipsi.config(1).target.rest),...
            vertcat(unit_data(unit).ipsi.config(1).target.prep), ...
            vertcat(unit_data(unit).ipsi.config(1).target.move)]];
        l_hand_holder = [l_hand_holder; ...
            [vertcat(unit_data(unit).contra.config(3).target.rest),...
            vertcat(unit_data(unit).contra.config(3).target.prep), ...
            vertcat(unit_data(unit).contra.config(3).target.move)]];
    end
    
    
    %% normalize unit firing rates
    
    % calculate rest period mean to subtract out
    for c = 1:3 % use all configs for this part
        rest_ipsi = [rest_ipsi; ...
            vertcat(unit_data(unit).ipsi.config(c).target.rest)];
        rest_contra = [rest_contra;
            vertcat(unit_data(unit).contra.config(c).target.rest)];
    end
    rest_ipsi = rest_ipsi(:, 11:26);
    rest_contra = rest_contra(:, 11:26);
    mu_rest(unit) = mean([rest_ipsi(:); rest_contra(:)]);
    
    switch norm_method
        
        case 'rest'
            % calculate rest period std for Z-scoring
            sigma_rest(unit) = ...
                sqrt(mean([var(rest_ipsi(:)), var(rest_contra(:))])) + 1;
            
            tmp = (l_hand_holder(:,11:26) - mu_rest(unit))/sigma_rest(unit);
            l_hand.rest(:,unit) = tmp(:);
            tmp = (l_hand_holder(:,37:52) - mu_rest(unit))/sigma_rest(unit);
            l_hand.prep(:,unit) = tmp(:);
            tmp = (l_hand_holder(:,63:78) - mu_rest(unit))/sigma_rest(unit);
            l_hand.move(:,unit) = tmp(:);
            tmp = (r_hand_holder(:,11:26) - mu_rest(unit))/sigma_rest(unit);
            r_hand.rest(:,unit) = tmp(:);
            tmp = (r_hand_holder(:,37:52) - mu_rest(unit))/sigma_rest(unit);
            r_hand.prep(:,unit) = tmp(:);
            tmp = (r_hand_holder(:,63:78) - mu_rest(unit))/sigma_rest(unit);
            r_hand.move(:,unit) = tmp(:);
            
        case 'range'
            % calculate overall firing rate range of trial-averaged data
            % for normalization
            fr_cat = nan(1,36*16*3);
            for config = 1:3
                for target = 1:6
                    idx = ((config-1)*6 + target - 1)*48 + 1;
                    % ipsi
                    if isempty(unit_data(unit).ipsi.config(config).target(target).rest)
                        fr_cat(idx:(idx+47)) = nan(1,48);
                    else
                        fr = [unit_data(unit).ipsi.config(config).target(target).rest,...
                            unit_data(unit).ipsi.config(config).target(target).prep,...
                            unit_data(unit).ipsi.config(config).target(target).move];
                        fr_cat(idx:(idx+47)) = mean(fr(:,[11:26,37:52,63:78]));
                    end
                    % contra
                    if isempty(unit_data(unit).contra.config(config).target(target).rest)
                        fr_cat((idx+864):(idx+911)) = nan(1,48);
                    else
                        fr = [unit_data(unit).contra.config(config).target(target).rest,...
                            unit_data(unit).contra.config(config).target(target).prep,...
                            unit_data(unit).contra.config(config).target(target).move];
                        fr_cat((idx+864):(idx+911)) = mean(fr(:,[11:26,37:52,63:78]));
                    end
                end
            end
            
            fr_range = range(fr_cat)+5;
            tmp = (l_hand_holder(:,11:26) - mu_rest(unit))/fr_range;
            l_hand.rest(:,unit) = tmp(:);
            tmp = (l_hand_holder(:,37:52) - mu_rest(unit))/fr_range;
            l_hand.prep(:,unit) = tmp(:);
            tmp = (l_hand_holder(:,63:78) - mu_rest(unit))/fr_range;
            l_hand.move(:,unit) = tmp(:);
            tmp = (r_hand_holder(:,11:26) - mu_rest(unit))/fr_range;
            r_hand.rest(:,unit) = tmp(:);
            tmp = (r_hand_holder(:,37:52) - mu_rest(unit))/fr_range;
            r_hand.prep(:,unit) = tmp(:);
            tmp = (r_hand_holder(:,63:78) - mu_rest(unit))/fr_range;
            r_hand.move(:,unit) = tmp(:);
            
        otherwise
            error('Error. Invalid normalization method name.')
            
    end
    
end

X_ecc.l_pref.l_hand.rest = l_hand.rest(:,~r_pref_idx.rest);
X_ecc.l_pref.l_hand.prep = l_hand.prep(:,~r_pref_idx.prep);
X_ecc.l_pref.l_hand.move = l_hand.move(:,~r_pref_idx.move);
X_ecc.l_pref.r_hand.rest = r_hand.rest(:,~r_pref_idx.rest);
X_ecc.l_pref.r_hand.prep = r_hand.prep(:,~r_pref_idx.prep);
X_ecc.l_pref.r_hand.move = r_hand.move(:,~r_pref_idx.move);

X_ecc.r_pref.l_hand.rest = l_hand.rest(:,r_pref_idx.rest);
X_ecc.r_pref.l_hand.prep = l_hand.prep(:,r_pref_idx.prep);
X_ecc.r_pref.l_hand.move = l_hand.move(:,r_pref_idx.move);
X_ecc.r_pref.r_hand.rest = r_hand.rest(:,r_pref_idx.rest);
X_ecc.r_pref.r_hand.prep = r_hand.prep(:,r_pref_idx.prep);
X_ecc.r_pref.r_hand.move = r_hand.move(:,r_pref_idx.move);


