function [X] = prep_pca(unit_data, configs, hem, norm_method)

%
% Collects firing rate data into time-locked and trial-separated data
% matrices for performing cross-validated PCA.
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
% hem - hemisphere to use. 'left', 'right', or 'both'
%
% norm_method - selected normalization method for unit firing rates. 
%               'rest' for resting std + 1hz
%               'range' for full firing rate range + 5hz
%
% OUTPUTS:
%
% X - for each hand, time-locked and trial-separated firing rate traces for
%     all simultaneously recorded units
%     (struct: 1x1, fields 'l_hand' and 'r_hand'
%      nested cell arrays: num_trials x 1
%      nested matrices: num_samples x num_units)
% 


%% Set up data matrices

% preallocate matrices for fine-timescale PCA
% isolate only non-repeated units
unit_data = unit_data([unit_data.repeat]==false);
num_units = length(unit_data);

num_trials_l_hand = 0;
num_trials_r_hand = 0;
for config = configs
    num_trials_l_hand = num_trials_l_hand + size(...
        vertcat(unit_data(1).ipsi.config(config).target.rest), 1);
    num_trials_r_hand = num_trials_r_hand + size(...
        vertcat(unit_data(1).contra.config(config).target.rest), 1);
end

X.l_hand = zeros([113, num_units, num_trials_l_hand]);
X.r_hand = zeros([113, num_units, num_trials_r_hand]);

l_hem_idx = [unit_data.hem]==0;
r_hem_idx = [unit_data.hem]==1;


%% Log the data

% loop over all units
for unit = num_units:-1:1
    
    %% log all the time-locked firing rates 
    % (not the most efficient way of doing so)
    l_hand = [];
    r_hand = [];
    rest_ipsi = [];
    rest_contra = [];
    
    % if left hemisphere
    if unit_data(unit).hem == 0
        for config = configs
            l_hand = [l_hand; ...
                [vertcat(unit_data(unit).ipsi.config(config).target.rest),...
                vertcat(unit_data(unit).ipsi.config(config).target.prep), ...
                vertcat(unit_data(unit).ipsi.config(config).target.move)]];
            r_hand = [r_hand; ...
                [vertcat(unit_data(unit).contra.config(config).target.rest),...
                vertcat(unit_data(unit).contra.config(config).target.prep), ...
                vertcat(unit_data(unit).contra.config(config).target.move)]];
        end
    
    % if right hemisphere
    elseif unit_data(unit).hem == 1
        for config = configs
            r_hand = [r_hand; ...
                [vertcat(unit_data(unit).ipsi.config(config).target.rest),...
                vertcat(unit_data(unit).ipsi.config(config).target.prep), ...
                vertcat(unit_data(unit).ipsi.config(config).target.move)]];
            l_hand = [l_hand; ...
                [vertcat(unit_data(unit).contra.config(config).target.rest),...
                vertcat(unit_data(unit).contra.config(config).target.prep), ...
                vertcat(unit_data(unit).contra.config(config).target.move)]];
        end
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
            
            X.l_hand(:,unit,:) = ...
                (l_hand' - mu_rest(unit))/sigma_rest(unit);
            X.r_hand(:,unit,:) = ...
                (r_hand' - mu_rest(unit))/sigma_rest(unit);
            
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
            X.l_hand(:,unit,:) = ...
                (l_hand' - mu_rest(unit))/fr_range;
            X.r_hand(:,unit,:) = ...
                (r_hand' - mu_rest(unit))/fr_range;
                     
        otherwise
            error('Error. Invalid normalization method name.')
            
    end
    
    
end


%% reorganize into output data structure

switch hem
    case 'left'
        X.l_hand = X.l_hand(:,l_hem_idx,:);
        X.r_hand = X.r_hand(:,l_hem_idx,:);
    case 'right'
        X.l_hand = X.l_hand(:,r_hem_idx,:);
        X.r_hand = X.r_hand(:,r_hem_idx,:);
    case 'both'
        X.l_hand = X.l_hand;
        X.r_hand = X.r_hand;
    otherwise
        error('Error. Invalid hand parameter name.')
end

X.l_hand = squeeze(num2cell(X.l_hand, [1,2]));
X.r_hand = squeeze(num2cell(X.r_hand, [1,2]));






