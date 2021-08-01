function [X, var_fr] = prep_pca_trial_avg(unit_data, norm_method)

%
% Collects firing rate data into time-locked and trial-averaged data
% matrices for performing cross-validated dPCA.
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
% OUTPUTS:
%
% X - for each hand, time-locked and trial-averaged firing rate traces for
%     all simultaneously recorded units
%     (struct: 1x1, fields 'l_hand' and 'r_hand'
%      nested struct: 1x3, field 'config'
%      nested struct: 1x1, fields 'rest', 'prep', 'move')
%      nested matrices: units x time_samples x targets
% 


%% loop over all units
for unit = length(unit_data):-1:1
    
    %% log all the time-locked, trial-averaged firing rates 
    
    % switch from ipsi/contra designation to left/right so that both
    % hemispheres can be combined with the appropriate trial sets
    if unit_data(unit).hem == 0
        l_hand = unit_data(unit).ipsi;
        r_hand = unit_data(unit).contra;
    elseif unit_data(unit).hem == 1
        l_hand = unit_data(unit).contra;
        r_hand = unit_data(unit).ipsi;
    end
    
    for config = 3:-1:1
        for target = 6:-1:1
            % collect the data into 4D trial-averaged matrices, indexed as:
            % (unit, time, target, config)
            % The time indices are as follows:
            % 1:16  - Rest, (Inst-300ms):(Inst), 16 samples at 50Hz
            % 17:42 - Prep, (Inst):(Inst+500ms), 26 samples at 50Hz
            % 43:78 - Move, (Move-200ms):(Move+ 500ms), 36 samples at 50Hz
            if isempty(l_hand.config(config).target(target).rest)
                l_hand_avg(unit,:,target,config) = nan(1,78);
            else
                l_hand_avg(unit,:,target,config) = [...
                    mean(l_hand.config(config).target(target).rest(:,end-15:end),1),...
                    mean(l_hand.config(config).target(target).prep(:,1:26),1),...
                    mean(l_hand.config(config).target(target).move(:,1:36),1)]';
            end
            if isempty(r_hand.config(config).target(target).rest)
                r_hand_avg(unit,:,target,config) = nan(1,78);
            else
                r_hand_avg(unit,:,target,config) = [...
                    mean(r_hand.config(config).target(target).rest(:,end-15:end),1),...
                    mean(r_hand.config(config).target(target).prep(:,1:26),1),...
                    mean(r_hand.config(config).target(target).move(:,1:36),1)]';                
            end
        end
        
    end
    
    lh = l_hand_avg(unit,:,:,3);
    rh = r_hand_avg(unit,:,:,1);
    var_fr(unit, 1) = var(lh(:));
    var_fr(unit, 2) = var(rh(:));
    
    %% compute normalizing constant, depending on the input either:
    % 'rest': the standard devation of NON-averaged Rest firing rates
    % 'range': the full range of trial-averaged firing rates
    if strcmp(norm_method, 'rest')
        rest_left = [vertcat(l_hand.config(1).target.rest);...
            vertcat(l_hand.config(2).target.rest);...
            vertcat(l_hand.config(3).target.rest)];
        rest_right = [vertcat(r_hand.config(1).target.rest);...
            vertcat(r_hand.config(2).target.rest);...
            vertcat(r_hand.config(3).target.rest)];
        rest_left = rest_left(:, 11:26);
        rest_right = rest_right(:, 11:26);
        mu = mean([mean2(rest_left), mean2(rest_right)]);
        norm_constant = mean([std2(rest_left), std2(rest_right)]) + 1;
    elseif strcmp(norm_method, 'range')
        avg_cat = [l_hand_avg(~isnan(l_hand_avg));...
            r_hand_avg(~isnan(r_hand_avg))];
        mu = mean(avg_cat);
        norm_constant = max(avg_cat,[],'all') - min(avg_cat,[],'all') + 5;
    end
    
    %% de-mean, normalize, and log in output structure
    for config = 3:-1:1
        X.l_hand.config(config).rest(unit,:,:) =...
            (l_hand_avg(unit,1:16,:,config) - mu)/norm_constant;
        X.l_hand.config(config).prep(unit,:,:) =...
            (l_hand_avg(unit,17:42,:,config) - mu)/norm_constant;
        X.l_hand.config(config).move(unit,:,:) =...
            (l_hand_avg(unit,43:78,:,config) - mu)/norm_constant;
        X.r_hand.config(config).rest(unit,:,:) =...
            (r_hand_avg(unit,1:16,:,config) - mu)/norm_constant;
        X.r_hand.config(config).prep(unit,:,:) =...
            (r_hand_avg(unit,17:42,:,config) - mu)/norm_constant;
        X.r_hand.config(config).move(unit,:,:) =...
            (r_hand_avg(unit,43:78,:,config) - mu)/norm_constant;
    end
     
end


