function [modulation]...
    = pref_nonpref_mod_across_time(unit_data)

%
% Collects modulation across time for each unit and separates into
% responses for the preferred arm or the non-preferred arm.
%
%
% INPUTS: 
%
% unit_data - unit-separated data struct containing hemisphere, brain area,
%             firing rate and other metrics
%             (struct: 1 x num_units)
%
%
% OUTPUTS:
%
% modulation - Modulation across time for all units
%              (struct: 1x1, fields 'pref', 'nonpref'
%               nested matrix: num_units*2 x num_samples)
%


%% Set up data matrices and helper variables

% isolate only non-repeated units
unit_data = unit_data([unit_data.repeat]==false);
num_units = length(unit_data);
% find division between hemispheres
left_hem_idx = [unit_data.hem]==0;
right_hem_idx = [unit_data.hem]==1;
% reorder unit_data in case left and right hem are intermixed
unit_data = [unit_data(left_hem_idx), unit_data(right_hem_idx)];
left_hem_idx = [unit_data.hem]==0;
right_hem_idx = [unit_data.hem]==1;
% compute arm preferences with the stationary hand center config
[arm_pref, ~, ~, ~] = calc_limb_dedication(unit_data, 0);
r_pref_idx.rest = arm_pref.rest>0;
r_pref_idx.prep = arm_pref.prep>0;
r_pref_idx.move = arm_pref.move>0;
% preallocate
ipsi = zeros(num_units,113);
contra = zeros(num_units,113);


%% Log average modulation across time for all units

for unit = num_units:-1:1
    % use the stationary hand eccentric config
    ipsi_config = 3-2*unit_data(unit).hem;
    contra_config = 1+2*unit_data(unit).hem;
    
    % log firing rates for each trial, separated by limb
    fr_ipsi = [];
    fr_contra = [];
    for target = 1:6
        if ~isempty(unit_data(unit).ipsi.config(ipsi_config)...
                .target(target).rest)
            fr_ipsi = [fr_ipsi;...
                [unit_data(unit).ipsi.config(ipsi_config)...
                .target(target).rest,...
                unit_data(unit).ipsi.config(ipsi_config)...
                .target(target).prep,...
                unit_data(unit).ipsi.config(ipsi_config)...
                .target(target).move]];
        end
        if ~isempty(unit_data(unit).contra.config(contra_config)...
                .target(target).rest)
            fr_contra = [fr_contra;...
                [unit_data(unit).contra.config(contra_config)...
                .target(target).rest,...
                unit_data(unit).contra.config(contra_config)...
                .target(target).prep,...
                unit_data(unit).contra.config(contra_config)...
                .target(target).move]];
        end
    end
    
    % calculate rest period mu and std for Z-scoring
    trials_to_use = min([size(fr_ipsi,1), size(fr_contra,1)]);
    ipsi_rest = fr_ipsi(1:trials_to_use, 11:26);
    contra_rest = fr_contra(1:trials_to_use, 11:26);
    mu_rest(unit) = mean([ipsi_rest(:); contra_rest(:)]);
    sigma_rest(unit) = ...
        sqrt(mean([var(ipsi_rest(:)), var(contra_rest(:))])) + 1;
    
    % Z-score
    z_ipsi = (fr_ipsi-mu_rest(unit))/sigma_rest(unit);
    z_contra = (fr_contra-mu_rest(unit))/sigma_rest(unit);
    
    % log average modulation across time (mean squared Z-score firing rate)
    contra(unit,:) = mean(z_contra.^2);
    ipsi(unit,:) = mean(z_ipsi.^2);
    
end

% Reorganize into left and right hand response matrices
lhand = [ipsi(left_hem_idx,:); contra(right_hem_idx,:)];
rhand = [contra(left_hem_idx,:); ipsi(right_hem_idx,:)];


%% Calculate epoch-specific modulation and separate based on arm preference

% replace non-preferring units with nan's to isolate preferring units
lhand_holder = lhand;
lhand_holder(r_pref_idx.rest, 1:26) = nan;
lhand_holder(r_pref_idx.prep, 27:52) = nan;
lhand_holder(r_pref_idx.move, 53:end) = nan;
rhand_holder = rhand;
rhand_holder(~r_pref_idx.rest, 1:26) = nan;
rhand_holder(~r_pref_idx.prep, 27:52) = nan;
rhand_holder(~r_pref_idx.move, 53:end) = nan;

modulation.pref = [lhand_holder; rhand_holder];

% replace preferring units with nan's to isolate non-preferring units
lhand_holder = lhand;
lhand_holder(~r_pref_idx.rest, 1:26) = nan;
lhand_holder(~r_pref_idx.prep, 27:52) = nan;
lhand_holder(~r_pref_idx.move, 53:end) = nan;
rhand_holder = rhand;
rhand_holder(r_pref_idx.rest, 1:26) = nan;
rhand_holder(r_pref_idx.prep, 27:52) = nan;
rhand_holder(r_pref_idx.move, 53:end) = nan;

modulation.nonpref = [lhand_holder; rhand_holder];




