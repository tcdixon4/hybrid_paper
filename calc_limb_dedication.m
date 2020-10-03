function [arm_pref, modulation, mu_rest, sigma_rest, p_mod]...
    = calc_limb_dedication(unit_data, config)

%
% Calculates metrics for each unit, including arm preference, modulation,
% etc
%
%
% INPUTS: 
%
% unit_data - unit-separated data struct containing hemisphere, brain area,
%             firing rate and other metrics
%             (struct: 1 x num_units)
%
% OUTPUTS:
%
% arm_pref - arm preference, scalar for each unit x phase
%            (struct: 1x1, fields 'rest', 'prep', 'move'
%             nested vector: num_units x 1)
% 
% modulation - modulation strength, scalar for each unit x phase x hand
%              (struct: 1x1, fields 'rest', 'prep', 'move'
%               nested matrix: num_units x 2, col 1 left hand, col 2 right)
% 
% mu_rest - mean firing rate at rest, scalar for each unit
%           (vector: num_units x 1)
% 
% sigma_rest - standard deviation of firing rate at rest, softened by 
%              adding 1, scalar for each unit
%              (vector: num_units x 1)
% 
% p_mod - p-values for modulation differences
%         (struct: 1x1, fields 'ipsi', 'contra', 'hand_discrim_rest',
%          'd_hand'
%          nested in ipsi/contra, struct: 1x1, fields 'prep', 'move'
%              nested vector: num_units x 1, p-vals for mod above rest)
%          hand_discrim_rest vector: num_units x 1, p-vals for mod
%              differences between hands at rest
%          d_hand vector: 1 x num_units, difference in mean firing rate 
%              between hands at rest, abs value
%


%% Set up data matrices and helper variables

num_units = length(unit_data);
% create indexing vectors for each hemisphere
left_hem_idx = [unit_data.hem]==0;
right_hem_idx = [unit_data.hem]==1;
% preallocate
arm_pref.rest = zeros(num_units,1);
modulation.rest = zeros(num_units,1);
arm_pref.prep = zeros(num_units,1);
modulation.prep = zeros(num_units,1);
p_mod.ipsi.prep = zeros(num_units,1);
p_mod.contra.prep = zeros(num_units,1);
arm_pref.move = zeros(num_units,1);
modulation.move = zeros(num_units,1);
p_mod.ipsi.move = zeros(num_units,1);
p_mod.contra.move = zeros(num_units,1);
p_mod.hand_discrim_rest = zeros(num_units,1);

ipsi = zeros(num_units,113);
contra = zeros(num_units,113);


%% Log average modulation across time for all units

for unit = num_units:-1:1
    % use the config with the stationary hand in the eccentric position or
    % the center position, depending on the input
    if config == 0
        ipsi_config = 2;
        contra_config = 2;
    else
        ipsi_config = 3-2*unit_data(unit).hem;
        contra_config = 1+2*unit_data(unit).hem;
    end
    
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
    ipsi_rest_fr = fr_ipsi(1:trials_to_use, 11:26);
    contra_rest_fr = fr_contra(1:trials_to_use, 11:26);
    
    mu_rest(unit) = mean([ipsi_rest_fr(:); contra_rest_fr(:)]);
    sigma_rest(unit) = ...
        sqrt(mean([var(ipsi_rest_fr(:)), var(contra_rest_fr(:))])) + 1;
    
    % Z-score
    z_ipsi = (fr_ipsi-mu_rest(unit))/sigma_rest(unit);
    z_contra = (fr_contra-mu_rest(unit))/sigma_rest(unit);
    
    % for each trial, compute average modulation in each epoch
    ipsi_rest = z_ipsi(1:trials_to_use, 11:26).^2;
    contra_rest = z_contra(1:trials_to_use, 11:26).^2;
    ipsi_prep = z_ipsi(1:trials_to_use, 37:52).^2;
    contra_prep = z_contra(1:trials_to_use, 37:52).^2;
    ipsi_move = z_ipsi(1:trials_to_use, 63:78).^2;
    contra_move = z_contra(1:trials_to_use, 63:78).^2;
    
    [~, p_mod.ipsi.prep(unit)] = ttest2(mean(ipsi_prep,2),...
        [mean(ipsi_rest,2);mean(contra_rest,2)]);
    [~, p_mod.contra.prep(unit)] = ttest2(mean(contra_prep,2),...
        [mean(ipsi_rest,2);mean(contra_rest,2)]);
    [~, p_mod.ipsi.move(unit)] = ttest2(mean(ipsi_move,2),...
        [mean(ipsi_rest,2);mean(contra_rest,2)]);
    [~, p_mod.contra.move(unit)] = ttest2(mean(contra_move,2),...
        [mean(ipsi_rest,2);mean(contra_rest,2)]);
    [~, p_mod.hand_discrim_rest(unit)] = ...
        ttest2(mean(fr_ipsi(:, 11:26),2), mean(fr_contra(:, 11:26),2));
    p_mod.d_hand(unit) = abs(...
        mean(mean(fr_ipsi(:, 11:26),2)) - ...
        mean(mean(fr_contra(:, 11:26),2)));
    
    % log average modulation across time (mean squared Z-score firing rate)
    contra(unit,:) = mean(z_contra.^2);
    ipsi(unit,:) = mean(z_ipsi.^2);
    
end

% Reorganize into left and right hand response matrices
lhand = [ipsi(left_hem_idx,:); contra(right_hem_idx,:)];
rhand = [contra(left_hem_idx,:); ipsi(right_hem_idx,:)];


%% Calculate epoch-specific modulation and arm preference for each unit

modulation.rest = [mean(lhand(:,11:26),2), mean(rhand(:,11:26),2)];
modulation.prep = [mean(lhand(:,37:52),2), mean(rhand(:,37:52),2)];
modulation.move = [mean(lhand(:,63:78),2), mean(rhand(:,63:78),2)];

% positive arm_pref means right hand preferring
arm_pref.rest = (modulation.rest(:,2)-modulation.rest(:,1))./ ...
    (modulation.rest(:,1)+modulation.rest(:,2));
arm_pref.prep = (modulation.prep(:,2)-modulation.prep(:,1))./ ...
    (modulation.prep(:,1)+modulation.prep(:,2));
arm_pref.move = (modulation.move(:,2)-modulation.move(:,1))./ ...
    (modulation.move(:,1)+modulation.move(:,2));



