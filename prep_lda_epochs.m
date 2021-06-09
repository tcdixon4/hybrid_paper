function [X, Y] = prep_lda_epochs(unit_data)

%
% Collects firing rate data phase-specific and trial-separated data
% matrices for performing cross-validated LDA, separating the population
% into left-preferring and right-preferring prior to analysis. Also
% collects target id's for each trial to train and test classifiers with.
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
% X - neural data (predictors)
%     (struct: 1x1, fields 'l_pref' and 'r_pref'
%      nested struct: 1x1, fields 'l_hand' and 'r_hand'
%      nested struct: 1x1, fields 'rest', 'prep', and 'move'
%      nested matrices: num_trials x num_units)
%
% Y - target id's (classes)
%     (struct: 1x1, fields 'l_hand' and 'r_hand'
%      nested matrices: num_trials)
% 


%% Compute arm preferences

% isolate only non-repeated units
unit_data = unit_data([unit_data.repeat]==false);
[arm_pref, ~, ~, ~]...
    = calc_limb_dedication(unit_data, 0);

r_pref_idx.rest = arm_pref.rest>0;
r_pref_idx.prep = arm_pref.prep>0;
r_pref_idx.move = arm_pref.move>0;


%% Set up data matrices

% preallocate class and neural data matrices
num_units = length(unit_data);

num_trials_l_hand = ...
    size(vertcat(unit_data(1).ipsi.config(3).target.rest), 1);
num_trials_r_hand = ...
    size(vertcat(unit_data(1).contra.config(1).target.rest), 1);

l_hand.rest = zeros(num_trials_l_hand, num_units);
l_hand.prep = zeros(num_trials_l_hand, num_units);
l_hand.move = zeros(num_trials_l_hand, num_units);
r_hand.rest = zeros(num_trials_r_hand, num_units);
r_hand.prep = zeros(num_trials_r_hand, num_units);
r_hand.move = zeros(num_trials_r_hand, num_units);


%% log class (target) vectors

Y.l_hand = [];
Y.r_hand = [];

for targ = 1:6
    Y.l_hand = [Y.l_hand; ...
        targ*ones(size(unit_data(1).ipsi.config(3).target(targ).rest, 1), 1)];
    Y.r_hand = [Y.r_hand; ...
        targ*ones(size(unit_data(1).contra.config(1).target(targ).rest, 1), 1)];
end


%% Log the data

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
    
    % calculate rest period mu and std for Z-scoring, doesn't really matter
    % for LDA
    for c = 1:3 % use all configs for this part
        rest_ipsi = [rest_ipsi; ...
            vertcat(unit_data(unit).ipsi.config(c).target.rest)];
        rest_contra = [rest_contra;
            vertcat(unit_data(unit).contra.config(c).target.rest)];
    end
    rest_ipsi = rest_ipsi(:, 11:26);
    rest_contra = rest_contra(:, 11:26);
    mu_rest(unit) = mean([rest_ipsi(:); rest_contra(:)]);
    sigma_rest(unit) = ...
        sqrt(mean([var(rest_ipsi(:)), var(rest_contra(:))])) + 1;
    
    l_hand.rest(:,unit) = ...
        (mean(l_hand_holder(:,11:26),2) - mu_rest(unit))/sigma_rest(unit);
    l_hand.prep(:,unit) = ...
        (mean(l_hand_holder(:,37:52),2) - mu_rest(unit))/sigma_rest(unit);
    l_hand.move(:,unit) = ...
        (mean(l_hand_holder(:,63:78),2) - mu_rest(unit))/sigma_rest(unit);
    r_hand.rest(:,unit) = ...
        (mean(r_hand_holder(:,11:26),2) - mu_rest(unit))/sigma_rest(unit);
    r_hand.prep(:,unit) = ...
        (mean(r_hand_holder(:,37:52),2) - mu_rest(unit))/sigma_rest(unit);
    r_hand.move(:,unit) = ...
        (mean(r_hand_holder(:,63:78),2) - mu_rest(unit))/sigma_rest(unit);
    
end

X.l_pref.l_hand.rest = l_hand.rest(:,~r_pref_idx.rest);
X.l_pref.l_hand.prep = l_hand.prep(:,~r_pref_idx.prep);
X.l_pref.l_hand.move = l_hand.move(:,~r_pref_idx.move);
X.l_pref.r_hand.rest = r_hand.rest(:,~r_pref_idx.rest);
X.l_pref.r_hand.prep = r_hand.prep(:,~r_pref_idx.prep);
X.l_pref.r_hand.move = r_hand.move(:,~r_pref_idx.move);

X.r_pref.l_hand.rest = l_hand.rest(:,r_pref_idx.rest);
X.r_pref.l_hand.prep = l_hand.prep(:,r_pref_idx.prep);
X.r_pref.l_hand.move = l_hand.move(:,r_pref_idx.move);
X.r_pref.r_hand.rest = r_hand.rest(:,r_pref_idx.rest);
X.r_pref.r_hand.prep = r_hand.prep(:,r_pref_idx.prep);
X.r_pref.r_hand.move = r_hand.move(:,r_pref_idx.move);


