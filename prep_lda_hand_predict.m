function [X, Y, unit_area] = prep_lda_hand_predict(unit_data, trials)

% make sure trials already contains only those which you would like to
% analyze. operates on a single session of `unit_data` and `trials`


%% Setup

num_units = length(unit_data);
% find division between hemispheres
left_hem_idx = [unit_data.hem]==0;
right_hem_idx = [unit_data.hem]==1;
% reorder unit_data in case left and right hem are intermixed
unit_data = [unit_data(left_hem_idx), unit_data(right_hem_idx)];

num_trials = length(trials);


%% log class (target x config x hand) vectors
% trial classes, coded as:
% 1: target 1, config 1, left hand     19: target 1, config 1, right hand
% 2: target 2, config 1, left hand     20: target 2, config 1, right hand
% 3: target 3, config 1, left hand     21: target 3, config 1, right hand
% 4: target 4, config 1, left hand     22: target 4, config 1, right hand
% 5: target 5, config 1, left hand     23: target 5, config 1, right hand
% 6: target 6, config 1, left hand X   24: target 6, config 1, right hand
% 7: target 1, config 2, left hand     25: target 1, config 2, right hand
% 8: target 2, config 2, left hand     26: target 2, config 2, right hand
% 9: target 3, config 2, left hand     27: target 3, config 2, right hand
% 10: target 4, config 2, left hand    28: target 4, config 2, right hand X
% 11: target 5, config 2, left hand    29: target 5, config 2, right hand
% 12: target 6, config 2, left hand X  30: target 6, config 2, right hand
% 13: target 1, config 3, left hand    31: target 1, config 3, right hand
% 14: target 2, config 3, left hand    32: target 2, config 3, right hand
% 15: target 3, config 3, left hand    33: target 3, config 3, right hand
% 16: target 4, config 3, left hand    34: target 4, config 3, right hand X
% 17: target 5, config 3, left hand    35: target 5, config 3, right hand
% 18: target 6, config 3, left hand    36: target 6, config 3, right hand

Y = zeros(num_trials,2);
% trial class assignments
Y(:,1) = [trials.target_id] + ...
    ([trials.config]+1)*6 + ...
    ([trials.reach_hand]-1)*18;
% switching trials
Y(:,2) = [false,(diff([trials.config])~=0 | diff([trials.reach_hand])~=0)];

% reorganize to match the trial order for the neural data (predictors) that
% will be logged in the next code block - chunking targets within configs
% within hand
reorder_idx = [];
for class_lbl = 1:36
    reorder_idx = [reorder_idx; find(Y(:,1)==class_lbl)];
end
Y = Y(reorder_idx,:);

clearvars trials


%% Log the neural data (predictors)

% loop over all units
for unit = num_units:-1:1
   
    % first, log which area each unit was in for model structure analysis
    % (PMd=0, M1=1)
    unit_area(unit,1) = unit_data(unit).area;
    
    l_hand_holder = [];
    r_hand_holder = [];
    
    for config = 1:3
        
        % if left hemisphere
        if unit_data(unit).hem == 0
            l_hand_holder = [l_hand_holder;...
                vertcat(unit_data(unit).ipsi.config(config).target.rest),...
                vertcat(unit_data(unit).ipsi.config(config).target.prep), ...
                vertcat(unit_data(unit).ipsi.config(config).target.move)];
            r_hand_holder = [r_hand_holder;...
                vertcat(unit_data(unit).contra.config(config).target.rest),...
                vertcat(unit_data(unit).contra.config(config).target.prep), ...
                vertcat(unit_data(unit).contra.config(config).target.move)];
            
        % if right hemisphere
        elseif unit_data(unit).hem == 1
            r_hand_holder = [r_hand_holder;...
                vertcat(unit_data(unit).ipsi.config(config).target.rest),...
                vertcat(unit_data(unit).ipsi.config(config).target.prep), ...
                vertcat(unit_data(unit).ipsi.config(config).target.move)];
            l_hand_holder = [l_hand_holder;...
                vertcat(unit_data(unit).contra.config(config).target.rest),...
                vertcat(unit_data(unit).contra.config(config).target.prep), ...
                vertcat(unit_data(unit).contra.config(config).target.move)];
        end
        
    end
    
    % Z-score at every time window so that coefficient analysis can be done
    % properly at the end (not Z-scoring will mean the coefficients are a
    % function of both the the dynamic range and signal:noise ratio.
    % Z-scoring may be done automatically in the MATLAB LDA implementation,
    % but do it here just to be safe)
    l_hand_holder = l_hand_holder(:,[11:52, 53:88]);
    r_hand_holder = r_hand_holder(:,[11:52, 53:88]);
    mu = mean([l_hand_holder; r_hand_holder]);
    sigma = std([l_hand_holder; r_hand_holder])+1e-4;
    
    l_hand(:,unit,:) = ...
        (l_hand_holder - repmat(mu, [size(l_hand_holder,1),1]))./...
        repmat(sigma, [size(l_hand_holder,1),1]);
    r_hand(:,unit,:) = ...
        (r_hand_holder - repmat(mu, [size(r_hand_holder,1),1]))./...
        repmat(sigma, [size(r_hand_holder,1),1]);
    
end

X = cat(1, l_hand, r_hand);




