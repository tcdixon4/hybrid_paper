function [p_val] = dedicated_mod_stats(unit_data)

% Multi-factorial permutation test for "dedicated" modulation
% Phase(3) X Pref(3)
% Units are crossed with Phase.
% Phase is nested in Pref
% ie each Unit has 3 different Phases that it is recorded during, each
% of which has its own Pref
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
% p_val - p-values for each effect type, organized in a nested data struct
%         (struct: p_val.[effect type].[factors])
% 

%% Prepare arm preference and mod data structures

% isolate only single units and non-repeated units
unit_data = unit_data(([unit_data.unit_type]==1) & ([unit_data.repeat]==false));
num_units = length(unit_data);

% col 1: Unit id
% col 2: Phase -           Rest=1, Instruct=2, Move=3
% col 3: Arm Pref Regime - Ipsi = -1, Neutral = 0, Contra = 1, excl = 9
% col 4: Hand -            Ipsi = 0, Contra = 1
% col 5: Modulation -      continuous depend. variable

for unit = num_units:-1:1
    
    %% unit id
    data_mat((6*unit-5):unit*6,1) = repmat(unit,[6,1]);
    
    
    %% phase
    data_mat((6*unit-5):unit*6,2) = [1;2;3;1;2;3];
    
    
    %% hand
    data_mat((6*unit-5):unit*6,4) = [0;0;0;1;1;1];
    
    
    %% arm pref regime (contra/neutral/ipsi) and modulation
    if unit_data(unit).hem == 0 % L hem
        
        % modulation
        data_mat((6*unit-5),5) = unit_data(unit).modulation.rest(1);
        data_mat((6*unit-4),5) = unit_data(unit).modulation.prep(1);
        data_mat((6*unit-3),5) = unit_data(unit).modulation.move(1);
        data_mat((6*unit-2),5) = unit_data(unit).modulation.rest(2);
        data_mat((6*unit-1),5) = unit_data(unit).modulation.prep(2);
        data_mat((6*unit),5) = unit_data(unit).modulation.move(2);
        
        % Rest arm pref regime
        if unit_data(unit).arm_pref.rest>0.4
            data_mat([(6*unit-5),(6*unit-2)],3) = 1; % contra
        elseif unit_data(unit).arm_pref.rest>-0.3 && ...
                unit_data(unit).arm_pref.rest<0.3
            data_mat([(6*unit-5),(6*unit-2)],3) = 0; % neutral
        elseif unit_data(unit).arm_pref.rest<0.4
            data_mat([(6*unit-5),(6*unit-2)],3) = -1; % ipsi
        else
            data_mat([(6*unit-5),(6*unit-2)],3) = 9; % exclude
        end
        
        % Instruct arm pref regime
        if unit_data(unit).arm_pref.prep>0.4
            data_mat([(6*unit-4),(6*unit-1)],3) = 1; % contra
        elseif unit_data(unit).arm_pref.prep>-0.3 && ...
                unit_data(unit).arm_pref.prep<0.3
            data_mat([(6*unit-4),(6*unit-1)],3) = 0; % neutral
        elseif unit_data(unit).arm_pref.prep<0.4
            data_mat([(6*unit-4),(6*unit-1)],3) = -1; % ipsi
        else
            data_mat([(6*unit-4),(6*unit-1)],3) = 9; % exclude
        end
        
        % Move arm pref regime
        if unit_data(unit).arm_pref.move>0.4
            data_mat([(6*unit-3),(6*unit)],3) = 1; % contra
        elseif unit_data(unit).arm_pref.move>-0.3 && ...
                unit_data(unit).arm_pref.move<0.3
            data_mat([(6*unit-3),(6*unit)],3) = 0; % neutral
        elseif unit_data(unit).arm_pref.move<0.4
            data_mat([(6*unit-3),(6*unit)],3) = -1; % ipsi
        else
            data_mat([(6*unit-1),(6*unit)],3) = 9; % exclude
        end
        
        
    elseif unit_data(unit).hem == 1 % R hem
        
        % modulation
        data_mat((6*unit-5),5) = unit_data(unit).modulation.rest(2);
        data_mat((6*unit-4),5) = unit_data(unit).modulation.prep(2);
        data_mat((6*unit-3),5) = unit_data(unit).modulation.move(2);
        data_mat((6*unit-2),5) = unit_data(unit).modulation.rest(1);
        data_mat((6*unit-1),5) = unit_data(unit).modulation.prep(1);
        data_mat((6*unit),5) = unit_data(unit).modulation.move(1);
        
        % Rest arm pref regime
        if unit_data(unit).arm_pref.rest>0.4
            data_mat([(6*unit-5),(6*unit-2)],3) = -1; % ipsi
        elseif unit_data(unit).arm_pref.rest>-0.3 && ...
                unit_data(unit).arm_pref.rest<0.3
            data_mat([(6*unit-5),(6*unit-2)],3) = 0; % neutral
        elseif unit_data(unit).arm_pref.rest<0.4
            data_mat([(6*unit-5),(6*unit-2)],3) = 1; % contra
        else
            data_mat([(6*unit-5),(6*unit-2)],3) = 9; % exclude
        end
        
        % Instruct arm pref regime
        if unit_data(unit).arm_pref.prep>0.4
            data_mat([(6*unit-4),(6*unit-1)],3) = -1; % ipsi
        elseif unit_data(unit).arm_pref.prep>-0.3 && ...
                unit_data(unit).arm_pref.prep<0.3
            data_mat([(6*unit-4),(6*unit-1)],3) = 0; % neutral
        elseif unit_data(unit).arm_pref.prep<0.4
            data_mat([(6*unit-4),(6*unit-1)],3) = 1; % contra
        else
            data_mat([(6*unit-4),(6*unit-1)],3) = 9; % exclude
        end
        
        % Move arm pref regime
        if unit_data(unit).arm_pref.move>0.4
            data_mat([(6*unit-3),(6*unit)],3) = -1; % ipsi
        elseif unit_data(unit).arm_pref.move>-0.3 && ...
                unit_data(unit).arm_pref.move<0.3
            data_mat([(6*unit-3),(6*unit)],3) = 0; % neutral
        elseif unit_data(unit).arm_pref.move<0.4
            data_mat([(6*unit-3),(6*unit)],3) = 1; % contra
        else
            data_mat([(6*unit-3),(6*unit)],3) = 9; % exclude
        end
        
    end
    
    
end


%% Permutation procedure

for perm = 10000:-1:1
    
    %% Main effect: Hand
    % Is the proportion of contra-dedicated modulation larger than the
    % proportion of ipsi-dedicated modulation, in general?
    %
    % Kind of a weird test because it will be totally dominated by Move
    % phase data since it is way higher modulation, but doing it anyway for
    % thoroughness
    
    % for each unit, flip ipsi and contra with probability p=0.5
    perm_data = data_mat;
    for unit = num_units:-1:1
        if randi(2,1)==1
            perm_data((6*unit-5):(6*unit),3) ...
                = -perm_data((6*unit-5):(6*unit),3);
            perm_data((6*unit-5):(6*unit),4) ...
                = 1 - perm_data((6*unit-5):(6*unit),4);
        end
    end
    
    % calculate proportional modulation of each new grouping
    contra_ded_idx = (perm_data(:,3)==1) & (perm_data(:,4)==1);
    ipsi_ded_idx = (perm_data(:,3)==-1) & (perm_data(:,4)==0);
    contra_full_idx = (perm_data(:,4)==1);
    ipsi_full_idx = (perm_data(:,4)==0);
    
    contra_ded = sum(perm_data(contra_ded_idx,5))/...
        sum(perm_data(contra_full_idx,5));
    ipsi_ded = sum(perm_data(ipsi_ded_idx,5))/...
        sum(perm_data(ipsi_full_idx,5));
    
    % compute contrast
    me_hand_null(perm,1) = abs(contra_ded - ipsi_ded);
    
    
    %% Interaction effect: Hand X Phase
    % Do the differences across phases depend on which hand we consider?
    
    % for each unit, flip ipsi and contra with probability p=0.5
    perm_data = data_mat;
    for unit = num_units:-1:1
        if randi(2,1)==1
            perm_data((6*unit-5):(6*unit),3) ...
                = -perm_data((6*unit-5):(6*unit),3);
            perm_data((6*unit-5):(6*unit),4) ...
                = 1 - perm_data((6*unit-5):(6*unit),4);
        end
    end
    
    % shuffle Phase labels within units
    for unit = num_units:-1:1
        new_phases = randperm(3);
        perm_data((6*unit-5):(6*unit),2) = [new_phases'; new_phases'];
    end
    
    % calculate proportional modulation of each new grouping
    rest_idx = perm_data(:,2)==1;
    prep_idx = perm_data(:,2)==2;
    move_idx = perm_data(:,2)==3;
    contra_ded_idx = (perm_data(:,3)==1) & (perm_data(:,4)==1);
    ipsi_ded_idx = (perm_data(:,3)==-1) & (perm_data(:,4)==0);
    contra_full_idx = (perm_data(:,4)==1);
    ipsi_full_idx = (perm_data(:,4)==0);
    
    contra_rest_ded = sum(perm_data(rest_idx & contra_ded_idx,5))/...
        sum(perm_data(rest_idx & contra_full_idx,5));
    contra_prep_ded = sum(perm_data(prep_idx & contra_ded_idx,5))/...
        sum(perm_data(prep_idx & contra_full_idx,5));
    contra_move_ded = sum(perm_data(move_idx & contra_ded_idx,5))/...
        sum(perm_data(move_idx & contra_full_idx,5));
    ipsi_rest_ded = sum(perm_data(rest_idx & ipsi_ded_idx,5))/...
        sum(perm_data(rest_idx & ipsi_full_idx,5));
    ipsi_prep_ded = sum(perm_data(prep_idx & ipsi_ded_idx,5))/...
        sum(perm_data(prep_idx & ipsi_full_idx,5));
    ipsi_move_ded = sum(perm_data(move_idx & ipsi_ded_idx,5))/...
        sum(perm_data(move_idx & ipsi_full_idx,5));
    
    % compute contrast
    ie_hand_phase_null(perm,1) = ...
        abs(var([contra_rest_ded, contra_prep_ded, contra_move_ded])...
        - var([ipsi_rest_ded, ipsi_prep_ded, ipsi_move_ded]));
    
    
    %% Simple effect: Phase within Hand
    % Does the proportion of Contra-dedicated modulation differ across task
    % phases? Ipsi-dedicated?
    
    % shuffle Phase labels within units
    perm_data = data_mat;
    for unit = num_units:-1:1
        new_phases = randperm(3);
        perm_data((6*unit-5):(6*unit),2) = [new_phases'; new_phases'];
    end
    
    % calculate proportional modulation of each new grouping
    rest_idx = perm_data(:,2)==1;
    prep_idx = perm_data(:,2)==2;
    move_idx = perm_data(:,2)==3;
    contra_ded_idx = (perm_data(:,3)==1) & (perm_data(:,4)==1);
    ipsi_ded_idx = (perm_data(:,3)==-1) & (perm_data(:,4)==0);
    contra_full_idx = (perm_data(:,4)==1);
    ipsi_full_idx = (perm_data(:,4)==0);
    
    contra_rest_ded = sum(perm_data(rest_idx & contra_ded_idx,5))/...
        sum(perm_data(rest_idx & contra_full_idx,5));
    contra_prep_ded = sum(perm_data(prep_idx & contra_ded_idx,5))/...
        sum(perm_data(prep_idx & contra_full_idx,5));
    contra_move_ded = sum(perm_data(move_idx & contra_ded_idx,5))/...
        sum(perm_data(move_idx & contra_full_idx,5));
    ipsi_rest_ded = sum(perm_data(rest_idx & ipsi_ded_idx,5))/...
        sum(perm_data(rest_idx & ipsi_full_idx,5));
    ipsi_prep_ded = sum(perm_data(prep_idx & ipsi_ded_idx,5))/...
        sum(perm_data(prep_idx & ipsi_full_idx,5));
    ipsi_move_ded = sum(perm_data(move_idx & ipsi_ded_idx,5))/...
        sum(perm_data(move_idx & ipsi_full_idx,5));
    
    % compute contrast
    se_phase_null.ipsi(perm,1) = ...
        var([ipsi_rest_ded, ipsi_prep_ded, ipsi_move_ded]);
    se_phase_null.contra(perm,1) = ...
        var([contra_rest_ded, contra_prep_ded, contra_move_ded]);
    
    
    %% Simple effect: Hand within Phase
    % Is the proportion of contra-dedicated modulation larger than the
    % proportion of ipsi-dedicated modulation within each phase?
    
    % for each unit, flip ipsi and contra with probability p=0.5
    perm_data = data_mat;
    for unit = num_units:-1:1
        if randi(2,1)==1
            perm_data((6*unit-5):(6*unit),3) ...
                = -perm_data((6*unit-5):(6*unit),3);
            perm_data((6*unit-5):(6*unit),4) ...
                = 1 - perm_data((6*unit-5):(6*unit),4);
        end
    end
    
    % calculate proportional modulation of each new grouping
    rest_idx = perm_data(:,2)==1;
    prep_idx = perm_data(:,2)==2;
    move_idx = perm_data(:,2)==3;
    contra_ded_idx = (perm_data(:,3)==1) & (perm_data(:,4)==1);
    ipsi_ded_idx = (perm_data(:,3)==-1) & (perm_data(:,4)==0);
    contra_full_idx = (perm_data(:,4)==1);
    ipsi_full_idx = (perm_data(:,4)==0);
    
    contra_rest_ded = sum(perm_data(rest_idx & contra_ded_idx,5))/...
        sum(perm_data(rest_idx & contra_full_idx,5));
    contra_prep_ded = sum(perm_data(prep_idx & contra_ded_idx,5))/...
        sum(perm_data(prep_idx & contra_full_idx,5));
    contra_move_ded = sum(perm_data(move_idx & contra_ded_idx,5))/...
        sum(perm_data(move_idx & contra_full_idx,5));
    ipsi_rest_ded = sum(perm_data(rest_idx & ipsi_ded_idx,5))/...
        sum(perm_data(rest_idx & ipsi_full_idx,5));
    ipsi_prep_ded = sum(perm_data(prep_idx & ipsi_ded_idx,5))/...
        sum(perm_data(prep_idx & ipsi_full_idx,5));
    ipsi_move_ded = sum(perm_data(move_idx & ipsi_ded_idx,5))/...
        sum(perm_data(move_idx & ipsi_full_idx,5));
    
    % compute contrast
    se_hand_null.rest(perm,1) = ...
        abs(contra_rest_ded - ipsi_rest_ded);
    se_hand_null.prep(perm,1) = ...
        abs(contra_prep_ded - ipsi_prep_ded);
    se_hand_null.move(perm,1) = ...
        abs(contra_move_ded - ipsi_move_ded);
    
    
end


%% Main effect: Hand p-value

% calculate proportional modulation of each grouping
contra_ded_idx = (data_mat(:,3)==1) & (data_mat(:,4)==1);
ipsi_ded_idx = (data_mat(:,3)==-1) & (data_mat(:,4)==0);
contra_full_idx = (data_mat(:,4)==1);
ipsi_full_idx = (data_mat(:,4)==0);

contra_ded = sum(data_mat(contra_ded_idx,5))/...
    sum(data_mat(contra_full_idx,5));
ipsi_ded = sum(data_mat(ipsi_ded_idx,5))/...
    sum(data_mat(ipsi_full_idx,5));

% compute observed contrast and p-value
me_hand_obs = abs(contra_ded - ipsi_ded);
p_val.me.hand = mean(me_hand_null>me_hand_obs);


%% Interaction effect: Hand X Phase p-value

% calculate proportional modulation of each grouping
rest_idx = data_mat(:,2)==1;
prep_idx = data_mat(:,2)==2;
move_idx = data_mat(:,2)==3;

contra_rest_ded = sum(data_mat(rest_idx & contra_ded_idx,5))/...
    sum(data_mat(rest_idx & contra_full_idx,5));
contra_prep_ded = sum(data_mat(prep_idx & contra_ded_idx,5))/...
    sum(data_mat(prep_idx & contra_full_idx,5));
contra_move_ded = sum(data_mat(move_idx & contra_ded_idx,5))/...
    sum(data_mat(move_idx & contra_full_idx,5));
ipsi_rest_ded = sum(data_mat(rest_idx & ipsi_ded_idx,5))/...
    sum(data_mat(rest_idx & ipsi_full_idx,5));
ipsi_prep_ded = sum(data_mat(prep_idx & ipsi_ded_idx,5))/...
    sum(data_mat(prep_idx & ipsi_full_idx,5));
ipsi_move_ded = sum(data_mat(move_idx & ipsi_ded_idx,5))/...
    sum(data_mat(move_idx & ipsi_full_idx,5));

% compute contrast and p-value
ie_hand_phase_obs = ...
    abs(var([contra_rest_ded, contra_prep_ded, contra_move_ded])...
    - var([ipsi_rest_ded, ipsi_prep_ded, ipsi_move_ded]));
p_val.ie.hand_phase = mean(ie_hand_phase_null>ie_hand_phase_obs);


%% Simple effect: Phase within Hand p-value

% compute contrast and p-value
se_phase_obs.ipsi = ...
    var([ipsi_rest_ded, ipsi_prep_ded, ipsi_move_ded]);
se_phase_obs.contra = ...
    var([contra_rest_ded, contra_prep_ded, contra_move_ded]);
p_val.se.phase.ipsi = mean(se_phase_null.ipsi>se_phase_obs.ipsi);
p_val.se.phase.contra = mean(se_phase_null.contra>se_phase_obs.contra);


%% Simple effect: Hand within Phase p-value

% compute contrast and p-value
se_hand_obs.rest = ...
    abs(contra_rest_ded - ipsi_rest_ded);
se_hand_obs.prep = ...
    abs(contra_prep_ded - ipsi_prep_ded);
se_hand_obs.move = ...
    abs(contra_move_ded - ipsi_move_ded);
p_val.se.hand.rest = mean(se_hand_null.rest>se_hand_obs.rest);
p_val.se.hand.prep = mean(se_hand_null.prep>se_hand_obs.prep);
p_val.se.hand.move = mean(se_hand_null.move>se_hand_obs.move);






    