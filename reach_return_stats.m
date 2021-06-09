function [p_val] = reach_return_stats(unit_data)

%
% Statistical testing for control analysis using concatenated reach and
% return data. This function first tests for significant non-zero slopes in 
% the regression lines fit to arm preference and modulation data, independently for ipsi- and contra-preferring units.
% It then tests whether the proportion of modulation represented within 
% arm-dedicated regimes of the arm-preference spectrum differs across 
% phases after extend phase data is replaced with Reach and Return. Again, 
% this is done independently for ipsi- and contra-dedicated regimes. 
% All of the testing here uses permutation approaches.
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

%col 1: Unit id
%col 2: Area -            Pmd=0, M1 = 1
%col 3: Phase -           Rest=1, Instruct=2, Reach&Return=3
%col 4: Preferred arm -   Ipsi=0, Contra=1
%col 5: Arm Pref Regime - Ipsi = -1, Neutral = 0, Contra = 1, excl = 9
%col 6: Arm Pref -        continuous depend. variable
%col 7: Hand -            Ipsi = 0, Contra = 1
%col 8: Modulation -      continuous depend. variable

for unit = num_units:-1:1
    %% unit id
    data_mat((6*unit-5):unit*6,1) = repmat(unit,[6,1]);
    
    %% area
    data_mat((6*unit-5):unit*6,2) = repmat(unit_data(unit).area,[6,1]);
    
    %% phase
    data_mat((6*unit-5):unit*6,3) = [1;2;3;1;2;3];  
    
    %% hand
    data_mat((6*unit-5):unit*6,7) = [0;0;0;1;1;1];
    
    %% arm pref (continuous depend. variable)
    data_mat((6*unit-[2,5]),6) = abs(unit_data(unit).arm_pref.rest);
    data_mat((6*unit-[1,4]),6) = abs(unit_data(unit).arm_pref.prep);
    data_mat(6*unit-[0,3],6) = abs(unit_data(unit).arm_pref.extend);   
    
    %% preferred arm, arm pref regime (contra/neutral/ipsi), and modulation
    if unit_data(unit).hem == 0 % L hem
        
        % preferred arm (ipsi/contra)
        data_mat(6*unit-[2,5],4) = unit_data(unit).arm_pref.rest>0;
        data_mat(6*unit-[1,4],4) = unit_data(unit).arm_pref.prep>0;
        data_mat(6*unit-[0,3],4) = unit_data(unit).arm_pref.extend>0;
        
        % modulation
        data_mat((6*unit-5),8) = unit_data(unit).modulation.rest(1);
        data_mat((6*unit-4),8) = unit_data(unit).modulation.prep(1);
        data_mat((6*unit-3),8) = unit_data(unit).modulation.extend(1);
        data_mat((6*unit-2),8) = unit_data(unit).modulation.rest(2);
        data_mat((6*unit-1),8) = unit_data(unit).modulation.prep(2);
        data_mat((6*unit),8) = unit_data(unit).modulation.extend(2);
        
        % Rest arm pref regime
        if unit_data(unit).arm_pref.rest>0.4
            data_mat(6*unit-[2,5],5) = 1; % contra
        elseif unit_data(unit).arm_pref.rest>-0.3 && ...
                unit_data(unit).arm_pref.rest<0.3
            data_mat(6*unit-[2,5],5) = 0; % neutral
        elseif unit_data(unit).arm_pref.rest<0.4
            data_mat(6*unit-[2,5],5) = -1; % ipsi
        else
            data_mat(6*unit-[2,5],5) = 9; % exclude
        end
        
        % Instruct arm pref regime
        if unit_data(unit).arm_pref.prep>0.4
            data_mat(6*unit-[1,4],5) = 1; % contra
        elseif unit_data(unit).arm_pref.prep>-0.3 && ...
                unit_data(unit).arm_pref.prep<0.3
            data_mat(6*unit-[1,4],5) = 0; % neutral
        elseif unit_data(unit).arm_pref.prep<0.4
            data_mat(6*unit-[1,4],5) = -1; % ipsi
        else
            data_mat(6*unit-[1,4],5) = 9; % exclude
        end
        
        % Reach&Return arm pref regime
        if unit_data(unit).arm_pref.extend>0.4
            data_mat(6*unit-[0,3],5) = 1; % contra
        elseif unit_data(unit).arm_pref.extend>-0.3 && ...
                unit_data(unit).arm_pref.extend<0.3
            data_mat(6*unit-[0,3],5) = 0; % neutral
        elseif unit_data(unit).arm_pref.extend<0.4
            data_mat(6*unit-[0,3],5) = -1; % ipsi
        else
            data_mat(6*unit-[0,3],5) = 9; % exclude
        end
        
        
    elseif unit_data(unit).hem == 1 % R hem
        
        % preferred arm (ipsi/contra)
        data_mat(6*unit-[2,5],4) = unit_data(unit).arm_pref.rest<0;
        data_mat(6*unit-[1,4],4) = unit_data(unit).arm_pref.prep<0;
        data_mat(6*unit-[0,3],4) = unit_data(unit).arm_pref.extend<0;
        
        % modulation
        data_mat((6*unit-5),8) = unit_data(unit).modulation.rest(2);
        data_mat((6*unit-4),8) = unit_data(unit).modulation.prep(2);
        data_mat((6*unit-3),8) = unit_data(unit).modulation.extend(2);
        data_mat((6*unit-2),8) = unit_data(unit).modulation.rest(1);
        data_mat((6*unit-1),8) = unit_data(unit).modulation.prep(1);
        data_mat((6*unit),8) = unit_data(unit).modulation.extend(1);
        
        % Rest arm pref regime
        if unit_data(unit).arm_pref.rest>0.4
            data_mat(6*unit-[2,5],5) = -1; % ipsi
        elseif unit_data(unit).arm_pref.rest>-0.3 && ...
                unit_data(unit).arm_pref.rest<0.3
            data_mat(6*unit-[2,5],5) = 0; % neutral
        elseif unit_data(unit).arm_pref.rest<0.4
            data_mat(6*unit-[2,5],5) = 1; % contra
        else
            data_mat(6*unit-[2,5],5) = 9; % exclude
        end
        
        % Instruct arm pref regime
        if unit_data(unit).arm_pref.prep>0.4
            data_mat(6*unit-[1,4],5) = -1; % ipsi
        elseif unit_data(unit).arm_pref.prep>-0.3 && ...
                unit_data(unit).arm_pref.prep<0.3
            data_mat(6*unit-[1,4],5) = 0; % neutral
        elseif unit_data(unit).arm_pref.prep<0.4
            data_mat(6*unit-[1,4],5) = 1; % contra
        else
            data_mat(6*unit-[1,4],5) = 9; % exclude
        end
        
        % Reach&Return arm pref regime
        if unit_data(unit).arm_pref.extend>0.4
            data_mat(6*unit-[0,3],5) = -1; % ipsi
        elseif unit_data(unit).arm_pref.extend>-0.3 && ...
                unit_data(unit).arm_pref.extend<0.3
            data_mat(6*unit-[0,3],5) = 0; % neutral
        elseif unit_data(unit).arm_pref.extend<0.4
            data_mat(6*unit-[0,3],5) = 1; % contra
        else
            data_mat(6*unit-[0,3],5) = 9; % exclude
        end
        
    end
    
    
end


%% Permutation procedure

for perm = 10000:-1:1
    
    %% Construct permuted null distribution of regression slopes for 
    %  arm pref & modulation relationship
    
    % separate ipsi-pref units from contra-pref units
    ipsi_data = data_mat(data_mat(:,4)==0,:);
    contra_data = data_mat(data_mat(:,4)==1,:);
    % isolate Reach&Return data
    ipsi_data = ipsi_data(ipsi_data(:,3)==3,:);
    contra_data = contra_data(contra_data(:,3)==3,:);
    % isolate preferred arm data
    ipsi_data = ipsi_data(ipsi_data(:,7)==0,:);
    contra_data = contra_data(contra_data(:,7)==1,:);
    % separate PMd units from M1 units
    ipsi_pmd_data = ipsi_data(ipsi_data(:,2)==0,:);
    contra_pmd_data = contra_data(contra_data(:,2)==0,:);
    ipsi_m1_data = ipsi_data(ipsi_data(:,2)==1,:);
    contra_m1_data = contra_data(contra_data(:,2)==1,:);
    
    % randomize arm pref and modulation pairings and compute a new slope
    % for constructing the null distributions    
    b_ipsi_pmd_null(perm,1) = ...
        compute_slope(ipsi_pmd_data(:,6),ipsi_pmd_data(:,8), true);
    b_contra_pmd_null(perm,1) = ...
        compute_slope(contra_pmd_data(:,6),contra_pmd_data(:,8), true);
    b_ipsi_m1_null(perm,1) = ...
        compute_slope(ipsi_m1_data(:,6),ipsi_m1_data(:,8), true);
    b_contra_m1_null(perm,1) = ...
        compute_slope(contra_m1_data(:,6),contra_m1_data(:,8), true);
    
    
   %% Construct permuted null distribution of between-phase contrast in the 
   %  proportion of dedicated variance
   
   % shuffle Phase labels within units
    perm_data = data_mat;
    for unit = num_units:-1:1
        new_phases = randperm(3);
        perm_data((6*unit-5):(6*unit),3) = [new_phases'; new_phases'];
    end
    
    % calculate proportional modulation of each new grouping
    rest_idx = perm_data(:,3)==1;
    prep_idx = perm_data(:,3)==2;
    extend_idx = perm_data(:,3)==3;
    contra_ded_idx = (perm_data(:,5)==1) & (perm_data(:,7)==1);
    ipsi_ded_idx = (perm_data(:,5)==-1) & (perm_data(:,7)==0);
    contra_full_idx = (perm_data(:,7)==1);
    ipsi_full_idx = (perm_data(:,7)==0);
    
    contra_rest_ded = sum(perm_data(rest_idx & contra_ded_idx,8))/...
        sum(perm_data(rest_idx & contra_full_idx,8));
    contra_prep_ded = sum(perm_data(prep_idx & contra_ded_idx,8))/...
        sum(perm_data(prep_idx & contra_full_idx,8));
    contra_extend_ded = sum(perm_data(extend_idx & contra_ded_idx,8))/...
        sum(perm_data(extend_idx & contra_full_idx,8));
    ipsi_rest_ded = sum(perm_data(rest_idx & ipsi_ded_idx,8))/...
        sum(perm_data(rest_idx & ipsi_full_idx,8));
    ipsi_prep_ded = sum(perm_data(prep_idx & ipsi_ded_idx,8))/...
        sum(perm_data(prep_idx & ipsi_full_idx,8));
    ipsi_extend_ded = sum(perm_data(extend_idx & ipsi_ded_idx,8))/...
        sum(perm_data(extend_idx & ipsi_full_idx,8));
    
    % compute contrast
    se_phase_null.ipsi(perm,1) = ...
        var([ipsi_rest_ded, ipsi_prep_ded, ipsi_extend_ded]);
    se_phase_null.contra(perm,1) = ...
        var([contra_rest_ded, contra_prep_ded, contra_extend_ded]);
    
    
end


%% p-value for arm pref vs modulation regression slopes

% compute true slopes
b_ipsi_pmd_true = ...
        compute_slope(ipsi_pmd_data(:,6),ipsi_pmd_data(:,8), false);
b_contra_pmd_true = ...
        compute_slope(contra_pmd_data(:,6),contra_pmd_data(:,8), false);
b_ipsi_m1_true = ...
        compute_slope(ipsi_m1_data(:,6),ipsi_m1_data(:,8), false);
b_contra_m1_true = ...
        compute_slope(contra_m1_data(:,6),contra_m1_data(:,8), false);
    
% compute p-value as proportion of null distribution at least as extreme
p_val.slope.ipsi.pmd = mean(abs(b_ipsi_pmd_null)>=b_ipsi_pmd_true);
p_val.slope.contra.pmd = mean(abs(b_contra_pmd_null)>=b_contra_pmd_true);
p_val.slope.ipsi.m1 = mean(abs(b_ipsi_m1_null)>=b_ipsi_m1_true);
p_val.slope.contra.m1 = mean(abs(b_contra_m1_null)>=b_contra_m1_true);

   
%% p-value for simple effects of phase

% calculate proportional modulation of each original grouping
rest_idx = data_mat(:,3)==1;
prep_idx = data_mat(:,3)==2;
extend_idx = data_mat(:,3)==3;
contra_ded_idx = (data_mat(:,5)==1) & (data_mat(:,7)==1);
ipsi_ded_idx = (data_mat(:,5)==-1) & (data_mat(:,7)==0);
contra_full_idx = (data_mat(:,7)==1);
ipsi_full_idx = (data_mat(:,7)==0);

contra_rest_ded = sum(data_mat(rest_idx & contra_ded_idx,8))/...
    sum(data_mat(rest_idx & contra_full_idx,8));
contra_prep_ded = sum(data_mat(prep_idx & contra_ded_idx,8))/...
    sum(data_mat(prep_idx & contra_full_idx,8));
contra_extend_ded = sum(data_mat(extend_idx & contra_ded_idx,8))/...
    sum(data_mat(extend_idx & contra_full_idx,8));
ipsi_rest_ded = sum(data_mat(rest_idx & ipsi_ded_idx,8))/...
    sum(data_mat(rest_idx & ipsi_full_idx,8));
ipsi_prep_ded = sum(data_mat(prep_idx & ipsi_ded_idx,8))/...
    sum(data_mat(prep_idx & ipsi_full_idx,8));
ipsi_extend_ded = sum(data_mat(extend_idx & ipsi_ded_idx,8))/...
    sum(data_mat(extend_idx & ipsi_full_idx,8));

% compute true contrast
se_phase_true.ipsi = ...
    var([ipsi_rest_ded, ipsi_prep_ded, ipsi_extend_ded]);
se_phase_true.contra = ...
    var([contra_rest_ded, contra_prep_ded, contra_extend_ded]);

% compute p-value as proportion of null distribution at least as extreme
p_val.phase.ipsi = mean(se_phase_null.ipsi>=se_phase_true.ipsi);
p_val.phase.contra = mean(se_phase_null.contra>=se_phase_true.contra);




%% nested functions

    function b = compute_slope(X,Y, is_perm)
        X = abs(X);
        Y = log10(Y);
        if is_perm
            X = X(randperm(length(X)));
        end
        rmv = find(Y<-1);
        X(rmv) = [];
        Y(rmv) = [];
        b = polyfit(X,Y,1);
        b = b(1);
    end

end
    
    