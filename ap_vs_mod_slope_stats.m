function [p_val] = ap_vs_mod_slope_stats(unit_data)

%
% Multi-factorial permutation test for arm preference and modulation
% strength relationships
% Area(2) X Phase(3) X Pref(2)
% Units are nested in Area and crossed with Phase.
% Phase is nested in Pref
% ie each Unit is located in a single Area and has 3 different Phases that
% it is recorded during, each of which has its own Pref
%
% This also contains a permutation test for positive regression slope
% relating arm preference and modulation in the Reach & Return data. This 
% p-value is logged as p_val.se.unit_label.extend
% Unit labels are permuted for the arm preference while modulation is held
% constant. Note that we are NOT permuting unit labels on the modulation
% values for each arm and re-computing the arm preference to construct the
% null slope samples. This means that the slope is impacted by the marginal 
% distributions of modulation values, which is done intentionally. We are
% not interested in assessing the relationship independent of the marginal
% distributions, as that is part of the phenomenon we are studying.
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
% unit_data = unit_data([unit_data.unit_type]==1);
num_units = length(unit_data);

% data mat for multi-factorial permutation testing
% col 1: Unit id
% col 2: Area -            Pmd=0, M1 = 1
% col 3: Phase -           Rest=1, Instruct=2, Move=3
% col 4: Preferred arm -   Ipsi=0, Contra=1
% col 5: Arm Pref -        continuous depend. variable
% col 6: Modulation -      continuous depend. variable

for unit = num_units:-1:1  
    % unit id
    data_mat((3*unit-2):unit*3,1) = repmat(unit,[3,1]);
    % area
    data_mat((3*unit-2):unit*3,2) = repmat(unit_data(unit).area,[3,1]);
    % phase
    data_mat((3*unit-2):unit*3,3) = [1;2;3];
    % arm pref (continuous depend. variable)
    data_mat((3*unit-2),5) = abs(unit_data(unit).arm_pref.rest);
    data_mat((3*unit-1),5) = abs(unit_data(unit).arm_pref.prep);
    data_mat(3*unit,5) = abs(unit_data(unit).arm_pref.move);
    extend_ap(unit) = abs(unit_data(unit).arm_pref.extend);
    extend_area(unit) = unit_data(unit).area;
    
    if unit_data(unit).hem == 0
        % preferred arm (ipsi/contra)
        data_mat((3*unit-2),4) = unit_data(unit).arm_pref.rest>0;
        data_mat((3*unit-1),4) = unit_data(unit).arm_pref.prep>0;
        data_mat(3*unit,4) = unit_data(unit).arm_pref.move>0;
        extend_pref(unit) = unit_data(unit).arm_pref.extend>0;
        % modulation of preferred arm
        data_mat((3*unit-2),6) = log10(...
            unit_data(unit).modulation.rest(data_mat((3*unit-2),4)+1));
        data_mat((3*unit-1),6) = log10(...
            unit_data(unit).modulation.prep(data_mat((3*unit-1),4)+1));
        data_mat(3*unit,6) = log10(...
            unit_data(unit).modulation.move(data_mat(3*unit,4)+1));  
        extend_mod(unit) = log10(...
            unit_data(unit).modulation.extend(extend_pref(unit)+1));
        
    elseif unit_data(unit).hem == 1
        % preferred arm (ipsi/contra)
        data_mat((3*unit-2),4) = unit_data(unit).arm_pref.rest<0;
        data_mat((3*unit-1),4) = unit_data(unit).arm_pref.prep<0;
        data_mat(3*unit,4) = unit_data(unit).arm_pref.move<0;
        extend_pref(unit) = unit_data(unit).arm_pref.extend<0;
        % modulation of preferred arm
        data_mat((3*unit-2),6) = log10(...
            unit_data(unit).modulation.rest(2-data_mat((3*unit-2),4)));
        data_mat((3*unit-1),6) = log10(...
            unit_data(unit).modulation.prep(2-data_mat((3*unit-1),4)));
        data_mat(3*unit,6) = log10(...
            unit_data(unit).modulation.move(2-data_mat(3*unit,4))); 
        extend_mod(unit) = log10(...
            unit_data(unit).modulation.extend(2-extend_pref(unit)));
        
    end    
    
end

% separate Reach & Return data into PMd/M1 and Ipsi-pref/Contra-pref
ipsipref_m1_ap = extend_ap((~extend_pref) & (extend_area));
ipsipref_pmd_ap = extend_ap((~extend_pref) & (~extend_area));
contrapref_m1_ap = extend_ap((extend_pref) & (extend_area));
contrapref_pmd_ap = extend_ap((extend_pref) & (~extend_area));

ipsipref_m1_mod = extend_mod((~extend_pref) & (extend_area));
ipsipref_pmd_mod = extend_mod((~extend_pref) & (~extend_area));
contrapref_m1_mod = extend_mod((extend_pref) & (extend_area));
contrapref_pmd_mod = extend_mod((extend_pref) & (~extend_area));


%% Permutation procedure

for perm = 10000:-1:1
    
    %% Main effect: Area
    % shuffle Area labels across units
    perm_data = data_mat;
    shuff_area = perm_data(1:3:end,2);
    shuff_area = shuff_area(randperm(num_units))';
    shuff_area = repmat(shuff_area,[3,1]);
    perm_data(:,2) = shuff_area(:);
    % remove negative outliers (log scale can inflate these)
    rmv = perm_data(:,6)<-1;
    perm_data(rmv,:) = [];
    % calculate slope of each new grouping
    pmd_idx = find(perm_data(:,2)==0);
    m1_idx = find(perm_data(:,2)==1);
    b_pmd = polyfit(perm_data(pmd_idx,5),perm_data(pmd_idx,6),1);
    b_m1 = polyfit(perm_data(m1_idx,5),perm_data(m1_idx,6),1);
    % compute contrast
    me_area_null(perm,1) = abs(b_m1(1) - b_pmd(1));
    
    
    %% Main effect: Phase
    % shuffle Phase labels within units
    perm_data = data_mat;
    for unit = num_units:-1:1
        perm_data((3*unit-2):unit*3,3) = randperm(3);
    end
    % remove negative outliers (log scale can inflate these)
    rmv = perm_data(:,6)<-1;
    perm_data(rmv,:) = [];
    % calculate slope of each new grouping
    rest_idx = find(perm_data(:,3)==1);
    prep_idx = find(perm_data(:,3)==2);
    move_idx = find(perm_data(:,3)==3);
    b_rest = polyfit(perm_data(rest_idx,5),perm_data(rest_idx,6),1);
    b_prep = polyfit(perm_data(prep_idx,5),perm_data(prep_idx,6),1);
    b_move = polyfit(perm_data(move_idx,5),perm_data(move_idx,6),1);
    % compute contrast
    me_phase_null(perm,1) = var([b_rest(1), b_prep(1), b_move(1)]);
    
    
    %% Main effect: Pref
    % shuffle Pref labels within levels of Phase
    phase_sorted_idx = [find(data_mat(:,3)==1),...
        find(data_mat(:,3)==2),find(data_mat(:,3)==3)];
    perm_data = data_mat(phase_sorted_idx,:);
    perm_data(1:num_units,4) = ...
        perm_data(randperm(num_units),4);
    perm_data((num_units+1):(num_units*2), 4) = ...
        perm_data(randperm(num_units)+num_units, 4);
    perm_data((2*num_units+1):(num_units*3), 4) = ...
        perm_data(randperm(num_units)+2*num_units, 4);
    % remove negative outliers (log scale can inflate these)
    rmv = perm_data(:,6)<-1;
    perm_data(rmv,:) = [];
    % calculate slope of each new grouping
    ipsi_idx = find(perm_data(:,4)==0);
    contra_idx = find(perm_data(:,4)==1);
    b_ipsi = polyfit(perm_data(ipsi_idx,5),perm_data(ipsi_idx,6),1);
    b_contra = polyfit(perm_data(contra_idx,5),perm_data(contra_idx,6),1);
    % compute contrast
    me_pref_null(perm,1) = abs(b_contra(1) - b_ipsi(1));
    
    
    %% Interaction effect: Area X Phase (all phases)
    % shuffle Area labels across units
    perm_data = data_mat;
    shuff_area = perm_data(1:3:end,2);
    shuff_area = shuff_area(randperm(num_units))';
    shuff_area = repmat(shuff_area,[3,1]);
    perm_data(:,2) = shuff_area(:);
    % shuffle Phase labels within units
    for unit = num_units:-1:1
        perm_data((3*unit-2):unit*3,3) = randperm(3);
    end
    % remove negative outliers (log scale can inflate these)
    rmv = perm_data(:,6)<-1;
    perm_data(rmv,:) = [];
    % calculate slope of each new grouping
    pmd_rest_idx = find(perm_data(:,2)==0 & perm_data(:,3)==1);
    pmd_prep_idx = find(perm_data(:,2)==0 & perm_data(:,3)==2);
    pmd_move_idx = find(perm_data(:,2)==0 & perm_data(:,3)==3);
    m1_rest_idx = find(perm_data(:,2)==1 & perm_data(:,3)==1);
    m1_prep_idx = find(perm_data(:,2)==1 & perm_data(:,3)==2);
    m1_move_idx = find(perm_data(:,2)==1 & perm_data(:,3)==3);
    b_pmd_rest = ...
        polyfit(perm_data(pmd_rest_idx,5),perm_data(pmd_rest_idx,6),1);
    b_pmd_prep = ...
        polyfit(perm_data(pmd_prep_idx,5),perm_data(pmd_prep_idx,6),1);
    b_pmd_move = ...
        polyfit(perm_data(pmd_move_idx,5),perm_data(pmd_move_idx,6),1);
    b_m1_rest = ...
        polyfit(perm_data(m1_rest_idx,5),perm_data(m1_rest_idx,6),1);
    b_m1_prep = ...
        polyfit(perm_data(m1_prep_idx,5),perm_data(m1_prep_idx,6),1);
    b_m1_move = ...
        polyfit(perm_data(m1_move_idx,5),perm_data(m1_move_idx,6),1);
    % compute contrast
    ie_area_phase_null(perm,1) = abs(...
        var([b_m1_rest(1),b_m1_prep(1),b_m1_move(1)]) - ...
        var([b_pmd_rest(1),b_pmd_prep(1),b_pmd_move(1)]) );
    
    
    %% Interaction effect: Area X Phase (Instruct, Move)
    % shuffle Area labels across units
    perm_data = data_mat;
    shuff_area = perm_data(1:3:end,2);
    shuff_area = shuff_area(randperm(num_units))';
    shuff_area = repmat(shuff_area,[3,1]);
    perm_data(:,2) = shuff_area(:);
    % isolate Instruct and Move, then shuffle Phase labels within units
    rmv = perm_data(:,3)==1;
    perm_data(rmv,:) = [];
    for unit = num_units:-1:1
        perm_data((2*unit-1):unit*2,3) = randperm(2);
    end
    % remove negative outliers (log scale can inflate these)
    rmv = perm_data(:,6)<-1;
    perm_data(rmv,:) = [];
    % calculate slope of each new grouping
    pmd_prep_idx = find(perm_data(:,2)==0 & perm_data(:,3)==2);
    pmd_move_idx = find(perm_data(:,2)==0 & perm_data(:,3)==3);
    m1_prep_idx = find(perm_data(:,2)==1 & perm_data(:,3)==2);
    m1_move_idx = find(perm_data(:,2)==1 & perm_data(:,3)==3);
    b_pmd_prep = ...
        polyfit(perm_data(pmd_prep_idx,5),perm_data(pmd_prep_idx,6),1);
    b_pmd_move = ...
        polyfit(perm_data(pmd_move_idx,5),perm_data(pmd_move_idx,6),1);
    b_m1_prep = ...
        polyfit(perm_data(m1_prep_idx,5),perm_data(m1_prep_idx,6),1);
    b_m1_move = ...
        polyfit(perm_data(m1_move_idx,5),perm_data(m1_move_idx,6),1);
    % compute contrast
    ie_area_phaseIM_null(perm,1) = abs(...
        abs(b_pmd_prep(1)-b_m1_prep(1)) - ...
        abs(b_pmd_move(1)-b_m1_move(1)) );
    
    
    %% Interaction effect: Area X Pref
    % shuffle Area labels across units
    perm_data = data_mat;
    shuff_area = perm_data(1:3:end,2);
    shuff_area = shuff_area(randperm(num_units))';
    shuff_area = repmat(shuff_area,[3,1]);
    perm_data(:,2) = shuff_area(:);
    % shuffle Pref labels within levels of Phase
    phase_sorted_idx = [find(perm_data(:,3)==1),...
        find(perm_data(:,3)==2),find(perm_data(:,3)==3)];
    perm_data = perm_data(phase_sorted_idx,:);
    perm_data(1:num_units,4) = ...
        perm_data(randperm(num_units),4);
    perm_data((num_units+1):(num_units*2), 4) = ...
        perm_data(randperm(num_units)+num_units, 4);
    perm_data((2*num_units+1):(num_units*3), 4) = ...
        perm_data(randperm(num_units)+2*num_units, 4);
    % remove negative outliers (log scale can inflate these)
    rmv = perm_data(:,6)<-1;
    perm_data(rmv,:) = [];
    % calculate slope of each new grouping
    pmd_ipsi_idx = find(perm_data(:,2)==0 & perm_data(:,4)==0);
    pmd_contra_idx = find(perm_data(:,2)==0 & perm_data(:,4)==1);
    m1_ipsi_idx = find(perm_data(:,2)==1 & perm_data(:,4)==0);
    m1_contra_idx = find(perm_data(:,2)==1 & perm_data(:,4)==1);
    b_pmd_ipsi = ...
        polyfit(perm_data(pmd_ipsi_idx,5),perm_data(pmd_ipsi_idx,6),1);
    b_pmd_contra = ...
        polyfit(perm_data(pmd_contra_idx,5),perm_data(pmd_contra_idx,6),1);
    b_m1_ipsi = ...
        polyfit(perm_data(m1_ipsi_idx,5),perm_data(m1_ipsi_idx,6),1);
    b_m1_contra = ...
        polyfit(perm_data(m1_contra_idx,5),perm_data(m1_contra_idx,6),1);
    % compute contrast
    ie_area_pref_null(perm,1) = abs(...
        (b_m1_contra(1) - b_m1_ipsi(1)) - ...
        (b_pmd_contra(1) - b_pmd_ipsi(1)) );
 
    
    %% Interaction effect: Phase X Pref
    % shuffle Phase labels within units
    % keep Prefs w/ original Phases so that proportions of each Pref will  
    % be maintained from their original eg Rest:50%Ipsi, Prep:40%, Move:20%
    perm_data = data_mat;
    for unit = num_units:-1:1
        perm_idx = randperm(3);
        perm_data((3*unit-2):unit*3,3) = perm_idx;
        unit_prefs = perm_data((3*unit-2):unit*3,4);
        perm_data((3*unit-2):unit*3,4) = unit_prefs(perm_idx);
    end
    % shuffle Pref labels within levels of Phase (this may be unnecessary)
    phase_sorted_idx = [find(perm_data(:,3)==1),...
        find(perm_data(:,3)==2),find(perm_data(:,3)==3)];
    perm_data = perm_data(phase_sorted_idx,:);
    perm_data(1:num_units,4) = ...
        perm_data(randperm(num_units),4);
    perm_data((num_units+1):(num_units*2), 4) = ...
        perm_data(randperm(num_units)+num_units, 4);
    perm_data((2*num_units+1):(num_units*3), 4) = ...
        perm_data(randperm(num_units)+2*num_units, 4);
    % remove negative outliers (log scale can inflate these)
    rmv = perm_data(:,6)<-1;
    perm_data(rmv,:) = [];
    % calculate slope of each new grouping
    rest_ipsi_idx = find(perm_data(:,3)==1 & perm_data(:,4)==0);
    prep_ipsi_idx = find(perm_data(:,3)==2 & perm_data(:,4)==0);
    move_ipsi_idx = find(perm_data(:,3)==3 & perm_data(:,4)==0);
    rest_contra_idx = find(perm_data(:,3)==1 & perm_data(:,4)==1);
    prep_contra_idx = find(perm_data(:,3)==2 & perm_data(:,4)==1);
    move_contra_idx = find(perm_data(:,3)==3 & perm_data(:,4)==1);
    b_rest_ipsi = ...
        polyfit(perm_data(rest_ipsi_idx,5),perm_data(rest_ipsi_idx,6),1);
    b_prep_ipsi = ...
        polyfit(perm_data(prep_ipsi_idx,5),perm_data(prep_ipsi_idx,6),1);
    b_move_ipsi = ...
        polyfit(perm_data(move_ipsi_idx,5),perm_data(move_ipsi_idx,6),1);
    b_rest_contra = ...
        polyfit(perm_data(rest_contra_idx,5),perm_data(rest_contra_idx,6),1);
    b_prep_contra = ...
        polyfit(perm_data(prep_contra_idx,5),perm_data(prep_contra_idx,6),1);
    b_move_contra = ...
        polyfit(perm_data(move_contra_idx,5),perm_data(move_contra_idx,6),1);
    % compute contrast
    ie_phase_pref_null(perm,1) = abs(...
        var([b_rest_contra(1), b_prep_contra(1), b_move_contra(1)]) - ...
        var([b_rest_ipsi(1), b_prep_ipsi(1), b_move_ipsi(1)]) );
    
    
    %% Simple effect: Phase within Pref
    % shuffle Phase labels across units for each Pref individually, since
    % separating Prefs no longer allows shuffling Phase labels within units
    ipsi_data = data_mat(data_mat(:,4)==0,:);
    ipsi_data(:,3) = ipsi_data(randperm(size(ipsi_data,1)),3);
    contra_data = data_mat(data_mat(:,4)==1,:);
    contra_data(:,3) = contra_data(randperm(size(contra_data,1)),3);
    % remove negative outliers (log scale can inflate these)
    rmv = ipsi_data(:,6)<-1;
    ipsi_data(rmv,:) = [];
    rmv = contra_data(:,6)<-1;
    contra_data(rmv,:) = [];
    % calculate slope of each new grouping for Ipsi
    rest_idx = find(ipsi_data(:,3)==1);
    prep_idx = find(ipsi_data(:,3)==2);
    move_idx = find(ipsi_data(:,3)==3);
    b_rest = polyfit(ipsi_data(rest_idx,5),ipsi_data(rest_idx,6),1);
    b_prep = polyfit(ipsi_data(prep_idx,5),ipsi_data(prep_idx,6),1);
    b_move = polyfit(ipsi_data(move_idx,5),ipsi_data(move_idx,6),1);
    % compute contrast for contra
    se_phase_ipsi_null(perm,1) = var([b_rest(1), b_prep(1), b_move(1)]);
    % calculate slope of each new grouping for Contra
    rest_idx = find(contra_data(:,3)==1);
    prep_idx = find(contra_data(:,3)==2);
    move_idx = find(contra_data(:,3)==3);
    b_rest = polyfit(contra_data(rest_idx,5),contra_data(rest_idx,6),1);
    b_prep = polyfit(contra_data(prep_idx,5),contra_data(prep_idx,6),1);
    b_move = polyfit(contra_data(move_idx,5),contra_data(move_idx,6),1);
    % compute contrast for ipsi
    se_phase_contra_null(perm,1) = var([b_rest(1), b_prep(1), b_move(1)]);
    
    
    %% Simple effect: Area within Phase
    % shuffle Area labels across units for Instruct and Move independently
    prep_data = data_mat(data_mat(:,3)==2,:);
    prep_data(:,2) = prep_data(randperm(size(prep_data,1)),2);
    move_data = data_mat(data_mat(:,3)==3,:);
    move_data(:,2) = move_data(randperm(size(move_data,1)),2);
    % remove negative outliers (log scale can inflate these)
    rmv = prep_data(:,6)<-1;
    prep_data(rmv,:) = [];
    rmv = move_data(:,6)<-1;
    move_data(rmv,:) = [];
    % calculate slope of each new grouping for Instruct
    pmd_idx = find(prep_data(:,2)==0);
    m1_idx = find(prep_data(:,2)==1);
    b_pmd = polyfit(prep_data(pmd_idx,5),prep_data(pmd_idx,6),1);
    b_m1 = polyfit(prep_data(m1_idx,5),prep_data(m1_idx,6),1);
    % compute contrast for contra
    se_area_prep_null(perm,1) = abs(b_m1(1) - b_pmd(1));
    % calculate slope of each new grouping for Move
    pmd_idx = find(move_data(:,2)==0);
    m1_idx = find(move_data(:,2)==1);
    b_pmd = polyfit(move_data(pmd_idx,5),move_data(pmd_idx,6),1);
    b_m1 = polyfit(move_data(m1_idx,5),move_data(m1_idx,6),1);
    % compute contrast for contra
    se_area_move_null(perm,1) = abs(b_m1(1) - b_pmd(1));
    
    
    %% Simple effect: Pref within Phase
    % shuffle Area labels across units for Instruct and Move independently
    prep_data = data_mat(data_mat(:,3)==2,:);
    prep_data(:,4) = prep_data(randperm(size(prep_data,1)),4);
    move_data = data_mat(data_mat(:,3)==3,:);
    move_data(:,4) = move_data(randperm(size(move_data,1)),4);
    % remove negative outliers (log scale can inflate these)
    rmv = prep_data(:,6)<-1;
    prep_data(rmv,:) = [];
    rmv = move_data(:,6)<-1;
    move_data(rmv,:) = [];
    % calculate slope of each new grouping for Instruct
    ipsi_idx = find(prep_data(:,4)==0);
    contra_idx = find(prep_data(:,4)==1);
    b_ipsi = polyfit(prep_data(ipsi_idx,5),prep_data(ipsi_idx,6),1);
    b_contra = polyfit(prep_data(contra_idx,5),prep_data(contra_idx,6),1);
    % compute contrast for contra
    se_pref_prep_null(perm,1) = abs(b_contra(1) - b_ipsi(1));
    % calculate slope of each new grouping for Move
    ipsi_idx = find(move_data(:,4)==0);
    contra_idx = find(move_data(:,4)==1);
    b_ipsi = polyfit(move_data(ipsi_idx,5),move_data(ipsi_idx,6),1);
    b_contra = polyfit(move_data(contra_idx,5),move_data(contra_idx,6),1);
    % compute contrast for contra
    se_pref_move_null(perm,1) = abs(b_contra(1) - b_ipsi(1));
    
    
    %% permutation test for positive slope in Reach & Return data
    % shuffle unit labels across units for arm preference values
    shuff_ipsipref_m1_ap = ipsipref_m1_ap(randperm(length(ipsipref_m1_ap)));
    shuff_ipsipref_pmd_ap = ipsipref_pmd_ap(randperm(length(ipsipref_pmd_ap)));
    shuff_contrapref_m1_ap = contrapref_m1_ap(randperm(length(contrapref_m1_ap)));
    shuff_contrapref_pmd_ap = contrapref_pmd_ap(randperm(length(contrapref_pmd_ap)));
    % calculate slope of the shuffled dataset and log in null distribution
    b_ipsi_m1 = polyfit(shuff_ipsipref_m1_ap(ipsipref_m1_mod>=-1),...
        ipsipref_m1_mod(ipsipref_m1_mod>=-1), 1);
    b_ipsi_pmd = polyfit(shuff_ipsipref_pmd_ap(ipsipref_pmd_mod>=-1),...
        ipsipref_pmd_mod(ipsipref_pmd_mod>=-1), 1);
    b_contra_m1 = polyfit(shuff_contrapref_m1_ap(contrapref_m1_mod>=-1),...
        contrapref_m1_mod(contrapref_m1_mod>=-1), 1);
    b_contra_pmd = polyfit(shuff_contrapref_pmd_ap(contrapref_pmd_mod>=-1),...
        contrapref_pmd_mod(contrapref_pmd_mod>=-1), 1);
    extend_b_ipsi_m1_null(perm,1) = b_ipsi_m1(1);
    extend_b_ipsi_pmd_null(perm,1) = b_ipsi_pmd(1);
    extend_b_contra_m1_null(perm,1) = b_contra_m1(1);
    extend_b_contra_pmd_null(perm,1) = b_contra_pmd(1);
    
    
end


%% Main effect: Area p-value
% remove negative outliers (log scale can inflate these)
rmv = data_mat(:,6)<-1;
data_mat(rmv,:) = [];
% calculate slope of each observed grouping
pmd_idx = find(data_mat(:,2)==0);
m1_idx = find(data_mat(:,2)==1);
b_pmd = polyfit(data_mat(pmd_idx,5),data_mat(pmd_idx,6),1);
b_m1 = polyfit(data_mat(m1_idx,5),data_mat(m1_idx,6),1);
% compute observed contrast and p-value
me_area_obs = abs(b_m1(1) - b_pmd(1));
p_val.me.area = mean(me_area_null>me_area_obs);


%% Main effect: Phase p-value
% calculate slope of each observed grouping
rest_idx = find(data_mat(:,3)==1);
prep_idx = find(data_mat(:,3)==2);
move_idx = find(data_mat(:,3)==3);
b_rest = polyfit(data_mat(rest_idx,5), data_mat(rest_idx,6),1);
b_prep = polyfit(data_mat(prep_idx,5), data_mat(prep_idx,6),1);
b_move = polyfit(data_mat(move_idx,5), data_mat(move_idx,6),1);
% compute observed contrast and p-value
me_phase_obs = var([b_rest(1), b_prep(1), b_move(1)]);
p_val.me.phase = mean(me_phase_null>me_phase_obs);


%% Main effect: Pref
% calculate slope of observed grouping
ipsi_idx = find(data_mat(:,4)==0);
contra_idx = find(data_mat(:,4)==1);
b_ipsi = polyfit(data_mat(ipsi_idx,5),data_mat(ipsi_idx,6),1);
b_contra = polyfit(data_mat(contra_idx,5),data_mat(contra_idx,6),1);
% compute observed contrast and p-value
me_pref_obs = abs(b_contra(1) - b_ipsi(1));
p_val.me.pref = mean(me_pref_null>me_pref_obs);


%% Interaction effect: Area X Phase (all phases)
% calculate slope of each observed grouping
pmd_rest_idx = find(data_mat(:,2)==0 & data_mat(:,3)==1);
pmd_prep_idx = find(data_mat(:,2)==0 & data_mat(:,3)==2);
pmd_move_idx = find(data_mat(:,2)==0 & data_mat(:,3)==3);
m1_rest_idx = find(data_mat(:,2)==1 & data_mat(:,3)==1);
m1_prep_idx = find(data_mat(:,2)==1 & data_mat(:,3)==2);
m1_move_idx = find(data_mat(:,2)==1 & data_mat(:,3)==3);
b_pmd_rest = ...
    polyfit(data_mat(pmd_rest_idx,5),data_mat(pmd_rest_idx,6),1);
b_pmd_prep = ...
    polyfit(data_mat(pmd_prep_idx,5),data_mat(pmd_prep_idx,6),1);
b_pmd_move = ...
    polyfit(data_mat(pmd_move_idx,5),data_mat(pmd_move_idx,6),1);
b_m1_rest = ...
    polyfit(data_mat(m1_rest_idx,5),data_mat(m1_rest_idx,6),1);
b_m1_prep = ...
    polyfit(data_mat(m1_prep_idx,5),data_mat(m1_prep_idx,6),1);
b_m1_move = ...
    polyfit(data_mat(m1_move_idx,5),data_mat(m1_move_idx,6),1);
% compute observed contrast and p-value
ie_area_phase_obs = abs(...
    var([b_m1_rest(1),b_m1_prep(1),b_m1_move(1)]) - ...
    var([b_pmd_rest(1),b_pmd_prep(1),b_pmd_move(1)]) );
p_val.ie.area_phase = mean(ie_area_phase_null>ie_area_phase_obs);


%% Interaction effect: Area X Phase (Instruct, Move)
% calculate slope of each observed grouping
pmd_prep_idx = find(data_mat(:,2)==0 & data_mat(:,3)==2);
pmd_move_idx = find(data_mat(:,2)==0 & data_mat(:,3)==3);
m1_prep_idx = find(data_mat(:,2)==1 & data_mat(:,3)==2);
m1_move_idx = find(data_mat(:,2)==1 & data_mat(:,3)==3);
b_pmd_prep = ...
    polyfit(data_mat(pmd_prep_idx,5),data_mat(pmd_prep_idx,6),1);
b_pmd_move = ...
    polyfit(data_mat(pmd_move_idx,5),data_mat(pmd_move_idx,6),1);
b_m1_prep = ...
    polyfit(data_mat(m1_prep_idx,5),data_mat(m1_prep_idx,6),1);
b_m1_move = ...
    polyfit(data_mat(m1_move_idx,5),data_mat(m1_move_idx,6),1);
% compute observed contrast and p-value
ie_area_phaseIM_obs = abs(...
    abs(b_pmd_prep(1)-b_m1_prep(1)) - ...
    abs(b_pmd_move(1)-b_m1_move(1)) );
p_val.ie.area_phaseIM = mean(ie_area_phaseIM_null>ie_area_phaseIM_obs);



%% Interaction effect: Area X Pref
% calculate slope of each observed grouping
pmd_ipsi_idx = find(data_mat(:,2)==0 & data_mat(:,4)==0);
pmd_contra_idx = find(data_mat(:,2)==0 & data_mat(:,4)==1);
m1_ipsi_idx = find(data_mat(:,2)==1 & data_mat(:,4)==0);
m1_contra_idx = find(data_mat(:,2)==1 & data_mat(:,4)==1);
b_pmd_ipsi = ...
    polyfit(data_mat(pmd_ipsi_idx,5),data_mat(pmd_ipsi_idx,6),1);
b_pmd_contra = ...
    polyfit(data_mat(pmd_contra_idx,5),data_mat(pmd_contra_idx,6),1);
b_m1_ipsi = ...
    polyfit(data_mat(m1_ipsi_idx,5),data_mat(m1_ipsi_idx,6),1);
b_m1_contra = ...
    polyfit(data_mat(m1_contra_idx,5),data_mat(m1_contra_idx,6),1);
% compute observed contrast and p-value
ie_area_pref_obs = abs(...
    (b_m1_contra(1) - b_m1_ipsi(1)) - ...
    (b_pmd_contra(1) - b_pmd_ipsi(1)) );
p_val.ie.area_pref = mean(ie_area_pref_null>ie_area_pref_obs);


%% Interaction effect: Phase X Pref
% calculate slope of each observed grouping
rest_ipsi_idx = find(data_mat(:,3)==1 & data_mat(:,4)==0);
prep_ipsi_idx = find(data_mat(:,3)==2 & data_mat(:,4)==0);
move_ipsi_idx = find(data_mat(:,3)==3 & data_mat(:,4)==0);
rest_contra_idx = find(data_mat(:,3)==1 & data_mat(:,4)==1);
prep_contra_idx = find(data_mat(:,3)==2 & data_mat(:,4)==1);
move_contra_idx = find(data_mat(:,3)==3 & data_mat(:,4)==1);
b_rest_ipsi = ...
    polyfit(data_mat(rest_ipsi_idx,5),data_mat(rest_ipsi_idx,6),1);
b_prep_ipsi = ...
    polyfit(data_mat(prep_ipsi_idx,5),data_mat(prep_ipsi_idx,6),1);
b_move_ipsi = ...
    polyfit(data_mat(move_ipsi_idx,5),data_mat(move_ipsi_idx,6),1);
b_rest_contra = ...
    polyfit(data_mat(rest_contra_idx,5),data_mat(rest_contra_idx,6),1);
b_prep_contra = ...
    polyfit(data_mat(prep_contra_idx,5),data_mat(prep_contra_idx,6),1);
b_move_contra = ...
    polyfit(data_mat(move_contra_idx,5),data_mat(move_contra_idx,6),1);
% compute observed contrast and p-value
ie_phase_pref_obs = abs(...
    var([b_rest_contra(1), b_prep_contra(1), b_move_contra(1)]) - ...
    var([b_rest_ipsi(1), b_prep_ipsi(1), b_move_ipsi(1)]) );
p_val.ie.phase_pref = mean(ie_phase_pref_null>ie_phase_pref_obs);


%% Simple effect: Phase within Pref
% separate Ipsi and Contra levels of Pref
ipsi_data = data_mat(data_mat(:,4)==0,:);
contra_data = data_mat(data_mat(:,4)==1,:);
% remove negative outliers (log scale can inflate these)
rmv = ipsi_data(:,6)<-1;
ipsi_data(rmv,:) = [];
rmv = contra_data(:,6)<-1;
contra_data(rmv,:) = [];
% calculate slope of each observed grouping for Ipsi
rest_idx = find(ipsi_data(:,3)==1);
prep_idx = find(ipsi_data(:,3)==2);
move_idx = find(ipsi_data(:,3)==3);
b_rest = polyfit(ipsi_data(rest_idx,5),ipsi_data(rest_idx,6),1);
b_prep = polyfit(ipsi_data(prep_idx,5),ipsi_data(prep_idx,6),1);
b_move = polyfit(ipsi_data(move_idx,5),ipsi_data(move_idx,6),1);
% compute observed contrast for Contra and p-value
se_phase_ipsi_obs = var([b_rest(1), b_prep(1), b_move(1)]);
p_val.se.phase.ipsi = mean(se_phase_ipsi_null>se_phase_ipsi_obs);
% calculate slope of each observed grouping for Contra
rest_idx = find(contra_data(:,3)==1);
prep_idx = find(contra_data(:,3)==2);
move_idx = find(contra_data(:,3)==3);
b_rest = polyfit(contra_data(rest_idx,5),contra_data(rest_idx,6),1);
b_prep = polyfit(contra_data(prep_idx,5),contra_data(prep_idx,6),1);
b_move = polyfit(contra_data(move_idx,5),contra_data(move_idx,6),1);
% compute observed contrast for ipsi
se_phase_contra_obs = var([b_rest(1), b_prep(1), b_move(1)]);
p_val.se.phase.contra = mean(se_phase_contra_null>se_phase_contra_obs);


%% Simple effect: Area within Phase
% Separate Instruct and Move levels of Phase
prep_data = data_mat(data_mat(:,3)==2,:);
move_data = data_mat(data_mat(:,3)==3,:);
% remove negative outliers (log scale can inflate these)
rmv = prep_data(:,6)<-1;
prep_data(rmv,:) = [];
rmv = move_data(:,6)<-1;
move_data(rmv,:) = [];
% calculate slope of each observed grouping for Instruct
pmd_idx = find(prep_data(:,2)==0);
m1_idx = find(prep_data(:,2)==1);
b_pmd = polyfit(prep_data(pmd_idx,5),prep_data(pmd_idx,6),1);
b_m1 = polyfit(prep_data(m1_idx,5),prep_data(m1_idx,6),1);
% compute contrast for Instruct
se_area_prep_obs = abs(b_m1(1) - b_pmd(1));
p_val.se.area.prep = mean(se_area_prep_null>se_area_prep_obs);
% calculate slope of each observed grouping for Move
pmd_idx = find(move_data(:,2)==0);
m1_idx = find(move_data(:,2)==1);
b_pmd = polyfit(move_data(pmd_idx,5),move_data(pmd_idx,6),1);
b_m1 = polyfit(move_data(m1_idx,5),move_data(m1_idx,6),1);
% compute contrast for Move
se_area_move_obs = abs(b_m1(1) - b_pmd(1));
p_val.se.area.move = mean(se_area_move_null>se_area_move_obs);


%% Simple effect: Pref within Phase
% Separate Instruct and Move levels of Phase
prep_data = data_mat(data_mat(:,3)==2,:);
move_data = data_mat(data_mat(:,3)==3,:);
% remove negative outliers (log scale can inflate these)
rmv = prep_data(:,6)<-1;
prep_data(rmv,:) = [];
rmv = move_data(:,6)<-1;
move_data(rmv,:) = [];
% calculate slope of each observed grouping for Instruct
ipsi_idx = find(prep_data(:,4)==0);
contra_idx = find(prep_data(:,4)==1);
b_ipsi = polyfit(prep_data(ipsi_idx,5),prep_data(ipsi_idx,6),1);
b_contra = polyfit(prep_data(contra_idx,5),prep_data(contra_idx,6),1);
% compute contrast for Instruct
se_pref_prep_obs = abs(b_contra(1) - b_ipsi(1));
p_val.se.pref.prep = mean(se_pref_prep_null>se_pref_prep_obs);
% calculate slope of each observed grouping for Move
ipsi_idx = find(move_data(:,4)==0);
contra_idx = find(move_data(:,4)==1);
b_ipsi = polyfit(move_data(ipsi_idx,5),move_data(ipsi_idx,6),1);
b_contra = polyfit(move_data(contra_idx,5),move_data(contra_idx,6),1);
% compute contrast for Move
se_pref_move_obs = abs(b_contra(1) - b_ipsi(1));
p_val.se.pref.move = mean(se_pref_move_null>se_pref_move_obs);


%% permutation test for positive slope in Reach & Return data
b_ipsi_m1 = polyfit(ipsipref_m1_ap(ipsipref_m1_mod>=-1),...
    ipsipref_m1_mod(ipsipref_m1_mod>=-1), 1);
b_ipsi_pmd = polyfit(ipsipref_pmd_ap(ipsipref_pmd_mod>=-1),...
    ipsipref_pmd_mod(ipsipref_pmd_mod>=-1), 1);
b_contra_m1 = polyfit(contrapref_m1_ap(contrapref_m1_mod>=-1),...
    contrapref_m1_mod(contrapref_m1_mod>=-1), 1);
b_contra_pmd = polyfit(contrapref_pmd_ap(contrapref_pmd_mod>=-1),...
    contrapref_pmd_mod(contrapref_pmd_mod>=-1), 1);
extend_b_ipsi_m1_obs = b_ipsi_m1(1);
extend_b_ipsi_pmd_obs = b_ipsi_pmd(1);
extend_b_contra_m1_obs = b_contra_m1(1);
extend_b_contra_pmd_obs = b_contra_pmd(1);
% log p-value as proportion of null distribution with slope at least as
% extreme as the observed value (one-sided)
p_val.se.unit_label.ipsi_pref.m1 = ...
    mean(abs(extend_b_ipsi_m1_null)>extend_b_ipsi_m1_obs);
p_val.se.unit_label.ipsi_pref.pmd = ...
    mean(abs(extend_b_ipsi_pmd_null)>extend_b_ipsi_pmd_obs);
p_val.se.unit_label.contra_pref.m1 = ...
    mean(abs(extend_b_contra_m1_null)>extend_b_contra_m1_obs);
p_val.se.unit_label.contra_pref.pmd = ...
    mean(abs(extend_b_contra_pmd_null)>extend_b_contra_pmd_obs);
    


    
    