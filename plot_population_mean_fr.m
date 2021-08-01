function [mean_mod, mean_ap] = plot_population_mean_fr(unit_data)

%
% Plots the mean firing rate for PMd and M1 populations separately for each 
% hand (ipsi/contra), locked to onset of instructional cue and movement 
% onset
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
% mean_mod - mean modulation strength across the population
%            (struct, fields: pmd, m1)
%                (nested struct, fields: rest, prep, move)
%                    each field contains a 1x2 vector with the ipsi 
%                    modulation in the first index and the contra  
%                    modulation in the second
%
% mean_ap - mean arm preference across the population
%           (struct, fields: pmd, m1)
%               (nested struct, fields: rest, prep, move)
%                   each field the scalar arm preference value (positive
%                   indicates contra preference)
% 


%% log the firing rates for each population

[pmd_fr_ipsi, pmd_fr_contra, pmd_mod, pmd_ap] =...
    collect_pop_fr(unit_data, 'pmd');
[m1_fr_ipsi, m1_fr_contra, m1_mod, m1_ap] =...
    collect_pop_fr(unit_data, 'm1');


%% compute the modulation and arm preference means

mean_mod.pmd = squeeze(nanmean(pmd_mod,1));
mean_mod.m1 = squeeze(nanmean(m1_mod,1));
mean_ap.pmd = squeeze(nanmean(pmd_ap,1));
mean_ap.m1 = squeeze(nanmean(m1_ap,1));
    

%% Plot

figure('Position', [100, 100, 1500, 500])

% Contra responses
subplot(1,2,1)
pbaspect([3,1,1])
sem_plot(pmd_fr_contra, 'g')
sem_plot(m1_fr_contra, 'k')
ylabel('Firing rate (hz)')
title('Contralateral')

% Ipsi responses
subplot(1,2,2)
pbaspect([3,1,1])
sem_plot(pmd_fr_ipsi, 'g')
sem_plot(m1_fr_ipsi, 'k')
title('Ipsilateral')
f = gcf;
line_objects = f.Children(2).Children;
legend(line_objects([7,3]), {'PMd', 'M1'})

end


%% helper functions

function [pop_fr_ipsi, pop_fr_contra, pop_mod, pop_ap] =...
    collect_pop_fr(unit_data, area)

% isolate the desired area and filter out only non-repeated SU's
area = strcmp(area, 'm1');
unit_data = unit_data(([unit_data.unit_type]==1) & ...
                      ([unit_data.repeat]==0) & ...
                      ([unit_data.area]==area));
                  
% collect the mean firing rate for all units in the desired area
pop_fr_ipsi = zeros(length(unit_data),78);
pop_fr_contra = zeros(length(unit_data),78);
pop_mod = zeros(length(unit_data),3,2);
pop_ap = zeros(length(unit_data),3);
for unit = 1:length(unit_data)
    [fr_ipsi, fr_contra, mod, ap] = compute_unit_mean_fr(unit_data(unit));
    pop_fr_ipsi(unit,:) = fr_ipsi;
    pop_fr_contra(unit,:) = fr_contra;
    pop_mod(unit,:,:) = mod;
    pop_ap(unit,:) = ap;
end

% % remove NaN rows
% rmv = isnan(pop_fr_ipsi(:,1));
% pop_fr_ipsi(rmv,:) = [];
% pop_fr_contra(rmv,:) = [];
% pop_mod(rmv,:,:) = [];
% pop_ap(rmv,:) = [];

end


function [fr_ipsi, fr_contra, mod, ap] =...
    compute_unit_mean_fr(unit_data_row)

% only include units that were significantly modulated for at least one of
% the two arms in at least one of the two phases. return NaN array if this
% is not the case
% if ~(unit_data_row.p_mod.ipsi.prep<0.05 || ...
%      unit_data_row.p_mod.ipsi.move<0.05 || ...
%      unit_data_row.p_mod.contra.prep<0.05 || ...
%      unit_data_row.p_mod.contra.move<0.05)
%  
%     fr_ipsi = nan(1,78);
%     fr_contra = nan(1,78);
%     mod = nan(3,2);
%     ap = nan(3,1);
%     return
% end

% log modulation and arm pref, then determine which config to use for 
% firing rate data
if unit_data_row.hem==0 % left hemisphere
    mod = [unit_data_row.modulation.rest;...
           unit_data_row.modulation.prep;...
           unit_data_row.modulation.move];   % [ipsi, contra]
    ap = [unit_data_row.arm_pref.rest;...
          unit_data_row.arm_pref.prep;...
          unit_data_row.arm_pref.move];
    ipsi_config = 3;
    contra_config = 1;
else % right hemisphere
    mod = fliplr([unit_data_row.modulation.rest;...
                  unit_data_row.modulation.prep;...
                  unit_data_row.modulation.move]);   % [ipsi, contra]
    ap = -[unit_data_row.arm_pref.rest;...
           unit_data_row.arm_pref.prep;...
           unit_data_row.arm_pref.move];
    ipsi_config = 1;
    contra_config = 3;
end

% log firing rates for each trial, separated by limb
fr_ipsi = cell(6,1);
fr_contra = cell(6,1);
for target = 1:6
    
    fr_ipsi{target} = ...
        [unit_data_row.ipsi.config(ipsi_config).target(target).rest,...
         unit_data_row.ipsi.config(ipsi_config).target(target).prep,...
         unit_data_row.ipsi.config(ipsi_config).target(target).move];
    fr_ipsi{target} = fr_ipsi{target}(:,[11:52, 53:88]);
    
    fr_contra{target} = ...
        [unit_data_row.contra.config(contra_config).target(target).rest,...
         unit_data_row.contra.config(contra_config).target(target).prep,...
         unit_data_row.contra.config(contra_config).target(target).move];
    fr_contra{target} = fr_contra{target}(:,[11:52, 53:88]);
    
end

fr_ipsi = mean(vertcat(fr_ipsi{:}),1);
fr_contra = mean(vertcat(fr_contra{:}),1);

end


function sem_plot(y, c)

shadedErrorBar(-300:20:500, mean(y(:,1:41),1), ...
    std(y(:,1:41),1)/sqrt(size(y,1)),...
    'lineprops',c)
shadedErrorBar(600:20:1300, mean(y(:,43:78),1), ...
    std(y(:,43:78),1)/sqrt(size(y,1)),...
    'lineprops',c)
xline(0,'r')
xline(800,'g')
xlim([-300,1300])
xticks([-300,0,300,800,1100])
xticklabels({'-300ms','Instruct','+300ms','Move','+300ms'})

end




