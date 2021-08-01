function [p] = ap_vs_depth(unit_data)

%
% Blah
%
%
% INPUTS: 
%
% unit_data - struct: 1 x num_units
%             unit-separated data struct containing hemisphere, brain area,
%             firing rate and other metrics
%
% OUTPUTS:
%
% p - significance testing results
%
% figures
%
%


%% log vectors for depth, AP, modulation, area, and hemisphere

% isolate only single units and non-repeated units
unit_data = unit_data(([unit_data.unit_type]==1) & ([unit_data.repeat]==false));
num_units = length(unit_data);

% some of the info is burried in nested structures, so just use for-loop
% even though it's less efficient for some things
for unit = num_units:-1:1
    
    depth(unit) = unit_data(unit).depth; % Shallowest:1, Deepest:32
    is_m1(unit) = unit_data(unit).area; % PMd:0, M1:1
    is_rhem(unit) = unit_data(unit).hem; % Left:0, Right:1
    probe(unit) = unit_data(unit).probe_id; % unique ID for each probe insertion
    rest_fr(unit) = unit_data(unit).mu_rest; % resting firing rate
    
    if is_rhem(unit) % Right hem
        % positive ap means contra pref: [Prep, Move]
        ap(unit,1:2) = -[unit_data(unit).arm_pref.prep,...
                         unit_data(unit).arm_pref.move];
        % mod, log scale: [Ipsi_Prep, Contra_Prep, Ipsi_Move, Contra_Move]
        ms(unit,1:4) = log10([fliplr(unit_data(unit).modulation.prep),...
                               fliplr(unit_data(unit).modulation.move)]); 
    else % Left hem
        % positive ap means contra pref: [Prep, Move]
        ap(unit,1:2) = [unit_data(unit).arm_pref.prep,...
                        unit_data(unit).arm_pref.move];
        % mod, log scale: [Ipsi_Prep, Contra_Prep, Ipsi_Move, Contra_Move]
        ms(unit,1:4) = log10([unit_data(unit).modulation.prep,...
                               unit_data(unit).modulation.move]);
    end
    
end

%% Index the same arm preference regimes from previous analysis
contra_idx = ap>=0.4;
ipsi_idx = ap<=-0.4;
neut_idx = (ap>-0.4) & (ap<0.4);


%% plot distribution of unit depths
figure('Name','Unit depth sampling distribution')

subplot(2,1,1)
histogram(depth(~is_m1), 0.5:1:32.5)
xticks(1:10:31)
xticklabels({'0','1','2','3'})
title('PMd')

subplot(2,1,2)
histogram(depth(is_m1), 0.5:1:32.5)
xticks(1:10:31)
xticklabels({'0','1','2','3'})
xlabel('Distance from most superficial electrode (mm)')
ylabel('Unit count')
title('M1')


%% plot scatter of arm pref vs depth
figure('Name','Arm pref vs recording depth scatter')

subplot(2,2,1)
plot(ap(~is_m1,1), -depth(~is_m1), 'o')
title('PMd - Instruct')
format_scatter_axis

subplot(2,2,2)
plot(ap(~is_m1,2), -depth(~is_m1), 'o')
title('PMd - move')
format_scatter_axis

subplot(2,2,3)
plot(ap(is_m1,1), -depth(is_m1), 'o')
title('M1 - Instruct')
format_scatter_axis
ylbl = ylabel('Distance from most superficial electrode (mm)');
ylbl.Position(1) = -1.4;
ylbl.Position(2) = 5;
xlabel('Arm preference')

subplot(2,2,4)
plot(ap(is_m1,2), -depth(is_m1), 'o')
title('M1 - Move')
format_scatter_axis


%% plot the unit depths from each arm preference regime
% figure('Name','Average unit depth in each arm pref regime')
% 
% subplot(2,2,1)
% plot_regime_depth_violin(...
%     ~is_m1, contra_idx(:,1), neut_idx(:,1), ipsi_idx(:,1), depth)
% title('PMd - Instruct')
% 
% subplot(2,2,2)
% plot_regime_depth_violin(...
%     ~is_m1, contra_idx(:,2), neut_idx(:,2), ipsi_idx(:,2), depth)
% title('PMd - Move')
% 
% subplot(2,2,3)
% plot_regime_depth_violin(...
%     is_m1, contra_idx(:,1), neut_idx(:,1), ipsi_idx(:,1), depth)
% title('M1 - Instruct')
% xlabel('Distance from most superficial electrode (mm)')
% ylabel('Unit count')
% 
% subplot(2,2,4)
% plot_regime_depth_violin(...
%     is_m1, contra_idx(:,2), neut_idx(:,2), ipsi_idx(:,2), depth)
% title('M1 - Move')


%% plot the unit depths from each arm preference regime as violin plot

figure('Name','Average unit depth in each arm pref regime')

subplot(2,2,1)
plot_regime_depth_violin(...
    ~is_m1, contra_idx(:,1), neut_idx(:,1), ipsi_idx(:,1), depth)
title('PMd - Instruct')

subplot(2,2,2)
plot_regime_depth_violin(...
    ~is_m1, contra_idx(:,2), neut_idx(:,2), ipsi_idx(:,2), depth)
title('PMd - Move')

subplot(2,2,3)
plot_regime_depth_violin(...
    is_m1, contra_idx(:,1), neut_idx(:,1), ipsi_idx(:,1), depth)
title('M1 - Instruct')
xlabel('Distance from most superficial electrode (mm)')
ylabel('Unit count')

subplot(2,2,4)
plot_regime_depth_violin(...
    is_m1, contra_idx(:,2), neut_idx(:,2), ipsi_idx(:,2), depth)
title('M1 - Move')


%% test significance of differences in mean unit depth between each regime

% first do a permutation-based one-way ANOVA (uses variance as the contrast
% metric, not the F statistic)

% PMd Instruct
p.pmd.instruct.main = perm1WayANOVA(...
    {depth((~is_m1) & (contra_idx(:,1)'))'...
     depth((~is_m1) & (neut_idx(:,1)'))',...
     depth((~is_m1) & (ipsi_idx(:,1)'))'}, 10000);
% PMd Move
p.pmd.move.main = perm1WayANOVA(...
    {depth((~is_m1) & (contra_idx(:,2)'))'...
     depth((~is_m1) & (neut_idx(:,2)'))',...
     depth((~is_m1) & (ipsi_idx(:,2)'))'}, 10000);

% M1 Instruct
p.m1.instruct.main = perm1WayANOVA(...
    {depth((is_m1) & (contra_idx(:,1)'))'...
     depth((is_m1) & (neut_idx(:,1)'))',...
     depth((is_m1) & (ipsi_idx(:,1)'))'}, 10000);
% M1 Move
p.m1.move.main = perm1WayANOVA(...
    {depth((is_m1) & (contra_idx(:,2)'))'...
     depth((is_m1) & (neut_idx(:,2)'))',...
     depth((is_m1) & (ipsi_idx(:,2)'))'}, 10000);


% second, perform all 3 pairwise comparisons using permutation test

% PMd Instruct
p.pmd.instruct.ic = ...
    two_sample_perm_test(depth((~is_m1) & (contra_idx(:,1)')),...
                         depth((~is_m1) & (ipsi_idx(:,1)')), 10000);
p.pmd.instruct.in = ...
    two_sample_perm_test(depth((~is_m1) & (ipsi_idx(:,1)')),...
                         depth((~is_m1) & (neut_idx(:,1)')), 10000);
p.pmd.instruct.cn = ...
    two_sample_perm_test(depth((~is_m1) & (contra_idx(:,1)')),...
                         depth((~is_m1) & (neut_idx(:,1)')), 10000);
% PMd Move
p.pmd.move.ic = ...
    two_sample_perm_test(depth((~is_m1) & (contra_idx(:,2)')),...
                         depth((~is_m1) & (ipsi_idx(:,2)')), 10000);
p.pmd.move.in = ...
    two_sample_perm_test(depth((~is_m1) & (ipsi_idx(:,2)')),...
                         depth((~is_m1) & (neut_idx(:,2)')), 10000);
p.pmd.move.cn = ...
    two_sample_perm_test(depth((~is_m1) & (contra_idx(:,2)')),...
                         depth((~is_m1) & (neut_idx(:,2)')), 10000);

% M1 Instruct
p.m1.instruct.ic = ...
    two_sample_perm_test(depth((is_m1) & (contra_idx(:,1)')),...
                         depth((is_m1) & (ipsi_idx(:,1)')), 10000);
p.m1.instruct.in = ...
    two_sample_perm_test(depth((is_m1) & (ipsi_idx(:,1)')),...
                         depth((is_m1) & (neut_idx(:,1)')), 10000);
p.m1.instruct.cn = ...
    two_sample_perm_test(depth((is_m1) & (contra_idx(:,1)')),...
                         depth((is_m1) & (neut_idx(:,1)')), 10000);
% M1 Move
p.m1.move.ic = ...
    two_sample_perm_test(depth((is_m1) & (contra_idx(:,2)')),...
                         depth((is_m1) & (ipsi_idx(:,2)')), 10000);
p.m1.move.in = ...
    two_sample_perm_test(depth((is_m1) & (ipsi_idx(:,2)')),...
                         depth((is_m1) & (neut_idx(:,2)')), 10000);
p.m1.move.cn = ...
    two_sample_perm_test(depth((is_m1) & (contra_idx(:,2)')),...
                         depth((is_m1) & (neut_idx(:,2)')), 10000);


end


%% Helper functions

function [] = format_scatter_axis()

axis([-1 1 -33 0])
xticks(-1:1:1)
xticklabels({'-1','','1'})
xlim([-1,1])
yticks(-31:10:1)
yticklabels({'3','2','1','0'})
ylim([-31,-1])
box on
grid on
ax = gca;
ax.Clipping = 'off';

end


function [] = plot_regime_depth_violin(...
    area_idx, contra_idx, neut_idx, ipsi_idx, depth)

data = {-depth((area_idx) & (contra_idx')),...
        -depth((area_idx) & (neut_idx')),...
        -depth((area_idx) & (ipsi_idx')),...
        -depth(area_idx)};

violinplot(cell2table(data), [], 'ShowData',false, 'ShowMean',true);
xticks(1:4)
xticklabels({'Contra','Neutral','Ipsi','Total'})
xtickangle(45)
xlim([0.5,4.5])
yticks(-31:10:-1)
yticklabels({'3','2','1','0'})
ylim([-32,-1])
grid on
box on

end


