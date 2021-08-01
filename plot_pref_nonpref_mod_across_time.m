function [] = plot_pref_nonpref_mod_across_time(unit_data)

%
% Plots mean modulation strength across time, separated by preferred arm 
% and non-preferred arm responses.
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
% plots
%


%% Log mean and sd across sessions and isolate intended time windows

modulation = pref_nonpref_mod_across_time(unit_data);

n = size(modulation.pref,1);
mu_pref = nanmean(modulation.pref(:, [11:52,53:88]));
sem_pref = nanstd(modulation.pref(:, [11:52,53:88]))/sqrt(n);
mu_nonpref = nanmean(modulation.nonpref(:, [11:52,53:88]));
sem_nonpref = nanstd(modulation.nonpref(:, [11:52,53:88]))/sqrt(n);


%% Plot

figure('Position', [100, 100, 650, 500])
pbaspect([3,1,1])

shadedErrorBar(-300:20:500, mu_pref(1:41), sem_pref(1:41), ...
    'lineprops','b')
shadedErrorBar(600:20:1300, mu_pref(43:78), sem_pref(43:78), ...
    'lineprops','b')

shadedErrorBar(-300:20:500, mu_nonpref(1:41), sem_nonpref(1:41), ...
    'lineprops','r')
shadedErrorBar(600:20:1300, mu_nonpref(43:78), sem_nonpref(43:78), ...
    'lineprops','r')

hold on
xline(0,'k') % Instruct
xline(800,'k') % Move
xlim([-300,1300])
xticks([-300,0,500,600,800,1300])
xticklabels({'-300ms','Instruct','+500ms','','Move','+500ms'})
ylim([0,8])
yticks(linspace(0,8,5))
ylabel({'Modulation', '(a.u.)'})
ax = gca;
ax.Clipping = 'off';

legend({'Preferred arm','','Non-preferred arm','',...
    'Instruct','Move'})
