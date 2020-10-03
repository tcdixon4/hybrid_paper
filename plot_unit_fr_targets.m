function [] = plot_unit_fr_targets(unit_data, unit)

%
% Plots trial-averaged firing rate traces for single-units, locked to onset
% of instructional cue and movement onset
%
%
% INPUTS: 
%
% unit_data - unit-separated data struct containing hemisphere, brain area,
%             firing rate and other metrics
%             (struct: 1 x num_units)
%
% unit - selected unit to plot, referenced as the index within unit_data
%
% OUTPUTS:
%
% Only plots
% 


%% Select the correct configs to use so all targets are reached

if unit_data(unit).hem==0 % left hemisphere
    ipsi_config = 3;
    contra_config = 1;
else % right hemisphere
    ipsi_config = 1;
    contra_config = 3;
end


%% Log average modulation across time for all units

% log firing rates for each trial, separated by limb
fr_ipsi = cell(6,1);
fr_contra = cell(6,1);
for target = 1:6
    
    fr_ipsi{target} = ...
        [unit_data(unit).ipsi.config(ipsi_config).target(target).rest,...
        unit_data(unit).ipsi.config(ipsi_config).target(target).prep,...
        unit_data(unit).ipsi.config(ipsi_config).target(target).move];
    fr_ipsi{target} = fr_ipsi{target}(:,[11:52, 53:88]);
    
    fr_contra{target} = ...
        [unit_data(unit).contra.config(contra_config).target(target).rest,...
        unit_data(unit).contra.config(contra_config).target(target).prep,...
        unit_data(unit).contra.config(contra_config).target(target).move];
    fr_contra{target} = fr_contra{target}(:,[11:52, 53:88]);
    
end
    

%% Plot all targets together, but separately for ipsi and contra

figure('Position', [100, 100, 1500, 500])
subplot(1,2,1)
pbaspect([3,1,1])
for target = 1:6
    switch target
        case 1
            c = 'y';
        case 2
            c = 'm';
        case 3
            c = 'c';
        case 4
            c = 'r';
        case 5
            c = 'g';
        case 6
            c = 'b';
    end
    y = fr_contra{target};
    shadedErrorBar(-300:20:500, mean(y(:,1:41),1), ...
        std(y(:,1:41),1)/sqrt(size(y,1)),...
        'lineprops',c)
    shadedErrorBar(600:20:1300, mean(y(:,43:78),1), ...
        std(y(:,43:78),1)/sqrt(size(y,1)),...
        'lineprops',c)
    xticks([-300,0,300,800,1100])
    xticklabels({'-300ms','Instruct','+300ms','Move','+300ms'})
end
hold on
plot([0,0],[0,50],'r')
plot([800,800],[0,50],'g')
xlim([-300,1300])
ylim([0,50])
ylabel('Firing rate (hz)')
title('Contralateral')

subplot(1,2,2)
pbaspect([3,1,1])
for target = 1:6
    switch target
        case 1
            c = 'y';
        case 2
            c = 'm';
        case 3
            c = 'c';
        case 4
            c = 'r';
        case 5
            c = 'g';
        case 6
            c = 'b';
    end
    y = fr_ipsi{target};
    shadedErrorBar(-300:20:500, mean(y(:,1:41),1), ...
        std(y(:,1:41),1)/sqrt(size(y,1)),...
        'lineprops',c)
    shadedErrorBar(600:20:1300, mean(y(:,43:78),1), ...
        std(y(:,43:78),1)/sqrt(size(y,1)),...
        'lineprops',c)
    xticks([-300,0,300,800,1100])
    xticklabels({'-300ms','Instruct','+300ms','Move','+300ms'})
end
hold on
plot([0,0],[0,50],'r')
plot([800,800],[0,50],'g')
xlim([-300,1300])
ylim([0,50])
title('Ipsilateral')
legend({'Target 1','','Target 2','','Target 3','',...
    'Target 4','','Target 5','','Target 6','',...
    'Instruct','Move'})



