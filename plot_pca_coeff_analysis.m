function [] = plot_pca_coeff_analysis(...
    explained_total, weighted_arm_pref_total)

%
% INPUTS: 
%
% explained_total - cell array, with each cell containing data for a single
%                  session indexed as: 
%                 explained_total{[session]}.[phase].[hand].[trained space]
%                  which contains the proportion of total variance captured 
%                  by each component in a 1xp vector
%
% weighted_arm_pref_total - cell array, with each cell containing data for
%                  a single session indexed as: 
%                 explained_total{[session]}.[phase].[hand].[trained space]
%                  which contains the coefficient-weighted arm preference 
%                  for each component in a 1xp vector
% OUTPUTS:
%
% plots
%

%% Set up helper variables and pre-allocate

num_sessions = length(explained_total);
var_ratio_rest = zeros(num_sessions,5);
var_ratio_prep = zeros(num_sessions,5);
var_ratio_move = zeros(num_sessions,5);
coeff_rest = zeros(num_sessions,5);
coeff_prep = zeros(num_sessions,5);
coeff_move = zeros(num_sessions,5);


%% Log the variance ratio and coefficient values for the top 5 PC's

for k = num_sessions:-1:1
    % var_ratio convention: right over left
    % first, models trained on left hand trials
    var_ratio_rest(k,:) = ...
        log10(explained_total{k}.rest.r_hand.l_space(1:5)./...
        explained_total{k}.rest.l_hand.l_space(1:5));
    var_ratio_prep(k,:) = ...
        log10(explained_total{k}.prep.r_hand.l_space(1:5)./...
        explained_total{k}.prep.l_hand.l_space(1:5));
    var_ratio_move(k,:) = ...
        log10(explained_total{k}.move.r_hand.l_space(1:5)./...
        explained_total{k}.move.l_hand.l_space(1:5));
    coeff_rest(k,:) = weighted_arm_pref_total{k}.rest.l_hand(1:5);
    coeff_prep(k,:) = weighted_arm_pref_total{k}.prep.l_hand(1:5);
    coeff_move(k,:) = weighted_arm_pref_total{k}.move.l_hand(1:5);
    
    % next, models trained on right hand trials
    var_ratio_rest(k+num_sessions,:) = ...
        log10(explained_total{k}.rest.r_hand.r_space(1:5)./...
        explained_total{k}.rest.l_hand.r_space(1:5));
    var_ratio_prep(k+num_sessions,:) = ...
        log10(explained_total{k}.prep.r_hand.r_space(1:5)./...
        explained_total{k}.prep.l_hand.r_space(1:5));
    var_ratio_move(k+num_sessions,:) = ...
        log10(explained_total{k}.move.r_hand.r_space(1:5)./...
        explained_total{k}.move.l_hand.r_space(1:5));
    coeff_rest(k+num_sessions,:) = ...
        weighted_arm_pref_total{k}.rest.r_hand(1:5);
    coeff_prep(k+num_sessions,:) = ...
        weighted_arm_pref_total{k}.prep.r_hand(1:5);
    coeff_move(k+num_sessions,:) = ...
        weighted_arm_pref_total{k}.move.r_hand(1:5);
end


%% Plot scatter of var_ratio vs Coeff-weighted arm pref

figure('Name','PCA coefficient analysis', 'Position',[100,100,1800,400])
% Rest phase
subplot(1,3,1)
for k = 5:-1:1
    scatter(coeff_rest(1:num_sessions,k), var_ratio_rest(1:num_sessions,k))
    hold on
    scatter(coeff_rest(num_sessions+1:end,k), ...
        var_ratio_rest(num_sessions+1:end,k), 'filled')
end
grid on
xlim([-0.8,0.8])
xticks(-0.8:0.4:0.8)
xticklabels({'-0.8','','0','','0.8'})
ylim([-2,2])
yticks(-2:1:2)
yticklabels({'10^{-2}','','1','','10^2'})
[r,~] = corr(coeff_rest(:), var_ratio_rest(:), 'rows','complete');
annotation('textbox',[.3 .25 .1 .1],'String',strcat('r = ',num2str(r)),...
    'FitBoxToText','on');
title('Rest')
pbaspect([1,0.9,1])
xlabel({'Coefficient weighted', 'arm preference'})
ylabel({'Right/Left variance', 'ratio (log scale)'})
ax = gca;
ax.Clipping = 'off';

% Instruct phase
subplot(1,3,2)
for k = 5:-1:1
    scatter(coeff_prep(1:num_sessions,k), var_ratio_prep(1:num_sessions,k))
    hold on
    scatter(coeff_prep(num_sessions+1:end,k), ...
        var_ratio_prep(num_sessions+1:end,k), 'filled')
end
grid on
xlim([-0.8,0.8])
xticks(-0.8:0.4:0.8)
xticklabels({'-0.8','','0','','0.8'})
ylim([-2,2])
yticks(-2:1:2)
yticklabels({'10^{-2}','','1','','10^2'})
[r,~] = corr(coeff_prep(:), var_ratio_prep(:), 'rows','complete');
annotation('textbox',[.6 .25 .1 .1],'String',strcat('r = ',num2str(r)),...
    'FitBoxToText','on');
title('Prep')
pbaspect([1,0.9,1])
ax = gca;
ax.Clipping = 'off';

% Move phase
subplot(1,3,3)
for k = 5:-1:1
    scatter(coeff_move(1:num_sessions,k), var_ratio_move(1:num_sessions,k))
    hold on
    scatter(coeff_move(num_sessions+1:end,k), ...
        var_ratio_move(num_sessions+1:end,k), 'filled')
end
grid on
xlim([-0.8,0.8])
xticks(-0.8:0.4:0.8)
xticklabels({'-0.8','','0','','0.8'})
ylim([-2,2])
yticks(-2:1:2)
yticklabels({'10^{-2}','','1','','10^2'})
[r,~] = corr(coeff_move(:), var_ratio_move(:), 'rows','complete');
annotation('textbox',[.9 .25 .1 .1],'String',strcat('r = ',num2str(r)),...
    'FitBoxToText','on');
title('Move')
pbaspect([1,0.9,1])
ax = gca;
ax.Clipping = 'off';
legend
% Fix colors outside matlab to show which arm model each component belongs
% to (right:purple or left:yellow) but leave this way for now in case there
% is a desire to see which number component each datapoint belongs to






