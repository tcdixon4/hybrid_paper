function [] = plot_arm_pref(unit_data)

%
% Makes several plots visualizing modulation strength and arm preference
% distributions and correlations.
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
% Plots:
% 1 - Emprical CDF's for modulation strength for each brain area x hand x
%     task phase
%
% 2:9 - Arm preference distributions visualized in different forms
%
% 10:12 - Modulation strength vs arm preference visualizations
%
% 13:22 - Cumulative modulation accounted for across arm preference
%         spectrum
%


%% Prepare arm preference and modulation data structures

% isolate only single units
unit_data = unit_data([unit_data.unit_type]==1);

m1 = unit_data([unit_data.area]==1);
pmd = unit_data([unit_data.area]==0);

clearvars unit_data

% arm_pref_ic negative means ipsi pref, positive means contra pref
% arm_pref_lr negative means left pref, positive means right pref
% first column of modulation_ic is ipsi, second is contra
% first column of modulation_lr is left, second is right

for unit = length(m1):-1:1
    modulation_m1_lr.rest(unit,:) = m1(unit).modulation.rest;
    modulation_m1_lr.prep(unit,:) = m1(unit).modulation.prep;
    modulation_m1_lr.move(unit,:) = m1(unit).modulation.move;
    
    arm_pref_m1_lr.rest(unit,:) = m1(unit).arm_pref.rest;
    arm_pref_m1_lr.prep(unit,:) = m1(unit).arm_pref.prep;
    arm_pref_m1_lr.move(unit,:) = m1(unit).arm_pref.move;
    
    p_mod_m1.ipsi.prep(unit) = m1(unit).p_mod.ipsi.prep;
    p_mod_m1.contra.prep(unit) = m1(unit).p_mod.contra.prep;
    p_mod_m1.ipsi.move(unit) = m1(unit).p_mod.ipsi.move;
    p_mod_m1.contra.move(unit) = m1(unit).p_mod.contra.move;
    
    if m1(unit).hem == 0
        modulation_m1_ic.rest(unit,:) = m1(unit).modulation.rest;
        modulation_m1_ic.prep(unit,:) = m1(unit).modulation.prep;
        modulation_m1_ic.move(unit,:) = m1(unit).modulation.move;
        
        arm_pref_m1_ic.rest(unit,:) = m1(unit).arm_pref.rest;
        arm_pref_m1_ic.prep(unit,:) = m1(unit).arm_pref.prep;
        arm_pref_m1_ic.move(unit,:) = m1(unit).arm_pref.move;
    elseif m1(unit).hem == 1
        arm_pref_m1_ic.rest(unit,:) = -m1(unit).arm_pref.rest;
        arm_pref_m1_ic.prep(unit,:) = -m1(unit).arm_pref.prep;
        arm_pref_m1_ic.move(unit,:) = -m1(unit).arm_pref.move;
        
        modulation_m1_ic.rest(unit,:) = fliplr(m1(unit).modulation.rest);
        modulation_m1_ic.prep(unit,:) = fliplr(m1(unit).modulation.prep);
        modulation_m1_ic.move(unit,:) = fliplr(m1(unit).modulation.move);
    end
end

for unit = length(pmd):-1:1
    modulation_pmd_lr.rest(unit,:) = pmd(unit).modulation.rest;
    modulation_pmd_lr.prep(unit,:) = pmd(unit).modulation.prep;
    modulation_pmd_lr.move(unit,:) = pmd(unit).modulation.move;
    
    arm_pref_pmd_lr.rest(unit,:) = pmd(unit).arm_pref.rest;
    arm_pref_pmd_lr.prep(unit,:) = pmd(unit).arm_pref.prep;
    arm_pref_pmd_lr.move(unit,:) = pmd(unit).arm_pref.move;
    
    p_mod_pmd.ipsi.prep(unit) = pmd(unit).p_mod.ipsi.prep;
    p_mod_pmd.contra.prep(unit) = pmd(unit).p_mod.contra.prep;
    p_mod_pmd.ipsi.move(unit) = pmd(unit).p_mod.ipsi.move;
    p_mod_pmd.contra.move(unit) = pmd(unit).p_mod.contra.move;
    
    if pmd(unit).hem == 0
        modulation_pmd_ic.rest(unit,:) = pmd(unit).modulation.rest;
        modulation_pmd_ic.prep(unit,:) = pmd(unit).modulation.prep;
        modulation_pmd_ic.move(unit,:) = pmd(unit).modulation.move;
        
        arm_pref_pmd_ic.rest(unit,:) = pmd(unit).arm_pref.rest;
        arm_pref_pmd_ic.prep(unit,:) = pmd(unit).arm_pref.prep;
        arm_pref_pmd_ic.move(unit,:) = pmd(unit).arm_pref.move;
    elseif pmd(unit).hem == 1
        arm_pref_pmd_ic.rest(unit,:) = -pmd(unit).arm_pref.rest;
        arm_pref_pmd_ic.prep(unit,:) = -pmd(unit).arm_pref.prep;
        arm_pref_pmd_ic.move(unit,:) = -pmd(unit).arm_pref.move;
        
        modulation_pmd_ic.rest(unit,:) = fliplr(pmd(unit).modulation.rest);
        modulation_pmd_ic.prep(unit,:) = fliplr(pmd(unit).modulation.prep);
        modulation_pmd_ic.move(unit,:) = fliplr(pmd(unit).modulation.move);
    end
end

% clear unit_data since it is very large
% clearvars m1 pmd


%% plot empirical CDF's of modulation for each arm during each task phase
figure('Name', ...
    'modulation distributions across brain area and task epoch')
set(gcf, 'Position',  [2000, 200, 500, 600])

% PMd
subplot(2,1,1)
[pmd_contra_move_ecdf, x_contra] = ecdf(-log10(modulation_pmd_ic.move(:,2)));
[pmd_ipsi_move_ecdf, x_ipsi] = ecdf(-log10(modulation_pmd_ic.move(:,1)));
plot(x_contra, pmd_contra_move_ecdf, 'b:')
hold on
plot(x_ipsi, pmd_ipsi_move_ecdf, 'r:')
[pmd_contra_prep_ecdf, x_contra] = ecdf(-log10(modulation_pmd_ic.prep(:,2)));
[pmd_ipsi_prep_ecdf, x_ipsi] = ecdf(-log10(modulation_pmd_ic.prep(:,1)));
plot(x_contra,pmd_contra_prep_ecdf, 'b--')
plot(x_ipsi, pmd_ipsi_prep_ecdf, 'r--')
[pmd_contra_rest_ecdf, x_contra] = ecdf(-log10(modulation_pmd_ic.rest(:,2)));
[pmd_ipsi_rest_ecdf, x_ipsi] = ecdf(-log10(modulation_pmd_ic.rest(:,1)));
plot(x_contra,pmd_contra_rest_ecdf, 'b')
plot(x_ipsi, pmd_ipsi_rest_ecdf, 'r')
xlim([-2, 1])
xlabel('modulation Depth')
ylabel('F')
title('PMd')

% M1
subplot(2,1,2)
[m1_contra_move_ecdf, x_contra] = ecdf(-log10(modulation_m1_ic.move(:,2)));
[m1_ipsi_move_ecdf, x_ipsi] = ecdf(-log10(modulation_m1_ic.move(:,1)));
plot(x_contra,m1_contra_move_ecdf, 'b:')
hold on
plot(x_ipsi, m1_ipsi_move_ecdf, 'r:')
[m1_contra_prep_ecdf, x_contra] = ecdf(-log10(modulation_m1_ic.prep(:,2)));
[m1_ipsi_prep_ecdf, x_ipsi] = ecdf(-log10(modulation_m1_ic.prep(:,1)));
plot(x_contra,m1_contra_prep_ecdf, 'b--')
plot(x_ipsi, m1_ipsi_prep_ecdf, 'r--')
[m1_contra_rest_ecdf, x_contra] = ecdf(-log10(modulation_m1_ic.rest(:,2)));
[m1_ipsi_rest_ecdf, x_ipsi] = ecdf(-log10(modulation_m1_ic.rest(:,1)));
plot(x_contra,m1_contra_rest_ecdf, 'b')
plot(x_ipsi, m1_ipsi_rest_ecdf, 'r')
xlim([-2, 1])
xlabel('modulation Depth')
ylabel('F')
title('M1')

legend({'Contra Move','Ipsi Move','Contra Instruct','Ipsi Instruct',...
    'Contra Rest','Ipsi Rest'})


%% plot 2x3 subplot of arm pref histograms, separating PMd/M1 and Epochs

figure('Name', ...
    'Arm preference distributions across brain area and task epoch')
set(gcf, 'Position',  [2000, 200, 1000, 600])

% PMd Rest
subplot(2,3,1)
histogram(arm_pref_pmd_ic.rest, -1:0.1:1, 'Normalization', 'probability')
hold on
plot([0,0],[0,0.25],'Color',[0.5,0.5,0.5])
plot(nanmean(arm_pref_pmd_ic.rest), 0,  'r.', 'MarkerSize',20)
xticks([-1,0,1])
yticks([0,0.25])
xlim([-1.1,1.1])
ylim([0,0.25])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
ylabel('PMd','fontweight','bold','fontsize',24)
title('Rest','fontsize',24)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')

% PMd Prep
subplot(2,3,2)
histogram(arm_pref_pmd_ic.prep, -1:0.1:1, 'Normalization', 'probability')
hold on
plot([0,0],[0,0.25],'Color',[0.5,0.5,0.5])
plot(nanmean(arm_pref_pmd_ic.prep), 0,  'r.', 'MarkerSize',20)
xticks([-1,0,1])
yticks([0,0.25])
xlim([-1.1,1.1])
ylim([0,0.25])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
title('Prep','fontsize',24)


% PMd Move
subplot(2,3,3)
histogram(arm_pref_pmd_ic.move, -1:0.1:1, 'Normalization', 'probability')
hold on
plot([0,0],[0,0.25],'Color',[0.5,0.5,0.5])
plot(nanmean(arm_pref_pmd_ic.move), 0,  'r.', 'MarkerSize',20)
xticks([-1,0,1])
yticks([0,0.25])
xlim([-1.1,1.1])
ylim([0,0.25])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
title('Move','fontsize',24)

% M1 Rest
subplot(2,3,4)
histogram(arm_pref_m1_ic.rest, -1:0.1:1, 'Normalization', 'probability')
hold on
plot([0,0],[0,0.25],'Color',[0.5,0.5,0.5])
plot(nanmean(arm_pref_m1_ic.rest), 0,  'r.', 'MarkerSize',20)
xticks([-1,0,1])
yticks([0,0.25])
xlim([-1.1,1.1])
ylim([0,0.25])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
ylabel('M1','fontweight','bold','fontsize',24)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
xlabel('Arm Preference','fontsize',16)

% M1 Prep
subplot(2,3,5)
histogram(arm_pref_m1_ic.prep, -1:0.1:1, 'Normalization', 'probability')
hold on
plot([0,0],[0,0.25],'Color',[0.5,0.5,0.5])
plot(nanmean(arm_pref_m1_ic.prep), 0,  'r.', 'MarkerSize',20)
xticks([-1,0,1])
yticks([0,0.25])
xlim([-1.1,1.1])
ylim([0,0.25])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)

% M1 Move
subplot(2,3,6)
histogram(arm_pref_m1_ic.move, -1:0.1:1, 'Normalization', 'probability')
hold on
plot([0,0],[0,0.25],'Color',[0.5,0.5,0.5])
plot(nanmean(arm_pref_m1_ic.move), 0,  'r.', 'MarkerSize',20)
xticks([-1,0,1])
yticks([0,0.25])
xlim([-1.1,1.1])
ylim([0,0.25])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)


%% plot empirical CDF's for arm pref during each task phase

figure('Name', 'M1 - Arm pref ECDF')
[m1_rest_ecdf, x_rest] = ecdf(arm_pref_m1_ic.rest);
rest_uq = quantile(arm_pref_m1_ic.rest, 0.75);
[m1_prep_ecdf, x_prep] = ecdf(arm_pref_m1_ic.prep);
prep_uq = quantile(arm_pref_m1_ic.prep, 0.75);
[m1_move_ecdf, x_move] = ecdf(arm_pref_m1_ic.move);
move_uq = quantile(arm_pref_m1_ic.move, 0.75);
plot(x_rest,m1_rest_ecdf, 'b')
hold on
plot(x_prep, m1_prep_ecdf, 'r')
plot(x_move, m1_move_ecdf, 'g')
plot([rest_uq,rest_uq],[0,0.75],'b-o');
plot([prep_uq,prep_uq],[0,0.75],'r-o');
plot([move_uq,move_uq],[0,0.75],'g-o');
xticks(-1:0.5:1)
yticks(0:0.25:1)
grid on
xlabel('Arm preference')
ylabel({'Cumulative','proportion of units'})
title('M1')

figure('Name', 'PMd - Arm pref ECDF')
[pmd_rest_ecdf, x_rest] = ecdf(arm_pref_pmd_ic.rest);
rest_uq = quantile(arm_pref_pmd_ic.rest, 0.75);
[pmd_prep_ecdf, x_prep] = ecdf(arm_pref_pmd_ic.prep);
prep_uq = quantile(arm_pref_pmd_ic.prep, 0.75);
[pmd_move_ecdf, x_move] = ecdf(arm_pref_pmd_ic.move);
move_uq = quantile(arm_pref_pmd_ic.move, 0.75);
plot(x_rest,pmd_rest_ecdf, 'b')
hold on
plot(x_prep, pmd_prep_ecdf, 'r')
plot(x_move, pmd_move_ecdf, 'g')
plot([rest_uq,rest_uq],[0,0.75],'b-o');
plot([prep_uq,prep_uq],[0,0.75],'r-o');
plot([move_uq,move_uq],[0,0.75],'g-o');
xticks(-1:0.5:1)
yticks(0:0.25:1)
grid on
xlabel('Arm preference')
ylabel({'Cumulative','proportion of units'})
title('PMd')


%% plot 2x3 subplot of ABS arm pref histograms, separating PMd/M1 and Epoch

figure('Name', ...
    'Absolute arm preference distributions across brain area and task epoch')
set(gcf, 'Position',  [2000, 200, 1000, 600])

% PMd Rest
subplot(2,3,1)
histogram(abs(arm_pref_pmd_ic.rest), 0:0.05:1, 'Normalization', 'probability')
hold on
plot(nanmean(abs(arm_pref_pmd_ic.rest)), 0,  'r.', 'MarkerSize',20)
xticks([0,0.5,1])
yticks([0,0.25])
xlim([-0.05,1.05])
ylim([0,0.25])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
ylabel('PMd','fontweight','bold','fontsize',24)
title('Rest','fontsize',24)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')

% PMd Prep
subplot(2,3,2)
histogram(abs(arm_pref_pmd_ic.prep), 0:0.05:1, 'Normalization', 'probability')
hold on
plot(nanmean(abs(arm_pref_pmd_ic.prep)), 0,  'r.', 'MarkerSize',20)
xticks([0,0.5,1])
yticks([0,0.25])
xlim([-0.05,1.05])
ylim([0,0.25])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
title('Prep','fontsize',24)


% PMd Move
subplot(2,3,3)
histogram(abs(arm_pref_pmd_ic.move), 0:0.05:1, 'Normalization', 'probability')
hold on
plot(nanmean(abs(arm_pref_pmd_ic.move)), 0,  'r.', 'MarkerSize',20)
xticks([0,0.5,1])
yticks([0,0.25])
xlim([-0.05,1.05])
ylim([0,0.25])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
title('Move','fontsize',24)

% M1 Rest
subplot(2,3,4)
histogram(abs(arm_pref_m1_ic.rest), 0:0.05:1, 'Normalization', 'probability')
hold on
plot(nanmean(abs(arm_pref_m1_ic.rest)), 0,  'r.', 'MarkerSize',20)
xticks([0,0.5,1])
yticks([0,0.25])
xlim([-0.05,1.05])
ylim([0,0.25])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
ylabel('M1','fontweight','bold','fontsize',24)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
xlabel('Absolute Arm Preference','fontsize',16)

% M1 Prep
subplot(2,3,5)
histogram(abs(arm_pref_m1_ic.prep), 0:0.05:1, 'Normalization', 'probability')
hold on
plot(nanmean(abs(arm_pref_m1_ic.prep)), 0,  'r.', 'MarkerSize',20)
xticks([0,0.5,1])
yticks([0,0.25])
xlim([-0.05,1.05])
ylim([0,0.25])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)

% M1 Move
subplot(2,3,6)
histogram(abs(arm_pref_m1_ic.move), 0:0.05:1, 'Normalization', 'probability')
hold on
plot(nanmean(abs(arm_pref_m1_ic.move)), 0,  'r.', 'MarkerSize',20)
xticks([0,0.5,1])
yticks([0,0.25])
xlim([-0.05,1.05])
ylim([0,0.25])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)


%% plot empirical CDF's for ABS arm pref during each task phase

figure('Name', 'M1 - ABS arm pref ECDF')
[m1_rest_ecdf, x_rest] = ecdf(abs(arm_pref_m1_ic.rest));
rest_uq = quantile(abs(arm_pref_m1_ic.rest), 0.75);
[m1_prep_ecdf, x_prep] = ecdf(abs(arm_pref_m1_ic.prep));
prep_uq = quantile(abs(arm_pref_m1_ic.prep), 0.75);
[m1_move_ecdf, x_move] = ecdf(abs(arm_pref_m1_ic.move));
move_uq = quantile(abs(arm_pref_m1_ic.move), 0.75);
plot(x_rest,m1_rest_ecdf, 'b')
hold on
plot(x_prep, m1_prep_ecdf, 'r')
plot(x_move, m1_move_ecdf, 'g')
plot([rest_uq,rest_uq],[0,0.75],'b-o');
plot([prep_uq,prep_uq],[0,0.75],'r-o');
plot([move_uq,move_uq],[0,0.75],'g-o');
xticks(0:0.25:1)
yticks(0:0.25:1)
grid on
xlabel({'Absolute','Arm Preference'})
ylabel({'Cumulative','proportion of units'})
title('M1')

figure('Name', 'PMd - ABS arm pref ECDF')
[pmd_rest_ecdf, x_rest] = ecdf(abs(arm_pref_pmd_ic.rest));
rest_uq = quantile(abs(arm_pref_pmd_ic.rest), 0.75);
[pmd_prep_ecdf, x_prep] = ecdf(abs(arm_pref_pmd_ic.prep));
prep_uq = quantile(abs(arm_pref_pmd_ic.prep), 0.75);
[pmd_move_ecdf, x_move] = ecdf(abs(arm_pref_pmd_ic.move));
move_uq = quantile(abs(arm_pref_pmd_ic.move), 0.75);
plot(x_rest,pmd_rest_ecdf, 'b')
hold on
plot(x_prep, pmd_prep_ecdf, 'r')
plot(x_move, pmd_move_ecdf, 'g')
plot([rest_uq,rest_uq],[0,0.75],'b-o');
plot([prep_uq,prep_uq],[0,0.75],'r-o');
plot([move_uq,move_uq],[0,0.75],'g-o');
xticks(0:0.25:1)
yticks(0:0.25:1)
grid on
xlabel({'Absolute','Arm Preference'})
ylabel({'Cumulative','proportion of units'})
title('PMd')


%% plot mean +/- bootstrapped 95% CI for the distributions at each epoch

% run bootstrapping procedure over 10,000 bootstraps of same sized samples
n_m1_ipsi_pref = length(arm_pref_m1_ic.move);
n_pmd_ipsi_pref = length(arm_pref_pmd_ic.move);

for boot = 10000:-1:1
    
    boot_sample = ...
        datasample((arm_pref_m1_ic.rest), n_m1_ipsi_pref, 'Replace',true);
    boot_mean_m1_rest(boot) = quantile(boot_sample, 0.75);
    boot_sample = ...
        datasample((arm_pref_m1_ic.prep), n_m1_ipsi_pref, 'Replace',true);
    boot_mean_m1_prep(boot) = quantile(boot_sample, 0.75);
    boot_sample = ...
        datasample((arm_pref_m1_ic.move), n_m1_ipsi_pref, 'Replace',true);
    boot_mean_m1_move(boot) = quantile(boot_sample, 0.75);
    
    boot_sample = ...
        datasample((arm_pref_pmd_ic.rest), n_pmd_ipsi_pref, 'Replace',true);
    boot_mean_pmd_rest(boot) = quantile(boot_sample, 0.75);
    boot_sample = ...
        datasample((arm_pref_pmd_ic.prep), n_pmd_ipsi_pref, 'Replace',true);
    boot_mean_pmd_prep(boot) = quantile(boot_sample, 0.75);
    boot_sample = ...
        datasample((arm_pref_pmd_ic.move), n_pmd_ipsi_pref, 'Replace',true);
    boot_mean_pmd_move(boot) = quantile(boot_sample, 0.75);
    
end

% sort results
boot_mean_m1_rest = sort(boot_mean_m1_rest);
boot_mean_m1_prep = sort(boot_mean_m1_prep);
boot_mean_m1_move = sort(boot_mean_m1_move);

boot_mean_pmd_rest = sort(boot_mean_pmd_rest);
boot_mean_pmd_prep = sort(boot_mean_pmd_prep);
boot_mean_pmd_move = sort(boot_mean_pmd_move);

% find center of the bootstrapped distribution
Y_m1 = [mean(boot_mean_m1_rest), mean(boot_mean_m1_prep), ...
    mean(boot_mean_m1_move)];
Y_pmd = [mean(boot_mean_pmd_rest), mean(boot_mean_pmd_prep), ...
    mean(boot_mean_pmd_move)];

% find 95% CI
neg_m1 = Y_m1 - [boot_mean_m1_rest(26), boot_mean_m1_prep(26), ...
    boot_mean_m1_move(26)];
neg_pmd = Y_pmd - [boot_mean_pmd_rest(26), boot_mean_pmd_prep(26), ...
    boot_mean_pmd_move(26)];
pos_m1 = [boot_mean_m1_rest(end-25), boot_mean_m1_prep(end-25), ...
    boot_mean_m1_move(end-25)] - Y_m1;
pos_pmd = [boot_mean_pmd_rest(end-25), boot_mean_pmd_prep(end-25), ...
    boot_mean_pmd_move(end-25)] - Y_pmd;

% plot

figure('Name', 'Lateral bias across epochs')
set(gcf, 'Position',  [2000, 200, 500, 250])
errorbar([1:3]-0.05, Y_pmd, neg_pmd, pos_pmd, 'o-', 'LineWidth',1.5)
hold on
errorbar([1:3]+0.05, Y_m1, neg_m1, pos_m1, 'o-', 'LineWidth',1.5)
xlim([0.5, 3.5])
xticks([1,2,3])
xticklabels({'Rest','Prep','Move'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
yticks([0,0.7])
ylabel({'Arm preference upper quartile', 'mean +/- 95% CI'},'fontsize',16)
ylim([0,0.7])
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
legend({'PMd','M1'})


%% plot mean +/- bootstrapped 95% CI for the ABS distributions at each epoch

% run bootstrapping procedure over 10,000 bootstraps of same sized samples
n_m1_ipsi_pref = length(arm_pref_m1_ic.move);
n_pmd_ipsi_pref = length(arm_pref_pmd_ic.move);

for boot = 10000:-1:1
    
    boot_sample = ...
        datasample(abs(arm_pref_m1_ic.rest), n_m1_ipsi_pref, 'Replace',true);
    boot_mean_m1_rest(boot) = quantile(boot_sample, 0.75);
    boot_sample = ...
        datasample(abs(arm_pref_m1_ic.prep), n_m1_ipsi_pref, 'Replace',true);
    boot_mean_m1_prep(boot) = quantile(boot_sample, 0.75);
    boot_sample = ...
        datasample(abs(arm_pref_m1_ic.move), n_m1_ipsi_pref, 'Replace',true);
    boot_mean_m1_move(boot) = quantile(boot_sample, 0.75);
    
    boot_sample = ...
        datasample(abs(arm_pref_pmd_ic.rest), n_pmd_ipsi_pref, 'Replace',true);
    boot_mean_pmd_rest(boot) = quantile(boot_sample, 0.75);
    boot_sample = ...
        datasample(abs(arm_pref_pmd_ic.prep), n_pmd_ipsi_pref, 'Replace',true);
    boot_mean_pmd_prep(boot) = quantile(boot_sample, 0.75);
    boot_sample = ...
        datasample(abs(arm_pref_pmd_ic.move), n_pmd_ipsi_pref, 'Replace',true);
    boot_mean_pmd_move(boot) = quantile(boot_sample, 0.75);
    
end

% sort results
boot_mean_m1_rest = sort(boot_mean_m1_rest);
boot_mean_m1_prep = sort(boot_mean_m1_prep);
boot_mean_m1_move = sort(boot_mean_m1_move);

boot_mean_pmd_rest = sort(boot_mean_pmd_rest);
boot_mean_pmd_prep = sort(boot_mean_pmd_prep);
boot_mean_pmd_move = sort(boot_mean_pmd_move);

% find center of the bootstrapped distribution
Y_m1 = [mean(boot_mean_m1_rest), mean(boot_mean_m1_prep), ...
    mean(boot_mean_m1_move)];
Y_pmd = [mean(boot_mean_pmd_rest), mean(boot_mean_pmd_prep), ...
    mean(boot_mean_pmd_move)];

% find 95% CI
neg_m1 = Y_m1 - [boot_mean_m1_rest(26), boot_mean_m1_prep(26), ...
    boot_mean_m1_move(26)];
neg_pmd = Y_pmd - [boot_mean_pmd_rest(26), boot_mean_pmd_prep(26), ...
    boot_mean_pmd_move(26)];
pos_m1 = [boot_mean_m1_rest(end-25), boot_mean_m1_prep(end-25), ...
    boot_mean_m1_move(end-25)] - Y_m1;
pos_pmd = [boot_mean_pmd_rest(end-25), boot_mean_pmd_prep(end-25), ...
    boot_mean_pmd_move(end-25)] - Y_pmd;

% plot

figure('Name', 'Arm dedication across epochs')
set(gcf, 'Position',  [2000, 200, 500, 250])
errorbar([1:3]-0.05, Y_pmd, neg_pmd, pos_pmd, 'o-', 'LineWidth',1.5)
hold on
errorbar([1:3]+0.05, Y_m1, neg_m1, pos_m1, 'o-', 'LineWidth',1.5)
xlim([0.5, 3.5])
xticks([1,2,3])
xticklabels({'Rest','Prep','Move'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
ylabel({'|Arm preference|', 'mean +/- 95% CI'},'fontsize',16)
ylim([0,0.7])
yticks([0,0.7])
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
legend({'PMd','M1'})


%% plot scatter of arm pref vs preferred limb MD

figure('Name', 'Arm dedicated units tend to be more deeply modulated')
set(gcf, 'Position',  [2000, 200, 1000, 600])

% PMd Rest
subplot(2,3,1)
% plot scatter
X = arm_pref_pmd_ic.rest;
Y = log10(max(modulation_pmd_ic.rest,[],2)); %%%%
scatter(X, Y)
hold on
% fit line (log-linear space) and plot on top of scatter
p = polyfit(X(X<0),Y(X<0),1);
f = polyval(p,[-1,0]);
plot([-1,0],f,'r', 'LineWidth',1.5)
p = polyfit(X(X>0),Y(X>0),1);
f = polyval(p,[0,1]);
plot([0,1],f,'r', 'LineWidth',1.5)
% format display
xticks([-1,0,1])
yticks([0,1,2])
xlim([-1.1,1.1])
ylim([-1,2.5])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
yticklabels({'10^0','10^1','10^2'})
ylabel('PMd','fontweight','bold','fontsize',24)
title('Rest','fontsize',24)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')

% PMd prep
subplot(2,3,2)
% plot scatter
X = arm_pref_pmd_ic.prep;
Y = log10(max(modulation_pmd_ic.prep,[],2));
scatter(X, Y)
hold on
% fit line (log-linear space) and plot on top of scatter
p = polyfit(X(X<0),Y(X<0),1);
f = polyval(p,[-1,0]);
plot([-1,0],f,'r', 'LineWidth',1.5)
p = polyfit(X(X>0),Y(X>0),1);
f = polyval(p,[0,1]);
plot([0,1],f,'r', 'LineWidth',1.5)
% format display
xticks([-1,0,1])
yticks([-2,0,2])
xlim([-1.1,1.1])
ylim([-2,3])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
yticklabels({'10^0','10^1','10^2'})
title('Prep','fontsize',24)

% PMd move
subplot(2,3,3)
% plot scatter
X = arm_pref_pmd_ic.move;
Y = log10(max(modulation_pmd_ic.move,[],2));
scatter(X, Y)
hold on
% fit line (log-linear space) and plot on top of scatter
p = polyfit(X(X<0),Y(X<0),1);
f = polyval(p,[-1,0]);
plot([-1,0],f,'r', 'LineWidth',1.5)
p = polyfit(X(X>0),Y(X>0),1);
f = polyval(p,[0,1]);
plot([0,1],f,'r', 'LineWidth',1.5)
% format display
xticks([-1,0,1])
yticks([-2,0,2])
xlim([-1.1,1.1])
ylim([-2,3])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
yticklabels({'10^0','10^1','10^2'})
title('Move','fontsize',24)

% m1 Rest
subplot(2,3,4)
% plot scatter
X = arm_pref_m1_ic.rest;
Y = log10(max(modulation_m1_ic.rest,[],2));
scatter(X, Y)
hold on
% fit line (log-linear space) and plot on top of scatter
p = polyfit(X(X<0),Y(X<0),1);
f = polyval(p,[-1,0]);
plot([-1,0],f,'r', 'LineWidth',1.5)
p = polyfit(X(X>0),Y(X>0),1);
f = polyval(p,[0,1]);
plot([0,1],f,'r', 'LineWidth',1.5)
% format display
xticks([-1,0,1])
yticks([-2,0,2])
xlim([-1.1,1.1])
ylim([-2,3])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
yticklabels({'10^0','10^1','10^2'})
ylabel('M1','fontweight','bold','fontsize',24)
title('Rest','fontsize',24)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
xlabel('Arm Preference','fontsize',16)

% M1 prep
subplot(2,3,5)
% plot scatter
X = arm_pref_m1_ic.prep;
Y = log10(max(modulation_m1_ic.prep,[],2));
scatter(X, Y)
hold on
% fit line (log-linear space) and plot on top of scatter
p = polyfit(X(X<0),Y(X<0),1);
f = polyval(p,[-1,0]);
plot([-1,0],f,'r', 'LineWidth',1.5)
p = polyfit(X(X>0),Y(X>0),1);
f = polyval(p,[0,1]);
plot([0,1],f,'r', 'LineWidth',1.5)
% format display
xticks([-1,0,1])
yticks([-2,0,2])
xlim([-1.1,1.1])
ylim([-2,3])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
yticklabels({'10^0','10^1','10^2'})
title('Prep','fontsize',24)

% m1 move
subplot(2,3,6)
% plot scatter
X = arm_pref_m1_ic.move;
Y = log10(max(modulation_m1_ic.move,[],2));
scatter(X, Y)
hold on
% fit line (log-linear space) and plot on top of scatter
p = polyfit(X(X<0),Y(X<0),1);
f = polyval(p,[-1,0]);
plot([-1,0],f,'r', 'LineWidth',1.5)
p = polyfit(X(X>0),Y(X>0),1);
f = polyval(p,[0,1]);
plot([0,1],f,'r', 'LineWidth',1.5)
% format display
xticks([-1,0,1])
yticks([-2,0,2])
xlim([-1.1,1.1])
ylim([-2,3])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
yticklabels({'10^0','10^1','10^2'})
title('Move','fontsize',24)


%% plot scatter of absolute arm pref vs preferred limb MD

figure('Name', 'Arm dedicated units tend to be more deeply modulated')
set(gcf, 'Position',  [2000, 200, 1000, 600])

% PMd Rest
subplot(2,3,1)
% plot scatter
X = abs(arm_pref_pmd_ic.rest);
Y = log10(max(modulation_pmd_ic.rest,[],2));
idx = ~isnan(X)&Y>-1;
scatter(X(idx), Y(idx))
hold on
% fit line (log-linear space) and plot on top of scatter
p = polyfit(X(idx),Y(idx),1);
f = polyval(p,[-1,1]);
plot([-1,1],f,'r', 'LineWidth',1.5)
% format display
xticks([0,0.5,1])
yticks([0,1,2])
xlim([-0.05,1.05])
ylim([-1,2.5])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
yticklabels({'10^0','10^1','10^2'})
ylabel('PMd','fontweight','bold','fontsize',24)
title('Rest','fontsize',24)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
grid on

% PMd prep
subplot(2,3,2)
% plot scatter
X = abs(arm_pref_pmd_ic.prep);
Y = log10(max(modulation_pmd_ic.prep,[],2));
idx = ~isnan(X)&Y>-1;
scatter(X(idx), Y(idx))
hold on
% fit line (log-linear space) and plot on top of scatter
p = polyfit(X(idx),Y(idx),1);
f = polyval(p,[-1,1]);
plot([-1,1],f,'r', 'LineWidth',1.5)
% format display
xticks([0,0.5,1])
yticks([0,1,2])
xlim([-0.05,1.05])
ylim([-1,2.5])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
yticklabels({'10^0','10^1','10^2'})
title('Prep','fontsize',24)
grid on

% PMd move
subplot(2,3,3)
% plot scatter
X = abs(arm_pref_pmd_ic.move);
Y = log10(max(modulation_pmd_ic.move,[],2));
idx = ~isnan(X)&Y>-1;
scatter(X(idx), Y(idx))
hold on
% fit line (log-linear space) and plot on top of scatter
p = polyfit(X(idx),Y(idx),1);
f = polyval(p,[-1,1]);
plot([-1,1],f,'r', 'LineWidth',1.5)
% format display
xticks([0,0.5,1])
yticks([0,1,2])
xlim([-0.05,1.05])
ylim([-1,2.5])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
yticklabels({'10^0','10^1','10^2'})
title('Move','fontsize',24)
grid on

% m1 Rest
subplot(2,3,4)
% plot scatter
X = abs(arm_pref_m1_ic.rest);
Y = log10(max(modulation_m1_ic.rest,[],2));
idx = ~isnan(X)&Y>-1;
scatter(X(idx), Y(idx))
hold on
% fit line (log-linear space) and plot on top of scatter
p = polyfit(X(idx),Y(idx),1);
f = polyval(p,[-1,1]);
plot([-1,1],f,'r', 'LineWidth',1.5)
% format display
xticks([0,0.5,1])
yticks([0,1,2])
xlim([-0.05,1.05])
ylim([-1,2.5])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
yticklabels({'10^0','10^1','10^2'})
ylabel('M1','fontweight','bold','fontsize',24)
title('Rest','fontsize',24)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
xlabel('Arm Preference','fontsize',16)
grid on

% M1 prep
subplot(2,3,5)
% plot scatter
X = abs(arm_pref_m1_ic.prep);
Y = log10(max(modulation_m1_ic.prep,[],2));
idx = ~isnan(X)&Y>-1;
scatter(X(idx), Y(idx))
hold on
% fit line (log-linear space) and plot on top of scatter
p = polyfit(X(idx),Y(idx),1);
f = polyval(p,[-1,1]);
plot([-1,1],f,'r', 'LineWidth',1.5)
% format display
xticks([0,0.5,1])
yticks([0,1,2])
xlim([-0.05,1.05])
ylim([-1,2.5])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
yticklabels({'10^0','10^1','10^2'})
title('Prep','fontsize',24)
grid on

% m1 move
subplot(2,3,6)
% plot scatter
X = abs(arm_pref_m1_ic.move);
Y = log10(max(modulation_m1_ic.move,[],2));
idx = ~isnan(X)&Y>-1;
scatter(X(idx), Y(idx))
hold on
% fit line (log-linear space) and plot on top of scatter
p = polyfit(X(idx),Y(idx),1);
f = polyval(p,[-1,1]);
plot([-1,1],f,'r', 'LineWidth',1.5)
% format display
xticks([0,0.5,1])
yticks([0,1,2])
xlim([-0.05,1.05])
ylim([-1,2.5])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',14)
yticklabels({'10^0','10^1','10^2'})
title('Move','fontsize',24)
grid on


%% plot slopes across epochs, bootstrap

% run bootstrapping procedure over 10,000 bootstraps of same sized samples

for boot = 10000:-1:1
    
    % pmd
    X = arm_pref_pmd_ic.rest;
    Y = log10(max(modulation_pmd_ic.rest,[],2));
    idx = X<0 & ~isnan(X) & Y>-1;
    Y = Y(idx);
    X = X(idx);
    [X, idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(idx);
    b = polyfit(X,Y,1);
    slope_pmd_ipsi(boot,1) = b(1);
    
    X = arm_pref_pmd_ic.rest;
    Y = log10(max(modulation_pmd_ic.rest,[],2));
    idx = X>0 & ~isnan(X) & Y>-1;
    Y = Y(idx);
    X = X(idx);
    [X, idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(idx);
    b = polyfit(X,Y,1);
    slope_pmd_contra(boot,1) = b(1);
    
    X = arm_pref_pmd_ic.prep;
    Y = log10(max(modulation_pmd_ic.prep,[],2));
    idx = X<0 & ~isnan(X) & Y>-1;
    Y = Y(idx);
    X = X(idx);
    [X, idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(idx);
    b = polyfit(X,Y,1);
    slope_pmd_ipsi(boot,2) = b(1);
    
    X = arm_pref_pmd_ic.prep;
    Y = log10(max(modulation_pmd_ic.prep,[],2));
    idx = X>0 & ~isnan(X) & Y>-1;
    Y = Y(idx);
    X = X(idx);
    [X, idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(idx);
    b = polyfit(X,Y,1);
    slope_pmd_contra(boot,2) = b(1);
    
    X = arm_pref_pmd_ic.move;
    Y = log10(max(modulation_pmd_ic.move,[],2));
    idx = X<0 & ~isnan(X) & Y>-1;
    Y = Y(idx);
    X = X(idx);
    [X, idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(idx);
    b = polyfit(X,Y,1);
    slope_pmd_ipsi(boot,3) = b(1);
    
    X = arm_pref_pmd_ic.move;
    Y = log10(max(modulation_pmd_ic.move,[],2));
    idx = X>0 & ~isnan(X) & Y>-1;
    Y = Y(idx);
    X = X(idx);
    [X, idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(idx);
    b = polyfit(X,Y,1);
    slope_pmd_contra(boot,3) = b(1);
    
    % m1
    X = arm_pref_m1_ic.rest;
    Y = log10(max(modulation_m1_ic.rest,[],2));
    idx = X<0 & ~isnan(X) & Y>-1;
    Y = Y(idx);
    X = X(idx);   
    [X, idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(idx);
    b = polyfit(X,Y,1);
    slope_m1_ipsi(boot,1) = b(1);
    
    X = arm_pref_m1_ic.rest;
    Y = log10(max(modulation_m1_ic.rest,[],2));
    idx = X>0 & ~isnan(X) & Y>-1;
    Y = Y(idx);
    X = X(idx);
    [X, idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(idx);
    b = polyfit(X,Y,1);
    slope_m1_contra(boot,1) = b(1);
    
    X = arm_pref_m1_ic.prep;
    Y = log10(max(modulation_m1_ic.prep,[],2));
    idx = X<0 & ~isnan(X) & Y>-1;
    Y = Y(idx);
    X = X(idx);
    [X, idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(idx);
    b = polyfit(X,Y,1);
    slope_m1_ipsi(boot,2) = b(1);
    
    X = arm_pref_m1_ic.prep;
    Y = log10(max(modulation_m1_ic.prep,[],2));
    idx = X>0 & ~isnan(X) & Y>-1;
    Y = Y(idx);
    X = X(idx);
    [X, idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(idx);
    b = polyfit(X,Y,1);
    slope_m1_contra(boot,2) = b(1);
    
    X = arm_pref_m1_ic.move;
    Y = log10(max(modulation_m1_ic.move,[],2));
    idx = X<0 & ~isnan(X) & Y>-1;
    Y = Y(idx);
    X = X(idx);
    [X, idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(idx);
    b = polyfit(X,Y,1);
    slope_m1_ipsi(boot,3) = b(1);
    
    X = arm_pref_m1_ic.move;
    Y = log10(max(modulation_m1_ic.move,[],2));
    idx = X>0 & ~isnan(X) & Y>-1;
    Y = Y(idx);
    X = X(idx);
    [X, idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(idx);
    b = polyfit(X,Y,1);
    slope_m1_contra(boot,3) = b(1);
    
end

% sort results
slope_pmd_ipsi = [sort(-slope_pmd_ipsi(:,1),1),...
    sort(-slope_pmd_ipsi(:,2),1),...
    sort(-slope_pmd_ipsi(:,3),1)];
slope_pmd_contra = [sort(slope_pmd_contra(:,1),1),...
    sort(slope_pmd_contra(:,2),1),...
    sort(slope_pmd_contra(:,3),1)];
slope_m1_ipsi = [sort(-slope_m1_ipsi(:,1),1),...
    sort(-slope_m1_ipsi(:,2),1),...
    sort(-slope_m1_ipsi(:,3),1)];
slope_m1_contra = [sort(slope_m1_contra(:,1),1),...
    sort(slope_m1_contra(:,2),1),...
    sort(slope_m1_contra(:,3),1)];

% find center of the bootstrapped distribution
Y_pmd_ipsi = mean(slope_pmd_ipsi);
Y_pmd_contra = mean(slope_pmd_contra);
Y_m1_ipsi = mean(slope_m1_ipsi);
Y_m1_contra = mean(slope_m1_contra);

% find 95% CI
neg_pmd_ipsi = Y_pmd_ipsi - quantile(slope_pmd_ipsi, 0.025);
neg_pmd_contra = Y_pmd_contra - quantile(slope_pmd_contra, 0.025);
neg_m1_ipsi = Y_m1_ipsi - quantile(slope_m1_ipsi, 0.025);
neg_m1_contra = Y_m1_contra - quantile(slope_m1_contra, 0.025);

pos_pmd_ipsi = -Y_pmd_ipsi + quantile(slope_pmd_ipsi, 0.975);
pos_pmd_contra = -Y_pmd_contra + quantile(slope_pmd_contra, 0.975);
pos_m1_ipsi = -Y_m1_ipsi + quantile(slope_m1_ipsi, 0.975);
pos_m1_contra = -Y_m1_contra + quantile(slope_m1_contra, 0.975);

% permutation test for testing differences across task phases
% permute random resamples of each pairwise grouping
slope_pmd_ipsi_rp = zeros(10000,1);
slope_pmd_ipsi_rm = zeros(10000,1);
slope_pmd_ipsi_pm = zeros(10000,1);
slope_pmd_contra_rp = zeros(10000,1);
slope_pmd_contra_rm = zeros(10000,1);
slope_pmd_contra_pm = zeros(10000,1);
slope_m1_ipsi_rp = zeros(10000,1);
slope_m1_ipsi_rm = zeros(10000,1);
slope_m1_ipsi_pm = zeros(10000,1);
slope_m1_contra_rp = zeros(10000,1);
slope_m1_contra_rm = zeros(10000,1);
slope_m1_contra_pm = zeros(10000,1);
parfor perm = 1:10000
    
    % pmd, ipsi
    X_rest = arm_pref_pmd_ic.rest;
    Y_rest = log10(max(modulation_pmd_ic.rest,[],2));
    idx = X_rest<0 & ~isnan(X_rest) & Y_rest>-1;
    X_rest = X_rest(idx);
    Y_rest = Y_rest(idx);
    n_rest = length(X_rest);
    
    X_prep = arm_pref_pmd_ic.prep;
    Y_prep = log10(max(modulation_pmd_ic.prep,[],2));
    idx = X_prep<0 & ~isnan(X_prep) & Y_prep>-1;
    X_prep = X_prep(idx);
    Y_prep = Y_prep(idx);
    n_prep = length(X_prep);
    
    X_move = arm_pref_pmd_ic.move;
    Y_move = log10(max(modulation_pmd_ic.move,[],2));
    idx = X_move<0 & ~isnan(X_move) & Y_move>-1;
    X_move = X_move(idx);
    Y_move = Y_move(idx);
    n_move = length(X_move);
    
    X_rp = [X_rest; X_prep];
    Y_rp = [Y_rest; Y_prep];
    r_idx = randi(length(X_rp), n_rest,1);
    X_rp_r = X_rp(r_idx);
    Y_rp_r = Y_rp(r_idx);
    X_rp_p = X_rp(setdiff(1:length(X_rp),r_idx));
    Y_rp_p = Y_rp(setdiff(1:length(Y_rp),r_idx));
    b_r = polyfit(X_rp_r, Y_rp_r, 1);
    b_p = polyfit(X_rp_p, Y_rp_p, 1);
    slope_pmd_ipsi_rp(perm,1) = b_p(1) - b_r(1);
    
    X_rm = [X_rest; X_move];
    Y_rm = [Y_rest; X_move];
    r_idx = randi(length(X_rm), n_rest,1);
    X_rm_r = X_rm(r_idx);
    Y_rm_r = Y_rm(r_idx);
    X_rm_m = X_rm(setdiff(1:length(X_rm),r_idx));
    Y_rm_m = Y_rm(setdiff(1:length(Y_rm),r_idx));
    b_r = polyfit(X_rm_r, Y_rm_r, 1);
    b_m = polyfit(X_rm_m, Y_rm_m, 1);
    slope_pmd_ipsi_rm(perm,1) = b_m(1) - b_r(1);
    
    X_pm = [X_prep; X_move];
    Y_pm = [Y_prep; X_move];
    p_idx = randi(length(X_pm), n_prep,1);
    X_pm_p = X_pm(p_idx);
    Y_pm_p = Y_pm(p_idx);
    X_pm_m = X_pm(setdiff(1:length(X_pm),p_idx));
    Y_pm_m = Y_pm(setdiff(1:length(Y_pm),p_idx));
    b_p = polyfit(X_pm_p, Y_pm_p, 1);
    b_m = polyfit(X_pm_m, Y_pm_m, 1);
    slope_pmd_ipsi_pm(perm,1) = b_m(1) - b_p(1);
  
    
    % pmd, contra
    X_rest = arm_pref_pmd_ic.rest;
    Y_rest = log10(max(modulation_pmd_ic.rest,[],2));
    idx = X_rest>0 & ~isnan(X_rest) & Y_rest>-1;
    X_rest = X_rest(idx);
    Y_rest = Y_rest(idx);
    n_rest = length(X_rest);
    
    X_prep = arm_pref_pmd_ic.prep;
    Y_prep = log10(max(modulation_pmd_ic.prep,[],2));
    idx = X_prep>0 & ~isnan(X_prep) & Y_prep>-1;
    X_prep = X_prep(idx);
    Y_prep = Y_prep(idx);
    n_prep = length(X_prep);
    
    X_move = arm_pref_pmd_ic.move;
    Y_move = log10(max(modulation_pmd_ic.move,[],2));
    idx = X_move>0 & ~isnan(X_move) & Y_move>-1;
    X_move = X_move(idx);
    Y_move = Y_move(idx);
    n_move = length(X_move);
    
    X_rp = [X_rest; X_prep];
    Y_rp = [Y_rest; Y_prep];
    r_idx = randi(length(X_rp), n_rest,1);
    X_rp_r = X_rp(r_idx);
    Y_rp_r = Y_rp(r_idx);
    X_rp_p = X_rp(setdiff(1:length(X_rp),r_idx));
    Y_rp_p = Y_rp(setdiff(1:length(Y_rp),r_idx));
    b_r = polyfit(X_rp_r, Y_rp_r, 1);
    b_p = polyfit(X_rp_p, Y_rp_p, 1);
    slope_pmd_contra_rp(perm,1) = b_p(1) - b_r(1);
    
    X_rm = [X_rest; X_move];
    Y_rm = [Y_rest; X_move];
    r_idx = randi(length(X_rm), n_rest,1);
    X_rm_r = X_rm(r_idx);
    Y_rm_r = Y_rm(r_idx);
    X_rm_m = X_rm(setdiff(1:length(X_rm),r_idx));
    Y_rm_m = Y_rm(setdiff(1:length(Y_rm),r_idx));
    b_r = polyfit(X_rm_r, Y_rm_r, 1);
    b_m = polyfit(X_rm_m, Y_rm_m, 1);
    slope_pmd_contra_rm(perm,1) = b_m(1) - b_r(1);
    
    X_pm = [X_prep; X_move];
    Y_pm = [Y_prep; X_move];
    p_idx = randi(length(X_pm), n_prep,1);
    X_pm_p = X_pm(p_idx);
    Y_pm_p = Y_pm(p_idx);
    X_pm_m = X_pm(setdiff(1:length(X_pm),p_idx));
    Y_pm_m = Y_pm(setdiff(1:length(Y_pm),p_idx));
    b_p = polyfit(X_pm_p, Y_pm_p, 1);
    b_m = polyfit(X_pm_m, Y_pm_m, 1);
    slope_pmd_contra_pm(perm,1) = b_m(1) - b_p(1);
    
    % m1, ipsi
    X_rest = arm_pref_m1_ic.rest;
    Y_rest = log10(max(modulation_m1_ic.rest,[],2));
    idx = X_rest<0 & ~isnan(X_rest) & Y_rest>-1;
    X_rest = X_rest(idx);
    Y_rest = Y_rest(idx);
    n_rest = length(X_rest);
    
    X_prep = arm_pref_m1_ic.prep;
    Y_prep = log10(max(modulation_m1_ic.prep,[],2));
    idx = X_prep<0 & ~isnan(X_prep) & Y_prep>-1;
    X_prep = X_prep(idx);
    Y_prep = Y_prep(idx);
    n_prep = length(X_prep);
    
    X_move = arm_pref_m1_ic.move;
    Y_move = log10(max(modulation_m1_ic.move,[],2));
    idx = X_move<0 & ~isnan(X_move) & Y_move>-1;
    X_move = X_move(idx);
    Y_move = Y_move(idx);
    n_move = length(X_move);
    
    X_rp = [X_rest; X_prep];
    Y_rp = [Y_rest; Y_prep];
    r_idx = randi(length(X_rp), n_rest,1);
    X_rp_r = X_rp(r_idx);
    Y_rp_r = Y_rp(r_idx);
    X_rp_p = X_rp(setdiff(1:length(X_rp),r_idx));
    Y_rp_p = Y_rp(setdiff(1:length(Y_rp),r_idx));
    b_r = polyfit(X_rp_r, Y_rp_r, 1);
    b_p = polyfit(X_rp_p, Y_rp_p, 1);
    slope_m1_ipsi_rp(perm,1) = b_p(1) - b_r(1);
    
    X_rm = [X_rest; X_move];
    Y_rm = [Y_rest; X_move];
    r_idx = randi(length(X_rm), n_rest,1);
    X_rm_r = X_rm(r_idx);
    Y_rm_r = Y_rm(r_idx);
    X_rm_m = X_rm(setdiff(1:length(X_rm),r_idx));
    Y_rm_m = Y_rm(setdiff(1:length(Y_rm),r_idx));
    b_r = polyfit(X_rm_r, Y_rm_r, 1);
    b_m = polyfit(X_rm_m, Y_rm_m, 1);
    slope_m1_ipsi_rm(perm,1) = b_m(1) - b_r(1);
    
    X_pm = [X_prep; X_move];
    Y_pm = [Y_prep; X_move];
    p_idx = randi(length(X_pm), n_prep,1);
    X_pm_p = X_pm(p_idx);
    Y_pm_p = Y_pm(p_idx);
    X_pm_m = X_pm(setdiff(1:length(X_pm),p_idx));
    Y_pm_m = Y_pm(setdiff(1:length(Y_pm),p_idx));
    b_p = polyfit(X_pm_p, Y_pm_p, 1);
    b_m = polyfit(X_pm_m, Y_pm_m, 1);
    slope_m1_ipsi_pm(perm,1) = b_m(1) - b_p(1);
  
    
    % m1, contra
    X_rest = arm_pref_m1_ic.rest;
    Y_rest = log10(max(modulation_m1_ic.rest,[],2));
    idx = X_rest>0 & ~isnan(X_rest) & Y_rest>-1;
    X_rest = X_rest(idx);
    Y_rest = Y_rest(idx);
    n_rest = length(X_rest);
    
    X_prep = arm_pref_m1_ic.prep;
    Y_prep = log10(max(modulation_m1_ic.prep,[],2));
    idx = X_prep>0 & ~isnan(X_prep) & Y_prep>-1;
    X_prep = X_prep(idx);
    Y_prep = Y_prep(idx);
    n_prep = length(X_prep);
    
    X_move = arm_pref_m1_ic.move;
    Y_move = log10(max(modulation_m1_ic.move,[],2));
    idx = X_move>0 & ~isnan(X_move) & Y_move>-1;
    X_move = X_move(idx);
    Y_move = Y_move(idx);
    n_move = length(X_move);
    
    X_rp = [X_rest; X_prep];
    Y_rp = [Y_rest; Y_prep];
    r_idx = randi(length(X_rp), n_rest,1);
    X_rp_r = X_rp(r_idx);
    Y_rp_r = Y_rp(r_idx);
    X_rp_p = X_rp(setdiff(1:length(X_rp),r_idx));
    Y_rp_p = Y_rp(setdiff(1:length(Y_rp),r_idx));
    b_r = polyfit(X_rp_r, Y_rp_r, 1);
    b_p = polyfit(X_rp_p, Y_rp_p, 1);
    slope_m1_contra_rp(perm,1) = b_p(1) - b_r(1);
    
    X_rm = [X_rest; X_move];
    Y_rm = [Y_rest; X_move];
    r_idx = randi(length(X_rm), n_rest,1);
    X_rm_r = X_rm(r_idx);
    Y_rm_r = Y_rm(r_idx);
    X_rm_m = X_rm(setdiff(1:length(X_rm),r_idx));
    Y_rm_m = Y_rm(setdiff(1:length(Y_rm),r_idx));
    b_r = polyfit(X_rm_r, Y_rm_r, 1);
    b_m = polyfit(X_rm_m, Y_rm_m, 1);
    slope_m1_contra_rm(perm,1) = b_m(1) - b_r(1);
    
    X_pm = [X_prep; X_move];
    Y_pm = [Y_prep; X_move];
    p_idx = randi(length(X_pm), n_prep,1);
    X_pm_p = X_pm(p_idx);
    Y_pm_p = Y_pm(p_idx);
    X_pm_m = X_pm(setdiff(1:length(X_pm),p_idx));
    Y_pm_m = Y_pm(setdiff(1:length(Y_pm),p_idx));
    b_p = polyfit(X_pm_p, Y_pm_p, 1);
    b_m = polyfit(X_pm_m, Y_pm_m, 1);
    slope_m1_contra_pm(perm,1) = b_m(1) - b_p(1);
    
end

% record observed differences in slope

X_rest = arm_pref_pmd_ic.rest;
Y_rest = log10(max(modulation_pmd_ic.rest,[],2));
idx = X_rest<0 & ~isnan(X_rest) & Y_rest>-1;
X_rest = X_rest(idx);
Y_rest = Y_rest(idx);
b_r_pmd_ipsi = polyfit(X_rest, Y_rest, 1);

X_prep = arm_pref_pmd_ic.prep;
Y_prep = log10(max(modulation_pmd_ic.prep,[],2));
idx = X_prep<0 & ~isnan(X_prep) & Y_prep>-1;
X_prep = X_prep(idx);
Y_prep = Y_prep(idx);
b_p_pmd_ipsi = polyfit(X_prep, Y_prep, 1);

X_move = arm_pref_pmd_ic.move;
Y_move = log10(max(modulation_pmd_ic.move,[],2));
idx = X_move<0 & ~isnan(X_move) & Y_move>-1;
X_move = X_move(idx);
Y_move = Y_move(idx);
b_m_pmd_ipsi = polyfit(X_move, Y_move, 1);

X_rest = arm_pref_pmd_ic.rest;
Y_rest = log10(max(modulation_pmd_ic.rest,[],2));
idx = X_rest>0 & ~isnan(X_rest) & Y_rest>-1;
X_rest = X_rest(idx);
Y_rest = Y_rest(idx);
b_r_pmd_contra = polyfit(X_rest, Y_rest, 1);

X_prep = arm_pref_pmd_ic.prep;
Y_prep = log10(max(modulation_pmd_ic.prep,[],2));
idx = X_prep>0 & ~isnan(X_prep) & Y_prep>-1;
X_prep = X_prep(idx);
Y_prep = Y_prep(idx);
b_p_pmd_contra = polyfit(X_prep, Y_prep, 1);

X_move = arm_pref_pmd_ic.move;
Y_move = log10(max(modulation_pmd_ic.move,[],2));
idx = X_move>0 & ~isnan(X_move) & Y_move>-1;
X_move = X_move(idx);
Y_move = Y_move(idx);
b_m_pmd_contra = polyfit(X_move, Y_move, 1);

X_rest = arm_pref_m1_ic.rest;
Y_rest = log10(max(modulation_m1_ic.rest,[],2));
idx = X_rest<0 & ~isnan(X_rest) & Y_rest>-1;
X_rest = X_rest(idx);
Y_rest = Y_rest(idx);
b_r_m1_ipsi = polyfit(X_rest, Y_rest, 1);

X_prep = arm_pref_m1_ic.prep;
Y_prep = log10(max(modulation_m1_ic.prep,[],2));
idx = X_prep<0 & ~isnan(X_prep) & Y_prep>-1;
X_prep = X_prep(idx);
Y_prep = Y_prep(idx);
b_p_m1_ipsi = polyfit(X_prep, Y_prep, 1);

X_move = arm_pref_m1_ic.move;
Y_move = log10(max(modulation_m1_ic.move,[],2));
idx = X_move<0 & ~isnan(X_move) & Y_move>-1;
X_move = X_move(idx);
Y_move = Y_move(idx);
b_m_m1_ipsi = polyfit(X_move, Y_move, 1);

X_rest = arm_pref_m1_ic.rest;
Y_rest = log10(max(modulation_m1_ic.rest,[],2));
idx = X_rest>0 & ~isnan(X_rest) & Y_rest>-1;
X_rest = X_rest(idx);
Y_rest = Y_rest(idx);
b_r_m1_contra = polyfit(X_rest, Y_rest, 1);

X_prep = arm_pref_m1_ic.prep;
Y_prep = log10(max(modulation_m1_ic.prep,[],2));
idx = X_prep>0 & ~isnan(X_prep) & Y_prep>-1;
X_prep = X_prep(idx);
Y_prep = Y_prep(idx);
b_p_m1_contra = polyfit(X_prep, Y_prep, 1);

X_move = arm_pref_m1_ic.move;
Y_move = log10(max(modulation_m1_ic.move,[],2));
idx = X_move>0 & ~isnan(X_move) & Y_move>-1;
X_move = X_move(idx);
Y_move = Y_move(idx);
b_m_m1_contra = polyfit(X_move, Y_move, 1);

% assign p-value based on number of permutations that have equal or more
% extreme differences in slope
p_pmd_ipsi_rp = sum(abs(slope_pmd_ipsi_rp)>...
    abs(b_r_pmd_ipsi(1)-b_p_pmd_ipsi(1)))/10000;
p_pmd_ipsi_rm = sum(abs(slope_pmd_ipsi_rm)>...
    abs(b_r_pmd_ipsi(1)-b_m_pmd_ipsi(1)))/10000;
p_pmd_ipsi_pm = sum(abs(slope_pmd_ipsi_pm)>...
    abs(b_p_pmd_ipsi(1)-b_m_pmd_ipsi(1)))/10000;

p_pmd_contra_rp = sum(abs(slope_pmd_contra_rp)>...
    abs(b_r_pmd_contra(1)-b_p_pmd_contra(1)))/10000;
p_pmd_contra_rm = sum(abs(slope_pmd_contra_rm)>...
    abs(b_r_pmd_contra(1)-b_m_pmd_contra(1)))/10000;
p_pmd_contra_pm = sum(abs(slope_pmd_contra_pm)>...
    abs(b_p_pmd_contra(1)-b_m_pmd_contra(1)))/10000;

p_m1_ipsi_rp = sum(abs(slope_m1_ipsi_rp)>...
    abs(b_r_m1_ipsi(1)-b_p_m1_ipsi(1)))/10000;
p_m1_ipsi_rm = sum(abs(slope_m1_ipsi_rm)>...
    abs(b_r_m1_ipsi(1)-b_m_m1_ipsi(1)))/10000;
p_m1_ipsi_pm = sum(abs(slope_m1_ipsi_pm)>...
    abs(b_p_m1_ipsi(1)-b_m_m1_ipsi(1)))/10000;

p_m1_contra_rp = sum(abs(slope_m1_contra_rp)>...
    abs(b_r_m1_contra(1)-b_p_m1_contra(1)))/10000;
p_m1_contra_rm = sum(abs(slope_m1_contra_rm)>...
    abs(b_r_m1_contra(1)-b_m_m1_contra(1)))/10000;
p_m1_contra_pm = sum(abs(slope_m1_contra_pm)>...
    abs(b_p_m1_contra(1)-b_m_m1_contra(1)))/10000;

% plot
figure('Name', 'Slope of arm dedication, MD relationship across epochs')
set(gcf, 'Position',  [200, 200, 350, 600])

subplot(2,1,1)
errorbar([1:3]-0.05, Y_pmd_ipsi, neg_pmd_ipsi, pos_pmd_ipsi, ...
    'o-', 'LineWidth',1.5)
hold on
errorbar([1:3]+0.05, Y_pmd_contra, neg_pmd_contra, pos_pmd_contra, ...
    'o-', 'LineWidth',1.5)
xlim([0.5, 3.5])
xticks([1,2,3])
xticklabels({'Rest','Prep','Move'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
yticks([0,1])
ylabel('Slope','fontsize',16)
ylim([-0.6,1.6])
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
legend({'Ipsi','Contra'})
title('PMd')

subplot(2,1,2)
errorbar([1:3]-0.05, Y_m1_ipsi, neg_m1_ipsi, pos_m1_ipsi, ...
    'o-', 'LineWidth',1.5)
hold on
errorbar([1:3]+0.05, Y_m1_contra, neg_m1_contra, pos_m1_contra, ...
    'o-', 'LineWidth',1.5)
xlim([0.5, 3.5])
xticks([1,2,3])
xticklabels({'Rest','Prep','Move'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
yticks([0,1])
ylabel('Slope','fontsize',16)
ylim([-0.6,1.6])
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
legend({'Ipsi','Contra'})
title('M1')


%% plot cumulative modulation accounted for, rest phase
% Come back to this to bootstrap confidence intervals

% sort by arm preference (using held out data for assignments)
[sorted_arm_pref, idx] = ...
    sort([arm_pref_pmd_ic.rest; arm_pref_m1_ic.rest]);
sorted_modulation = ...
    [modulation_pmd_ic.rest; modulation_m1_ic.rest];
sorted_modulation = sorted_modulation(idx,:);
idx = ~isnan(sorted_arm_pref);
sorted_arm_pref = sorted_arm_pref(idx);
sorted_modulation = sorted_modulation(idx,:);

boot_cumsum_ipsi = zeros(10000,length(sorted_arm_pref)+2);
boot_cumsum_contra = zeros(10000,length(sorted_arm_pref)+2);
parfor boot = 1:10000
    
    X = sorted_arm_pref;
    Y = sorted_modulation;
    [X, boot_idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(boot_idx,:);
    [X, sort_idx] = sort(X);
    Y = Y(sort_idx,:);
    X = [-1;X;1];
    Y = [cumsum(Y(:,1))/sum(Y(:,1)), cumsum(Y(:,2))/sum(Y(:,2))];
    Y = [0,0; Y; 1,1];
    
    [resampled_Y1, ~] = resample(Y(:,1),X);
    [resampled_Y2, ~] = resample(Y(:,2),X);
    
    boot_cumsum_ipsi(boot,:) = [resampled_Y1];
    boot_cumsum_contra(boot,:) = [resampled_Y2];
    
end
X = linspace(-1,1,size(boot_cumsum_ipsi,2));

% open figure window
figure('Name', ...
    'Ipsi unit-level modulation separation (Rest)')
set(gcf, 'Position',  [2000, 200, 600, 250])

% plot data
shadedErrorBar(X,boot_cumsum_ipsi,{@mean,@std},'lineprops','b')
hold on
shadedErrorBar(X,boot_cumsum_contra,{@mean,@std},'lineprops','r')

% format view
xlim([-1.05, 1.05])
xlabel('Arm preference','fontsize',16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
ylim([-0.025, 1.025])
yticks([0,0.25,0.50,0.75,1.00])
yticklabels({'0','','50','','100'}) 
ylabel({'Cumulative modulation', 'accounted for'},'fontsize',16)
box on


% open figure window
figure('Name', ...
    'Contra unit-level modulation separation (Rest)')
set(gcf, 'Position',  [2000, 200, 600, 250])

% plot data
shadedErrorBar(-X,1-boot_cumsum_ipsi,{@mean,@std},'lineprops','b')
hold on
shadedErrorBar(-X,1-boot_cumsum_contra,{@mean,@std},'lineprops','r')

% format view
xlim([-1.05, 1.05])
xticks([-1,0,1])
xticklabels({'1','','-1'})
xlabel('Arm preference','fontsize',16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
ylim([-0.025, 1.025])
yticks([0,0.25,0.50,0.75,1.00])
yticklabels({'0','','50','','100'}) 
ylabel({'% Modulation', 'accounted for'},'fontsize',16)
box on


%% plot VAF in each of the thirds of the arm preference distribution (rest)

% find indices for +/-0.3 and +/-0.4 arm pref index
ipsi_idx = find(X>-0.4,1);
contra_idx = find(X>0.4,1)-1;
share_idx = [find(X>-0.3,1), find(X>0.3,1)-1];

% collect data
ipsi_in_ipsi = sort(boot_cumsum_ipsi(:,ipsi_idx));
contra_in_ipsi = sort(boot_cumsum_contra(:,ipsi_idx));

contra_in_contra = sort(1-boot_cumsum_contra(:,contra_idx));
ipsi_in_contra = sort(1-boot_cumsum_ipsi(:,contra_idx));

ipsi_in_shared = ...
    sort(boot_cumsum_ipsi(:,share_idx(2)) - boot_cumsum_ipsi(:,share_idx(1)));
contra_in_shared = ...
    sort(boot_cumsum_contra(:,share_idx(2)) - boot_cumsum_contra(:,share_idx(1)));

ipsi_means(1,:) = ...
    [mean(ipsi_in_contra), mean(ipsi_in_shared), mean(ipsi_in_ipsi)];
contra_means(1,:) = ...
    [mean(contra_in_contra), mean(contra_in_shared), mean(contra_in_ipsi)];

ipsi_neg(1,:) = ipsi_means(1,:) - ...
    [quantile(ipsi_in_contra, 0.025), quantile(ipsi_in_shared, 0.025), quantile(ipsi_in_ipsi, 0.025)];
contra_neg(1,:) = contra_means(1,:) - ...
    [quantile(contra_in_contra, 0.025), quantile(contra_in_shared, 0.025), quantile(contra_in_ipsi, 0.025)];

ipsi_pos(1,:) = -ipsi_means(1,:) + ...
    [quantile(ipsi_in_contra, 0.975), quantile(ipsi_in_shared, 0.975), quantile(ipsi_in_ipsi, 0.975)];
contra_pos(1,:) = -contra_means(1,:) + ...
    [quantile(contra_in_contra, 0.975), quantile(contra_in_shared, 0.975), quantile(contra_in_ipsi, 0.975)];


% plot
figure('Name', 'Modulation captured in each regime of arm preferences (Rest)')
set(gcf, 'Position',  [200, 200, 350, 300])

errorbar([1:3]-0.05, ipsi_means(1,:), ipsi_neg(1,:), ipsi_pos(1,:), ...
    'o', 'LineWidth',1.5)
hold on
errorbar([1:3]+0.05, contra_means(1,:), contra_neg(1,:), contra_pos(1,:), ...
    'o', 'LineWidth',1.5)
xlim([0.5, 3.5])
xticks([1,2,3])
xticklabels({'Contra','Neutral','Ipsi'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
yticks([0, 0.25, 0.5, 0.75, 1])
yticklabels({'0','','','','1'})
ylabel({'% Modulation', 'accounted for'},'fontsize',16)
ylim([0,1])
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
legend({'Ipsi','Contra'})
title('Rest')


%% plot cumulative modulation accounted for, instruct phase
% Come back to this to bootstrap confidence intervals

% sort by arm preference (using held out data for assignments)
[sorted_arm_pref, idx] = ...
    sort([arm_pref_pmd_ic.prep; arm_pref_m1_ic.prep]);
sorted_modulation = ...
    [modulation_pmd_ic.prep; modulation_m1_ic.prep];
sorted_modulation = sorted_modulation(idx,:);
idx = ~isnan(sorted_arm_pref);
sorted_arm_pref = sorted_arm_pref(idx);
sorted_modulation = sorted_modulation(idx,:);

boot_cumsum_ipsi = zeros(10000,length(sorted_arm_pref)+2);
boot_cumsum_contra = zeros(10000,length(sorted_arm_pref)+2);
parfor boot = 1:10000
    
    X = sorted_arm_pref;
    Y = sorted_modulation;
    [X, boot_idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(boot_idx,:);
    [X, sort_idx] = sort(X);
    Y = Y(sort_idx,:);
    X = [-1;X;1];
    Y = [cumsum(Y(:,1))/sum(Y(:,1)), cumsum(Y(:,2))/sum(Y(:,2))];
    Y = [0,0; Y; 1,1];
    
    [resampled_Y1, ~] = resample(Y(:,1),X);
    [resampled_Y2, ~] = resample(Y(:,2),X);
    
    boot_cumsum_ipsi(boot,:) = [resampled_Y1];
    boot_cumsum_contra(boot,:) = [resampled_Y2];
    
end
X = linspace(-1,1,size(boot_cumsum_ipsi,2));

% open figure window
figure('Name', ...
    'Ipsi unit-level modulation separation (Instruct)')
set(gcf, 'Position',  [2000, 200, 600, 250])

% plot data
shadedErrorBar(X,boot_cumsum_ipsi,{@mean,@std},'lineprops','b')
hold on
shadedErrorBar(X,boot_cumsum_contra,{@mean,@std},'lineprops','r')

% format view
xlim([-1.05, 1.05])
xlabel('Arm preference','fontsize',16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
ylim([-0.025, 1.025])
yticks([0,0.25,0.50,0.75,1.00])
yticklabels({'0','','50','','100'}) 
ylabel({'Cumulative modulation', 'accounted for'},'fontsize',16)
box on


% open figure window
figure('Name', ...
    'Contra unit-level modulation separation (Instruct)')
set(gcf, 'Position',  [2000, 200, 600, 250])

% plot data
shadedErrorBar(-X,1-boot_cumsum_ipsi,{@mean,@std},'lineprops','b')
hold on
shadedErrorBar(-X,1-boot_cumsum_contra,{@mean,@std},'lineprops','r')

% format view
xlim([-1.05, 1.05])
xticks([-1,0,1])
xticklabels({'1','','-1'})
xlabel('Arm preference','fontsize',16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
ylim([-0.025, 1.025])
yticks([0,0.25,0.50,0.75,1.00])
yticklabels({'0','','50','','100'}) 
ylabel({'% Modulation', 'accounted for'},'fontsize',16)
box on


%% plot VAF in each of the thirds of the arm preference distribution (Instruct)

% find indices for +/-0.3 and +/-0.4 arm pref index
ipsi_idx = find(X>-0.4,1);
contra_idx = find(X>0.4,1)-1;
share_idx = [find(X>-0.3,1), find(X>0.3,1)-1];

% collect data
ipsi_in_ipsi = sort(boot_cumsum_ipsi(:,ipsi_idx));
contra_in_ipsi = sort(boot_cumsum_contra(:,ipsi_idx));

contra_in_contra = sort(1-boot_cumsum_contra(:,contra_idx));
ipsi_in_contra = sort(1-boot_cumsum_ipsi(:,contra_idx));

ipsi_in_shared = ...
    sort(boot_cumsum_ipsi(:,share_idx(2)) - boot_cumsum_ipsi(:,share_idx(1)));
contra_in_shared = ...
    sort(boot_cumsum_contra(:,share_idx(2)) - boot_cumsum_contra(:,share_idx(1)));

ipsi_means(2,:) = ...
    [mean(ipsi_in_contra), mean(ipsi_in_shared), mean(ipsi_in_ipsi)];
contra_means(2,:) = ...
    [mean(contra_in_contra), mean(contra_in_shared), mean(contra_in_ipsi)];

ipsi_neg(2,:) = ipsi_means(2,:) - ...
    [quantile(ipsi_in_contra, 0.025), quantile(ipsi_in_shared, 0.025), quantile(ipsi_in_ipsi, 0.025)];
contra_neg(2,:) = contra_means(2,:) - ...
    [quantile(contra_in_contra, 0.025), quantile(contra_in_shared, 0.025), quantile(contra_in_ipsi, 0.025)];

ipsi_pos(2,:) = -ipsi_means(2,:) + ...
    [quantile(ipsi_in_contra, 0.975), quantile(ipsi_in_shared, 0.975), quantile(ipsi_in_ipsi, 0.975)];
contra_pos(2,:) = -contra_means(2,:) + ...
    [quantile(contra_in_contra, 0.975), quantile(contra_in_shared, 0.975), quantile(contra_in_ipsi, 0.975)];


% plot
figure('Name', 'Modulation captured in each regime of arm preferences (Instruct)')
set(gcf, 'Position',  [200, 200, 350, 300])

errorbar([1:3]-0.05, ipsi_means(2,:), ipsi_neg(2,:), ipsi_pos(2,:), ...
    'o', 'LineWidth',1.5)
hold on
errorbar([1:3]+0.05, contra_means(2,:), contra_neg(2,:), contra_pos(2,:), ...
    'o', 'LineWidth',1.5)
xlim([0.5, 3.5])
xticks([1,2,3])
xticklabels({'Contra','Neutral','Ipsi'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
yticks([0, 0.25, 0.5, 0.75, 1])
yticklabels({'0','','','','1'})
ylabel({'% Modulation', 'accounted for'},'fontsize',16)
ylim([0,1])
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
legend({'Ipsi','Contra'})
title('Instruct')


%% plot cumulative modulation accounted for, move phase
% Come back to this to bootstrap confidence intervals

% sort by arm preference (using held out data for assignments)
[sorted_arm_pref, idx] = ...
    sort([arm_pref_pmd_ic.move; arm_pref_m1_ic.move]);
sorted_modulation = ...
    [modulation_pmd_ic.move; modulation_m1_ic.move];
sorted_modulation = sorted_modulation(idx,:);
idx = ~isnan(sorted_arm_pref);
sorted_arm_pref = sorted_arm_pref(idx);
sorted_modulation = sorted_modulation(idx,:);

boot_cumsum_ipsi = zeros(10000,length(sorted_arm_pref)+2);
boot_cumsum_contra = zeros(10000,length(sorted_arm_pref)+2);
parfor boot = 1:10000
    
    X = sorted_arm_pref;
    Y = sorted_modulation;
    [X, boot_idx] = datasample(X, length(X), 'Replace',true);
    Y = Y(boot_idx,:);
    [X, sort_idx] = sort(X);
    Y = Y(sort_idx,:);
    X = [-1;X;1];
    Y = [cumsum(Y(:,1))/sum(Y(:,1)), cumsum(Y(:,2))/sum(Y(:,2))];
    Y = [0,0; Y; 1,1];
    
    [resampled_Y1, ~] = resample(Y(:,1),X);
    [resampled_Y2, ~] = resample(Y(:,2),X);
    
    boot_cumsum_ipsi(boot,:) = [resampled_Y1];
    boot_cumsum_contra(boot,:) = [resampled_Y2];
    
end
X = linspace(-1,1,size(boot_cumsum_ipsi,2));

% open figure window
figure('Name', ...
    'Ipsi unit level modulation separation (Move)')
set(gcf, 'Position',  [2000, 200, 600, 250])

% plot data
shadedErrorBar(X,boot_cumsum_ipsi,{@mean,@std},'lineprops','b')
hold on
shadedErrorBar(X,boot_cumsum_contra,{@mean,@std},'lineprops','r')

% format view
xlim([-1.05, 1.05])
xlabel('Arm preference','fontsize',16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
ylim([-0.025, 1.025])
yticks([0,0.25,0.50,0.75,1.00])
yticklabels({'0','','50','','100'}) 
ylabel({'Cumulative modulation', 'accounted for'},'fontsize',16)
box on


% open figure window
figure('Name', ...
    'Contra unit level modulation separation (Move)')
set(gcf, 'Position',  [2000, 200, 600, 250])

% plot data
shadedErrorBar(-X,1-boot_cumsum_ipsi,{@mean,@std},'lineprops','b')
hold on
shadedErrorBar(-X,1-boot_cumsum_contra,{@mean,@std},'lineprops','r')

% format view
xlim([-1.05, 1.05])
xticks([-1,0,1])
xticklabels({'1','','-1'})
xlabel('Arm preference','fontsize',16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
ylim([-0.025, 1.025])
yticks([0,0.25,0.50,0.75,1.00])
yticklabels({'0','','50','','100'}) 
ylabel({'% Modulation', 'accounted for'},'fontsize',16)
box on


%% plot VAF in each of the thirds of the arm preference distribution (Move)

% find indices for +/-0.3 and +/-0.4 arm pref index
ipsi_idx = find(X>-0.4,1);
contra_idx = find(X>0.4,1)-1;
share_idx = [find(X>-0.3,1), find(X>0.3,1)-1];

% collect data
ipsi_in_ipsi = sort(boot_cumsum_ipsi(:,ipsi_idx));
contra_in_ipsi = sort(boot_cumsum_contra(:,ipsi_idx));

contra_in_contra = sort(1-boot_cumsum_contra(:,contra_idx));
ipsi_in_contra = sort(1-boot_cumsum_ipsi(:,contra_idx));

ipsi_in_shared = ...
    sort(boot_cumsum_ipsi(:,share_idx(2)) - boot_cumsum_ipsi(:,share_idx(1)));
contra_in_shared = ...
    sort(boot_cumsum_contra(:,share_idx(2)) - boot_cumsum_contra(:,share_idx(1)));

ipsi_means(3,:) = ...
    [mean(ipsi_in_contra), mean(ipsi_in_shared), mean(ipsi_in_ipsi)];
contra_means(3,:) = ...
    [mean(contra_in_contra), mean(contra_in_shared), mean(contra_in_ipsi)];

ipsi_neg(3,:) = ipsi_means(3,:) - ...
    [quantile(ipsi_in_contra, 0.025), quantile(ipsi_in_shared, 0.025), quantile(ipsi_in_ipsi, 0.025)];
contra_neg(3,:) = contra_means(3,:) - ...
    [quantile(contra_in_contra, 0.025), quantile(contra_in_shared, 0.025), quantile(contra_in_ipsi, 0.025)];

ipsi_pos(3,:) = -ipsi_means(3,:) + ...
    [quantile(ipsi_in_contra, 0.975), quantile(ipsi_in_shared, 0.975), quantile(ipsi_in_ipsi, 0.975)];
contra_pos(3,:) = -contra_means(3,:) + ...
    [quantile(contra_in_contra, 0.975), quantile(contra_in_shared, 0.975), quantile(contra_in_ipsi, 0.975)];


% plot
figure('Name', 'Modulation captured in each regime of arm preferences (Move)')
set(gcf, 'Position',  [200, 200, 350, 300])

errorbar([1:3]-0.05, ipsi_means(3,:), ipsi_neg(3,:), ipsi_pos(3,:), ...
    'o', 'LineWidth',1.5)
hold on
errorbar([1:3]+0.05, contra_means(3,:), contra_neg(3,:), contra_pos(3,:), ...
    'o', 'LineWidth',1.5)
xlim([0.5, 3.5])
xticks([1,2,3])
xticklabels({'Contra','Neutral','Ipsi'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
yticks([0, 0.25, 0.5, 0.75, 1])
yticklabels({'0','','','','1'})
ylabel({'% Modulation', 'accounted for'},'fontsize',16)
ylim([0,1])
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
legend({'Ipsi','Contra'})
title('Move')


%% plot all phases together
figure('Name', 'Modulation captured in each regime of arm preferences (All phases)')
set(gcf, 'Position',  [200, 200, 500, 250])

errorbar([1:3]+0.1, ipsi_means(:,1), ipsi_neg(:,1), ipsi_pos(:,1), ...
    'bo-', 'LineWidth',1.5)
hold on
errorbar([1:3]-0.1, contra_means(:,1), contra_neg(:,1), contra_pos(:,1), ...
    'ro-', 'LineWidth',1.5)

errorbar([5:7]+0.1, ipsi_means(:,2), ipsi_neg(:,2), ipsi_pos(:,2), ...
    'bo-', 'LineWidth',1.5)
hold on
errorbar([5:7]-0.1, contra_means(:,2), contra_neg(:,2), contra_pos(:,2), ...
    'ro-', 'LineWidth',1.5)

errorbar([9:11]+0.1, ipsi_means(:,3), ipsi_neg(:,3), ipsi_pos(:,3), ...
    'bo-', 'LineWidth',1.5)
hold on
errorbar([9:11]-0.1, contra_means(:,3), contra_neg(:,3), contra_pos(:,3), ...
    'ro-', 'LineWidth',1.5)

xlim([0, 12])
xticks([2,6,10])
xticklabels({'Contra','Neutral','Ipsi'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)
yticks([0, 0.25, 0.5, 0.75, 1])
yticklabels({'0','','','','1'})
ylabel({'% Modulation', 'accounted for'},'fontsize',16)
ylim([0,1])
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
legend({'Ipsi','Contra'})
title('MD dist across phases')
