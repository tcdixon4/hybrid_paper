function [speed_profiles] = plot_speed_profiles(trials, target)

%
% Plots example speed profiles for reaches to the specified target. The
% output contains data for all targets, but the code currently only plots
% one. This is a bit weird, and should be updated in the future to plot all
% targets.
%
%
% INPUTS: 
%
% trials - trial-separated data struct (struct: 1 x num_trials)
%
% target - selected target to plot, 1:6
%
% OUTPUTS:
%
% speed_profiles - nested data struct containing speed profiles, each hand 
%                  (struct: 1x6, fields 'left' and 'right', row=target. 
%                   nested struct 1x1: fields 'reach_hand' and 'stat_hand'
%                   nested cell array 1x trials_per_target: row=trial)
%                   nested matrix num_samples x 1: speed profile
% 


%% prepare data structures for both outward and return portions of reaches

% preallocate speed profiles for outward portion of the reach
speed_profiles = struct(...
    'left',struct('reach_hand',[],'stat_hand',[]),...
    'right',struct('reach_hand',[],'stat_hand',[]));
speed_profiles(6).left.reach_hand = cell(1,1);
speed_profiles(6).left.stat_hand = cell(1,1);
speed_profiles(6).right.reach_hand = cell(1,1);
speed_profiles(6).right.stat_hand = cell(1,1);

% preallocate speed profiles for return portion of the reach
return_speed_profiles = struct(...
    'left',struct('reach_hand',[],'stat_hand',[]),...
    'right',struct('reach_hand',[],'stat_hand',[]));
return_speed_profiles(6).left.reach_hand = cell(1,1);
return_speed_profiles(6).left.stat_hand = cell(1,1);
return_speed_profiles(6).right.reach_hand = cell(1,1);
return_speed_profiles(6).right.stat_hand = cell(1,1);

% separate into left and right trials
left_trials = trials(([trials.reach_hand]==1) & ([trials.config]==1));
right_trials = trials(([trials.reach_hand]==2) & ([trials.config]==-1));

clearvars trials 


%% log the trajectories, separated by target

for t = 6:-1:1
    % left hand reaches
    target_subset = find([left_trials.target_id]==t);
    trials_per_target = length(target_subset);
    for k = trials_per_target:-1:1
        % outward portion
        start = find(left_trials(target_subset(k)).l_hand.time > ...
            left_trials(target_subset(k)).reach_onset, 1) - 96;
        stop = find(left_trials(target_subset(k)).l_hand.time > ...
            left_trials(target_subset(k)).enter_targ+0.2, 1);
        speed_profiles(t).left.reach_hand{k} = ...
            left_trials(target_subset(k)).l_hand.speed(start:stop,:);
        speed_profiles(t).left.stat_hand{k} = ...
            left_trials(target_subset(k)).r_hand.speed(start:stop,:);
        % return portion
        start = find(left_trials(target_subset(k)).l_hand.time >= ...
            left_trials(target_subset(k)).return_peak, 1) - 72;
        stop = start + 144;
        return_speed_profiles(t).left.reach_hand{k} = ...
            left_trials(target_subset(k)).l_hand.speed(start:stop,:);
        return_speed_profiles(t).left.stat_hand{k} = ...
            left_trials(target_subset(k)).r_hand.speed(start:stop,:);
    end
    % right hand reaches
    target_subset = find([right_trials.target_id]==t);
    trials_per_target = length(target_subset);
    for k = trials_per_target:-1:1
        % outward portion
        start = find(right_trials(target_subset(k)).l_hand.time > ...
            right_trials(target_subset(k)).reach_onset, 1) - 96;
        stop = find(right_trials(target_subset(k)).l_hand.time > ...
            right_trials(target_subset(k)).enter_targ+0.2, 1);
        speed_profiles(t).right.reach_hand{k} = ...
            right_trials(target_subset(k)).r_hand.speed(start:stop,:);
        speed_profiles(t).right.stat_hand{k} = ...
            right_trials(target_subset(k)).l_hand.speed(start:stop,:);
        % return portion
        start = find(right_trials(target_subset(k)).l_hand.time >= ...
            right_trials(target_subset(k)).return_peak, 1) - 72;
        stop = start+144;
        return_speed_profiles(t).right.reach_hand{k} = ...
            right_trials(target_subset(k)).r_hand.speed(start:stop,:);
        return_speed_profiles(t).right.stat_hand{k} = ...
            right_trials(target_subset(k)).l_hand.speed(start:stop,:);
    end
end


%% LEFT HAND OUTWARD plotting
% Calculate the mean and SD of the velocity profiles, then plot

% log number of trials to the current target/hand combo
num_trials = length(speed_profiles(target).left.reach_hand);
% log the maximum trial length and pad a matrix for logging each trial
num_samples = 0;
for k = 1:num_trials
    num_samples = max(num_samples, ...
        length(speed_profiles(target).left.reach_hand{k}));
end
% Fill the data matrices for both the reach hand and stationary hand
reach_speed_mat = nan(num_trials, num_samples);
stat_speed_mat = nan(num_trials, num_samples);
for trial = num_trials:-1:1
    trial_length = length(speed_profiles(target).left.reach_hand{trial});
    reach_speed_mat(trial, 1:trial_length) = ...
        speed_profiles(target).left.reach_hand{trial};
    stat_speed_mat(trial, 1:trial_length) = ...
        speed_profiles(target).left.stat_hand{trial};
end
mu_reach = nanmean(reach_speed_mat);
sd_reach = nanstd(reach_speed_mat);
mu_stat = nanmean(stat_speed_mat);
sd_stat = nanstd(stat_speed_mat);
t = ((0:(num_samples-1)) - 96)/.24;

% plot
figure('Position', [100,100,500,500])
shadedErrorBar(t,mu_reach,sd_reach)
hold on
shadedErrorBar(t,mu_stat,sd_stat)
plot([0,0],[0,1], 'r')

xlim([-200, 400])
ylim([-0.1, 1.2])
yticks([0,0.5,1])

xlabel('Time (ms)')
ylabel('Speed (m/s)')
title(['Outward, Left Hand, Target: ', num2str(target)])


%% RIGHT HAND OUTWARD plotting
% Calculate the mean and SD of the velocity profiles, then plot

% log number of trials to the current target/hand combo
num_trials = length(speed_profiles(target).right.reach_hand);
% log the maximum trial length and pad a matrix for logging each trial
num_samples = 0;
for k = 1:num_trials
    num_samples = max(num_samples, ...
        length(speed_profiles(target).right.reach_hand{k}));
end
% Fill the data matrices for both the reach hand and stationary hand
reach_speed_mat = nan(num_trials, num_samples);
stat_speed_mat = nan(num_trials, num_samples);
for trial = num_trials:-1:1
    trial_length = length(speed_profiles(target).right.reach_hand{trial});
    reach_speed_mat(trial, 1:trial_length) = ...
        speed_profiles(target).right.reach_hand{trial};
    stat_speed_mat(trial, 1:trial_length) = ...
        speed_profiles(target).right.stat_hand{trial};
end
mu_reach = nanmean(reach_speed_mat);
sd_reach = nanstd(reach_speed_mat);
mu_stat = nanmean(stat_speed_mat);
sd_stat = nanstd(stat_speed_mat);
t = ((0:(num_samples-1)) - 96)/.24;

% plot
figure('Position', [100,100,500,500])
shadedErrorBar(t,mu_reach,sd_reach)
hold on
shadedErrorBar(t,mu_stat,sd_stat)
plot([0,0],[0,1], 'r')

xlim([-200, 400])
ylim([-0.1, 1.2])
yticks([0,0.5,1])

xlabel('Time (ms)')
ylabel('Speed (m/s)')
title(['Outward, Right Hand, Target: ', num2str(target)])


% %% LEFT HAND RETURN plotting 
% % Calculate the mean and SD of the velocity profiles, then plot
% 
% % convert to a matrix. Logging it as a cell array to begin with is a bit 
% % redundant, but it was important for the variable length vectors in the
% % outward portion. This should be fixed later.
% reach_speed_mat = return_speed_profiles(target).left.reach_hand;
% reach_speed_mat = ...
%     reach_speed_mat(~cellfun('isempty',reach_speed_mat));
% reach_speed_mat = cell2mat(reach_speed_mat)';
% stat_speed_mat = return_speed_profiles(target).left.stat_hand;
% stat_speed_mat = ...
%     stat_speed_mat(~cellfun('isempty',stat_speed_mat));
% stat_speed_mat = cell2mat(stat_speed_mat)';
% % remove rows(trials) where the stationary hand exceeded 9cm/s during the
% % first 300ms of the return
% max_stat_speed = max(stat_speed_mat(:,49:121), [], 2);
% idx = max_stat_speed < 0.09;
% stat_speed_mat = stat_speed_mat(idx,:);
% reach_speed_mat = reach_speed_mat(idx,:);
% % compute mean and SD
% mu_reach = mean(reach_speed_mat);
% sd_reach = std(reach_speed_mat);
% mu_stat = mean(stat_speed_mat);
% sd_stat = std(stat_speed_mat);
% num_samples = length(mu_reach);
% t = ((0:(num_samples-1)) - 72)/.24;
% 
% % plot
% figure('Position', [100,100,500,500])
% shadedErrorBar(t,mu_reach,sd_reach)
% hold on
% shadedErrorBar(t,mu_stat,sd_stat)
% plot([0,0],[0,1], 'r')
% 
% xlim([-300, 300])
% ylim([0, 1])
% yticks([0,0.5,1])
% 
% xlabel('Time (ms)')
% ylabel('Speed (m/s)')
% annotation('textbox', [0.5, 0.8, 0.1, 0.1], 'String', ...
%     ['Accepted trials: ', num2str(sum(idx)),'/', num2str(length(idx))])
% title(['Return, Left Hand, Target: ', num2str(target)])
% 
% 
% %% RIGHT HAND RETURN plotting 
% % Calculate the mean and SD of the velocity profiles, then plot
% 
% % convert to a matrix. Logging it as a cell array to begin with is a bit 
% % redundant, but it was important for the variable length vectors in the
% % outward portion. This should be fixed later.
% reach_speed_mat = return_speed_profiles(target).right.reach_hand;
% reach_speed_mat = ...
%     reach_speed_mat(~cellfun('isempty',reach_speed_mat));
% reach_speed_mat = cell2mat(reach_speed_mat)';
% stat_speed_mat = return_speed_profiles(target).right.stat_hand;
% stat_speed_mat = ...
%     stat_speed_mat(~cellfun('isempty',stat_speed_mat));
% stat_speed_mat = cell2mat(stat_speed_mat)';
% % remove rows(trials) where the stationary hand exceeded 9cm/s during the
% % first 300ms of the return
% max_stat_speed = max(stat_speed_mat(:,49:121), [], 2);
% idx = max_stat_speed < 0.09;
% stat_speed_mat = stat_speed_mat(idx,:);
% reach_speed_mat = reach_speed_mat(idx,:);
% % compute mean and SD
% sd_reach = std(reach_speed_mat);
% mu_stat = mean(stat_speed_mat);
% sd_stat = std(stat_speed_mat);
% num_samples = length(mu_reach);
% t = ((0:(num_samples-1)) - 72)/.24;
% 
% % plot
% figure('Position', [100,100,500,500])
% shadedErrorBar(t,mu_reach,sd_reach)
% hold on
% shadedErrorBar(t,mu_stat,sd_stat)
% plot([0,0],[0,1], 'r')
% 
% xlim([-300, 300])
% ylim([0, 1])
% yticks([0,0.5,1])
% 
% xlabel('Time (ms)')
% ylabel('Speed (m/s)')
% annotation('textbox', [0.5, 0.8, 0.1, 0.1], 'String', ...
%     ['Accepted trials: ', num2str(sum(idx)),'/', num2str(length(idx))])
% title(['Return, Right Hand, Target: ', num2str(target)])
% 
% 
% 
