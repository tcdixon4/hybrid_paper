function [speed_profiles] = plot_speed_profiles(trials, target)

%
% Plots example speed profiles for reaches to each target.
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


%% log the trajectories

speed_profiles = struct(...
    'left',struct('reach_hand',[],'stat_hand',[]),...
    'right',struct('reach_hand',[],'stat_hand',[]));
speed_profiles(6).left.reach_hand = cell(1,1);
speed_profiles(6).left.stat_hand = cell(1,1);
speed_profiles(6).right.reach_hand = cell(1,1);
speed_profiles(6).right.stat_hand = cell(1,1);

left_trials = trials(([trials.reach_hand]==1) & ([trials.config]==1));
right_trials = trials(([trials.reach_hand]==2) & ([trials.config]==-1));

for t = 6:-1:1
    % left hand reaches
    trial = find([left_trials.target_id]==t);
    trials_per_target = length(trial);
    for k = trials_per_target:-1:1
        start = find(left_trials(trial(k)).l_hand.time > ...
            left_trials(trial(k)).reach_onset, 1) - 96;
        stop = find(left_trials(trial(k)).l_hand.time > ...
            left_trials(trial(k)).enter_targ+0.2, 1);
        speed_profiles(t).left.reach_hand{k} = ...
            left_trials(trial(k)).l_hand.speed(start:stop,:);
        speed_profiles(t).left.stat_hand{k} = ...
            left_trials(trial(k)).r_hand.speed(start:stop,:);
    end
    % right hand reaches
    trial = find([right_trials.target_id]==t);
    trials_per_target = length(trial);
    for k = trials_per_target:-1:1
        start = find(right_trials(trial(k)).l_hand.time > ...
            right_trials(trial(k)).reach_onset, 1) - 96;
        stop = find(right_trials(trial(k)).l_hand.time > ...
            right_trials(trial(k)).enter_targ+0.2, 1);
        speed_profiles(t).right.reach_hand{k} = ...
            right_trials(trial(k)).r_hand.speed(start:stop,:);
        speed_profiles(t).right.stat_hand{k} = ...
            right_trials(trial(k)).l_hand.speed(start:stop,:);
    end
end

%% Calculate the mean and SD of the velocity profiles, then plot
 
%%%% Calculate mean and SD %%%%

% left hand reaches
% log number of trials to the current target/hand combo
num_trials = length(speed_profiles(target).left.reach_hand);
% log the maximum trial length and pad a matrix for logging each trial
num_samples = 0;
for k = 1:length(speed_profiles(target).left.reach_hand)
    num_samples = max(num_samples, ...
        length(speed_profiles(target).left.reach_hand{k}));
end
velocity_mat = nan(num_trials, num_samples);
% log each trial
for trial = num_trials:-1:1
    trial_length = length(speed_profiles(target).left.reach_hand{trial});
    velocity_mat(trial, 1:trial_length) = ...
        speed_profiles(target).left.reach_hand{trial};
end
mu_reach = nanmean(velocity_mat);
sd_reach = nanstd(velocity_mat);

velocity_mat = nan(num_trials, num_samples);
% log each trial
for trial = num_trials:-1:1
    trial_length = length(speed_profiles(target).left.stat_hand{trial});
    velocity_mat(trial, 1:trial_length) = ...
        speed_profiles(target).left.stat_hand{trial};
end
mu_stat = nanmean(velocity_mat);
sd_stat = nanstd(velocity_mat);
t = ((0:(num_samples-1)) - 96)/.24;

% plot
figure('Position', [100,100,500,500])
shadedErrorBar(t,mu_reach,sd_reach)
hold on
shadedErrorBar(t,mu_stat,sd_stat)

xlim([-200, 400])
ylim([-0.1, 1.2])

xlabel('Time (ms)')
ylabel('Speed (m/s)')
title(['Left Hand, Target: ', num2str(target)])

% right hand reaches
% log number of trials to the current target/hand combo
num_trials = length(speed_profiles(target).right.reach_hand);
% log the maximum trial length and pad a matrix for logging each trial
num_samples = 0;
for k = 1:length(speed_profiles(target).right.reach_hand)
    num_samples = max(num_samples, ...
        length(speed_profiles(target).right.reach_hand{k}));
end
velocity_mat = nan(num_trials, num_samples);
% log each trial
for trial = 1:num_trials
    trial_length = length(speed_profiles(target).right.reach_hand{trial});
    velocity_mat(trial, 1:trial_length) = ...
        speed_profiles(target).right.reach_hand{trial};
end
mu_reach = nanmean(velocity_mat);
sd_reach = nanstd(velocity_mat);

velocity_mat = nan(num_trials, num_samples);
% log each trial
for trial = 1:num_trials
    trial_length = length(speed_profiles(target).right.stat_hand{trial});
    velocity_mat(trial, 1:trial_length) = ...
        speed_profiles(target).right.stat_hand{trial};
end
mu_stat = nanmean(velocity_mat);
sd_stat = nanstd(velocity_mat);
t = ((0:(num_samples-1)) - 96)/.24;

% plot
figure('Position', [100,100,500,500])
shadedErrorBar(t,mu_reach,sd_reach)
hold on
shadedErrorBar(t,mu_stat,sd_stat)

xlim([-200, 400])
ylim([-0.1, 1.2])

xlabel('Time (ms)')
ylabel('Speed (m/s)')
title(['Right Hand, Target: ', num2str(target)])


