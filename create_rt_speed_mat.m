function rt_speed_mat = create_rt_speed_mat(trials)

%
% Creates a simple matrix containing behavioral data.
%
%
% INPUTS: 
%
% trials - trial-separated data struct containing behavioral data
%          (struct: 1 x num_trials)
%
% OUTPUTS:
%
% rt_speed_mat - matrix summarizing features of behavior. see next section 
%                for column descriptions
%                (matrix: num_trials x 12)
% 


%% rt_speed_mat column descriptions

% col 1: Reach hand -              Left=0, Right=1
% col 2: Target -                  1:6
% col 3: RT -                      continuous variable
% col 4: Speed,reach hand,rest -   continuous variable, mean over 300ms
% col 5: Speed,stat hand,rest -    continuous variable, mean over 300ms
% col 6: Speed,reach hand,instr -  continuous variable, mean over 300ms
% col 7: Speed,stat hand,instr -   continuous variable, mean over 300ms
% col 8: Speed,reach hand,move -   continuous variable, mean over 300ms
% col 9: Speed,stat hand,move -    continuous variable, mean over 300ms
% col 10: Reach duration -         continuous variable, move->target entry
% col 11: Reach duration -         continuous variable, move->stop
% col 12: Minimum movement speed - continuous variable


%% Iterate over trials and log each variable
for trial = length(trials):-1:1
    
    %% find kinematic timestamps for Rest, Prep, and Move phases   
    targ_onset = find(trials(trial).l_hand.time >...
        trials(trial).targ_onset,1)-1;
    rest_idx = (targ_onset-71):targ_onset; % -300ms:targ_onset
    prep_idx = (targ_onset+48):(targ_onset+119); % targ_onset+200ms:+500ms
    
    reach_onset = find(trials(trial).l_hand.time >...
        trials(trial).reach_onset,1);
    move_idx = reach_onset:(reach_onset+71); % reach_onset:+300ms
    
    %% Log each variable
    % Reach hand
    rt_speed_mat(trial,1) = trials(trial).reach_hand-1;
    % Target
    rt_speed_mat(trial,2) = trials(trial).target_id;
    % RT
    rt_speed_mat(trial,3) = trials(trial).rt;
    % Mean speeds
    % also log depth velocity traces to determine reach end for duration
    if trials(trial).reach_hand==1 % left hand reaching        
        rt_speed_mat(trial,4) = ...
            mean(trials(trial).l_hand.speed(rest_idx));
        rt_speed_mat(trial,5) = ...
            mean(trials(trial).r_hand.speed(rest_idx));
        rt_speed_mat(trial,6) = ...
            mean(trials(trial).l_hand.speed(prep_idx));
        rt_speed_mat(trial,7) = ...
            mean(trials(trial).r_hand.speed(prep_idx));
        rt_speed_mat(trial,8) = ...
            mean(trials(trial).l_hand.speed(move_idx));
        rt_speed_mat(trial,9) = ...
            mean(trials(trial).r_hand.speed(move_idx));
        
        depth_vel = trials(trial).l_hand.velocity(...
            (reach_onset+24):(reach_onset+264),2); % reach_onset+100ms:+1.25s
        speed = trials(trial).l_hand.speed(...
            (reach_onset+24):(reach_onset+264)); % reach_onset+100ms:+1.25s
        
    else % right hand reaching
        rt_speed_mat(trial,4) = ...
            mean(trials(trial).r_hand.speed(rest_idx));
        rt_speed_mat(trial,5) = ...
            mean(trials(trial).l_hand.speed(rest_idx));
        rt_speed_mat(trial,6) = ...
            mean(trials(trial).r_hand.speed(prep_idx));
        rt_speed_mat(trial,7) = ...
            mean(trials(trial).l_hand.speed(prep_idx));
        rt_speed_mat(trial,8) = ...
            mean(trials(trial).r_hand.speed(move_idx));
        rt_speed_mat(trial,9) = ...
            mean(trials(trial).l_hand.speed(move_idx));
        
        depth_vel = trials(trial).r_hand.velocity(...
            (reach_onset+24):(reach_onset+264),2); % reach_onset+100ms:+1.25s
        speed = trials(trial).r_hand.speed(...
            (reach_onset+24):(reach_onset+264)); % reach_onset+100ms:+1.25s
        
    end   
    % Reach duration
    rt_speed_mat(trial,10) = trials(trial).task_msgs.times(6)...
        - trials(trial).reach_onset;
    reach_end = find(depth_vel<=0 & speed<0.2,1);
    if ~isempty(reach_end)
        rt_speed_mat(trial,11) = (reach_end+23)/240;
        rt_speed_mat(trial,12) = min(speed);
    else
        rt_speed_mat(trial,11) = nan;
    end
    
end




