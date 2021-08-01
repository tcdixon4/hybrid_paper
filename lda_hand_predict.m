function [Y, Y_pred, correct, coeffs, unit_area] = ...
    lda_hand_predict(unit_data, trials)


% make sure trials already contains only those which you would like to
% analyze. operates on a single session of `unit_data` and `trials`

tic

%% Prepare feature and class data

[X, Y, unit_area] = prep_lda_hand_predict(unit_data, trials);
Y_non_switch = Y(Y(:,2)==0,1);
Y_switch = Y(Y(:,2)==1,1);
X_non_switch = X(Y(:,2)==0,:,:);
X_switch = X(Y(:,2)==1,:,:);
num_trials_non_switch = size(Y_non_switch,1);
num_trials_switch = size(Y_switch,1);
num_units = length(unit_data);

clearvars unit_data trials


%% prepare for assessing accuracy within each class subtype (hand, target) 
% this is effectively marginalizing over the other class subtypes

% create vectors that combine all classes that fall into a certain subtype
hand{1} = 1:18; % left hand classes, regardless of target, config
hand{2} = 19:36; % right hand classes, regardless of target, config
for t = 6:-1:1 % classes for each target, regardless of hand, config
    target{t} = (0:6:35)+t;
end


%%  perform leave-one-out cross-validation on switching trials

Y_pred_non_switch = zeros(num_trials_non_switch, 78);
correct_hand = zeros(num_trials_non_switch, 78);
correct_target = zeros(num_trials_non_switch, 78);

for test_trial = 1:num_trials_non_switch
    
    % define train trials and target classes
    train_trials = setdiff(1:num_trials_non_switch,test_trial);
    Y_train = Y_non_switch(train_trials);
    
    
    % train a separate model on each window
    parfor train_win = 1:78
        
        % define training feature data
        X_train = X_non_switch(train_trials, :, train_win);
        
        % train the model
        lda_model = fitcdiscr(X_train, Y_train,...
            'Prior','uniform', 'discrimType','pseudoLinear');
               
        % predict the test trial
        X_test = squeeze(X_non_switch(test_trial, :, train_win));
        Y_pred_non_switch(test_trial, train_win) = ...
            predict(lda_model, X_test);  
        
        % determine whether the prediction was correct for each class
        % sub-type
        actual_hand = (Y_non_switch(test_trial)>18) + 1;
        correct_hand(test_trial, train_win) = ismember(...
            Y_pred_non_switch(test_trial, train_win),...
            hand{actual_hand});
        actual_target = mod((Y_non_switch(test_trial)-1),6)+1;
        correct_target(test_trial, train_win) = ismember(...
            Y_pred_non_switch(test_trial, train_win),...
            target{actual_target});
        
    end
    
end

correct.non_switch.hand = correct_hand;
correct.non_switch.target = correct_target;


%% predict switching trials from model trained on all non-switching trials

Y_pred_switch = zeros(num_trials_switch, 78);
correct_hand = zeros(num_trials_switch, 78);
correct_target = zeros(num_trials_switch, 78);

coeffs = cell(78,1);

% train a separate model on each window
for train_win = 1:78
    
    % train the model on non-switch trials, log coeffs
    X_train = X_non_switch(:, :, train_win);
    lda_model = fitcdiscr(X_train, Y_non_switch,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    coeffs{train_win} = lda_model.Coeffs;
    
    % predict all switching trials
    X_test = squeeze(X_switch(:, :, train_win));
    Y_pred_switch(:, train_win) = predict(lda_model, X_test);
    
    % determine whether the prediction was correct for each class
    % sub-type
    for trial = length(Y_switch):-1:1
        actual_hand = (Y_switch(trial)>18) + 1;
        correct_hand(trial, train_win) = ismember(...
            Y_pred_switch(trial, train_win),...
            hand{actual_hand});
        actual_target = mod((Y_switch(trial)-1),6)+1;
        correct_target(trial, train_win) = ismember(...
            Y_pred_switch(trial, train_win),...
            target{actual_target});
    end
end

correct.switch.hand = correct_hand;
correct.switch.target = correct_target;

Y_pred.switch = Y_pred_switch;
Y_pred.non_switch = Y_pred_non_switch;

toc
disp([num2str(num_units), ' units, ', ...
    num2str(num_trials_non_switch+num_trials_switch), ' trials'])


