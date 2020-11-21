% Master script for time varying LDA analyses on the PP datasets
function [Y, Y_pred] = pp_lda_pref(unit_data)

%% Prepare feature and class data and preallocate output

[X, Y] = prep_lda(unit_data);
num_trials_l_hand = length(Y.l_hand);
num_trials_r_hand = length(Y.r_hand);


%%  preallocate predicted and actual matrices

Y_pred_l_pref_l_hand = zeros(78, 78, num_trials_l_hand);
Y_pred_l_pref_r_hand = zeros(78, 78, num_trials_r_hand);
Y_pred_r_pref_l_hand = zeros(78, 78, num_trials_l_hand);
Y_pred_r_pref_r_hand = zeros(78, 78, num_trials_r_hand);


%%  perform leave-one-out cross-validation on left hand trials

for test_trial = 1:num_trials_l_hand
    
    %% define train trials and target classes
    train_trials = setdiff(1:num_trials_l_hand,test_trial);
    Y_l_hand_train = Y.l_hand(train_trials);
    
    
    %% train a separate model on each window
    parfor train_win = 1:78
        %% left preferring model
        % define training feature data
        X_train = X.l_pref.l_hand(train_trials, :, train_win);
        % train the model and predict the test trial
        lda_model_l_pref = fitcdiscr(X_train, Y_l_hand_train,...
            'Prior','uniform', 'discrimType','pseudoLinear');

        
        %% right preferring model
        % define training and testing feature data
        X_train = X.r_pref.l_hand(train_trials, :, train_win);               
        % train the model and predict the test trial
        lda_model_r_pref = fitcdiscr(X_train, Y_l_hand_train,...
            'Prior','uniform', 'discrimType','pseudoLinear');
        
        
        %% predict the test trial for each window
        X_test = squeeze(X.l_pref.l_hand(test_trial, :, :))';
        Y_pred_l_pref_l_hand(train_win, :, test_trial) = ...
            predict(lda_model_l_pref, X_test);
        
        X_test = squeeze(X.r_pref.l_hand(test_trial, :, :))';
        Y_pred_r_pref_l_hand(train_win, :, test_trial) = ...
            predict(lda_model_r_pref, X_test);
        
        
    end
    
end



%% perform leave-one-out cross-validation on right hand trials

for test_trial = 1:num_trials_r_hand
    
    %% define train trials and target classes
    train_trials = setdiff(1:num_trials_r_hand,test_trial);
    Y_r_hand_train = Y.r_hand(train_trials);
    
    
    %% train a separate model on each window
    parfor train_win = 1:78
        %% left preferring model
        % define training feature data
        X_train = X.l_pref.r_hand(train_trials, :, train_win);
        % train the model and predict the test trial
        lda_model_l_pref = fitcdiscr(X_train, Y_r_hand_train,...
            'Prior','uniform', 'discrimType','pseudoLinear');

        
        %% right preferring model
        % define training and testing feature data
        X_train = X.r_pref.r_hand(train_trials, :, train_win);               
        % train the model and predict the test trial
        lda_model_r_pref = fitcdiscr(X_train, Y_r_hand_train,...
            'Prior','uniform', 'discrimType','pseudoLinear');
        
        
        %% predict the test trial for each window
        X_test = squeeze(X.l_pref.r_hand(test_trial, :, :))';
        Y_pred_l_pref_r_hand(train_win, :, test_trial) = ...
            predict(lda_model_l_pref, X_test);
        
        X_test = squeeze(X.r_pref.r_hand(test_trial, :, :))';
        Y_pred_r_pref_r_hand(train_win, :, test_trial) = ...
            predict(lda_model_r_pref, X_test);
        
        
    end
       
end


%% convert predictions to structure array
Y_pred.l_pref.l_hand = Y_pred_l_pref_l_hand;
Y_pred.l_pref.r_hand = Y_pred_l_pref_r_hand;
Y_pred.r_pref.l_hand = Y_pred_r_pref_l_hand;
Y_pred.r_pref.r_hand = Y_pred_r_pref_r_hand;




