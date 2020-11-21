% Master script for time varying LDA analyses on the PP datasets
function [Y, Y_pred] = pp_lda_pref_epochs(unit_data)

%% Prepare feature and class data

[X, Y] = prep_lda_epochs(unit_data);
num_trials_l_hand = length(Y.l_hand);
num_trials_r_hand = length(Y.r_hand);


%%  preallocate prediction matrices

Y_pred.l_pref.l_hand.rest = zeros(num_trials_l_hand, 1);
Y_pred.l_pref.l_hand.prep = zeros(num_trials_l_hand, 1);
Y_pred.l_pref.l_hand.move = zeros(num_trials_l_hand, 1);
Y_pred.l_pref.r_hand.rest = zeros(num_trials_r_hand, 1);
Y_pred.l_pref.r_hand.prep = zeros(num_trials_r_hand, 1);
Y_pred.l_pref.r_hand.move = zeros(num_trials_r_hand, 1);
Y_pred.r_pref = Y_pred.l_pref;


%%  perform leave-one-out cross-validation on left hand trials

for test_idx = 1:num_trials_l_hand
    
    %% define training indices and classes
    train_idx = setdiff(1:num_trials_l_hand,test_idx);
    Y_l_hand_train = Y.l_hand(train_idx);
    
    
    %% Rest, l hand, l pref    
    % prepare feature data 
    X_train = X.l_pref.l_hand.rest(train_idx,:);
    X_test = X.l_pref.l_hand.rest(test_idx,:);
    
    % train the model and predict the test trial
    lda_model = fitcdiscr(X_train, Y_l_hand_train,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    Y_pred.l_pref.l_hand.rest(test_idx) = predict(lda_model, X_test);    
    
    %% Rest, l hand, r pref
    % prepare feature data 
    X_train = X.r_pref.l_hand.rest(train_idx,:);
    X_test = X.r_pref.l_hand.rest(test_idx,:);
    
    % train the model and predict the test trial
    lda_model = fitcdiscr(X_train, Y_l_hand_train,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    Y_pred.r_pref.l_hand.rest(test_idx) = predict(lda_model, X_test);
    
    
    %% Prep, l hand, l pref    
    % prepare feature data 
    X_train = X.l_pref.l_hand.prep(train_idx,:);
    X_test = X.l_pref.l_hand.prep(test_idx,:);
    
    % train the model and predict the test trial
    lda_model = fitcdiscr(X_train, Y_l_hand_train,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    Y_pred.l_pref.l_hand.prep(test_idx) = predict(lda_model, X_test);    
    
    %% Rest, l hand, r pref
    % prepare feature data 
    X_train = X.r_pref.l_hand.prep(train_idx,:);
    X_test = X.r_pref.l_hand.prep(test_idx,:);
    
    % train the model and predict the test trial
    lda_model = fitcdiscr(X_train, Y_l_hand_train,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    Y_pred.r_pref.l_hand.prep(test_idx) = predict(lda_model, X_test);
    
    
    %% Move, l hand, l pref    
    % prepare feature data 
    X_train = X.l_pref.l_hand.move(train_idx,:);
    X_test = X.l_pref.l_hand.move(test_idx,:);
    
    % train the model and predict the test trial
    lda_model = fitcdiscr(X_train, Y_l_hand_train,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    Y_pred.l_pref.l_hand.move(test_idx) = predict(lda_model, X_test);    
    
    %% Rest, l hand, r pref
    % prepare feature data 
    X_train = X.r_pref.l_hand.move(train_idx,:);
    X_test = X.r_pref.l_hand.move(test_idx,:);
    
    % train the model and predict the test trial
    lda_model = fitcdiscr(X_train, Y_l_hand_train,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    Y_pred.r_pref.l_hand.move(test_idx) = predict(lda_model, X_test);
    

end


%% perform leave-one-out cross-validation on right hand trials

for test_idx = 1:num_trials_r_hand
    
    %% define training indices and classes
    train_idx = setdiff(1:num_trials_r_hand,test_idx);
    Y_r_hand_train = Y.r_hand(train_idx);
    
    
    %% Rest, r hand, l pref
    % prepare feature data 
    X_train = X.l_pref.r_hand.rest(train_idx,:);
    X_test = X.l_pref.r_hand.rest(test_idx,:);
    
    % train the model and predict the test trial
    lda_model = fitcdiscr(X_train, Y_r_hand_train,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    Y_pred.l_pref.r_hand.rest(test_idx) = predict(lda_model, X_test);
    
    %% Rest, r hand, r pref
    % prepare feature data 
    X_train = X.r_pref.r_hand.rest(train_idx,:);
    X_test = X.r_pref.r_hand.rest(test_idx,:);
    
    % train the model and predict the test trial
    lda_model = fitcdiscr(X_train, Y_r_hand_train,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    Y_pred.r_pref.r_hand.rest(test_idx) = predict(lda_model, X_test);
    
    
    %% Prep, r hand, l pref
    % prepare feature data 
    X_train = X.l_pref.r_hand.prep(train_idx,:);
    X_test = X.l_pref.r_hand.prep(test_idx,:);
    
    % train the model and predict the test trial
    lda_model = fitcdiscr(X_train, Y_r_hand_train,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    Y_pred.l_pref.r_hand.prep(test_idx) = predict(lda_model, X_test);
    
    %% Prep, r hand, r pref
    % prepare feature data 
    X_train = X.r_pref.r_hand.prep(train_idx,:);
    X_test = X.r_pref.r_hand.prep(test_idx,:);
    
    % train the model and predict the test trial
    lda_model = fitcdiscr(X_train, Y_r_hand_train,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    Y_pred.r_pref.r_hand.prep(test_idx) = predict(lda_model, X_test);

    
    %% Move, r hand, l pref
    % prepare feature data 
    X_train = X.l_pref.r_hand.move(train_idx,:);
    X_test = X.l_pref.r_hand.move(test_idx,:);
    
    % train the model and predict the test trial
    lda_model = fitcdiscr(X_train, Y_r_hand_train,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    Y_pred.l_pref.r_hand.move(test_idx) = predict(lda_model, X_test);
    
    %% Move, r hand, r pref
    % prepare feature data 
    X_train = X.r_pref.r_hand.move(train_idx,:);
    X_test = X.r_pref.r_hand.move(test_idx,:);
    
    % train the model and predict the test trial
    lda_model = fitcdiscr(X_train, Y_r_hand_train,...
        'Prior','uniform', 'discrimType','pseudoLinear');
    Y_pred.r_pref.r_hand.move(test_idx) = predict(lda_model, X_test);
    

end

