function [err, dim_est, total_ss] = estimate_dim_pca(neural_cell, phase)

%
% Estimates the dimensionality of PCA subspaces fit to data during trials
% using a single arm.
%
% approach based on response from user @amoeba on StackExchange post below:
% https://stats.stackexchange.com/questions/93845/how-to-perform-cross-validation-for-pca-to-determine-the-number-of-principal-com
%
%
% INPUTS: 
%
% neural_cell - time-locked and trial-separated firing rate traces for all
%               simultaneously recorded units
%               (cell array: num_trials x 1
%                nested matrices: num_samples x num_units)
%
% phase - selected task phase to use. 'rest', 'prep', 'move', or 'all'
%
%
% OUTPUTS:
%
% err - vector: 12x1, mean-centered and total sum of squares then
%       prediction error using 1:10 PC's
% 
% dim_est - integer: dimensionality estimate that minimizes error
%
% total_ss - scalar: total sum of squares
%


%% Extract just the data from the phase of interest

switch phase
    case 'rest'
        idx = 11:26;
    case 'prep'
        idx = 37:52;
    case 'move'
        idx = 63:78;
    case 'all'
        idx = [];
end

if ~isempty(idx)
    for trial = 1:length(neural_cell)
        neural_cell{trial} = neural_cell{trial}(idx,:);
    end
end


%% calculate cross-validated and native predictions

% full data matrix
X_true = cell2mat(neural_cell);

% iterate over trials, training leave-one-trial-out CV PCA models
for trial = 1:length(neural_cell)
    warning off
    % assign training trials
    train_idx = setdiff(1:length(neural_cell),trial);
    X_train = cell2mat(neural_cell(train_idx));
    % train PCA model using all but one trial
    [P, ~,~,~,~, ~] = pca(X_train, 'Centered', false);
    X_test = neural_cell{trial};
    
    % iterate over model sizes (number of components)
    for p = 10:-1:1
        
        % iterate over units, estimating each from the PC projection of the
        % rest
        for unit = size(X_true,2):-1:1
            
            % assign predictor units for imputing the held out unit
            impute_idx = setdiff(1:size(X_true,2),unit);
            % Create a new loading matrix that ignores the imputed unit by
            % removing the coefficients for the imputed unit, taking the
            % pseudoinverse to create a new orthogonal matrix, and then
            % taking the transpose
            P_impute = pinv(P(impute_idx,1:p))';
            % Compute the PC score with the leftover units and predict the
            % imputed unit using the transpose of the original loading 
            % matrix
            T = X_test(:,impute_idx)*P_impute;
            X_test_est = T*(P(:,1:p)');
            X_est_cv{trial,p}(:,unit) = X_test_est(:,unit);
            
        end
        
    end
    
end



%% calculate cross-validated prediction error and native R^2

for p = 10:-1:1
    
    X_est_cv_p = cell2mat(X_est_cv(:,p));
    err(p+2,1) = sum( (X_est_cv_p(:)-X_true(:)).^2 );
    
end


X_mean_centered = X_true - mean(X_true);
err(1,1) = sum(X_mean_centered(:).^2);
total_ss = sum(X_true(:).^2);
err(2,1) = total_ss;

[~,dim_est] = min(err(3:end));









