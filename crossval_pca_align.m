function [alignment, explained_a, explained_b] = ...
    crossval_pca_align(X_a, X_b, same_cond, num_comp, iterations)

%
% Computes alignment of p-dimensional PCA subspaces fit to datasets A and B
% using Monte Carlo cross-validation.
%
%
% INPUTS: 
%
% X_a - unit-separated data struct containing hemisphere, brain area,
%             firing rate and other metrics
%             (struct: 1 x num_units)
%
% same_cond - specify whether comparing trials using the same or different
%             hands. binary
%
% num_comp - number of principle components to use. integer
%
% iterations - number of iterations to randomly divide trials and compute
%              alignment (this is stochastic for both within and across 
%              hand comparisons). integer
%
%
% OUTPUTS:
%
% alignment - alignment across time for trials using the same hand
%             (vector: 1 x iterations) 
%
% explained_a - NON-cross-validated variance accounted for by each PC in
%               model fit to dataset A
%               (matrix: num_units x num_iterations)
% 
% explained_b - NON-cross-validated variance accounted for by each PC in
%               model fit to dataset B
%               (matrix: num_units x num_iterations)
% 


%% Setup 

num_trials_a = size(X_a,1);
num_trials_b = size(X_b,1);
num_units = size(X_a,2);

alignment = zeros(1,iterations);
explained_a = zeros(num_units,iterations);
explained_b = zeros(num_units,iterations);

%% Divide up data, train PC model on both portions, calculate alignment of 
% the two. 

if same_cond
    parfor iter = 1:iterations
        warning off;
        a_idx = randperm(num_trials_a, floor(num_trials_a/2)); %randomly select half the trials rounded down
        b_idx = setdiff(1:num_trials_a, a_idx); %the remainder are assigned to the other group
        b_idx = b_idx(1:length(a_idx)); %make groups equal size
        %train PC models on each dataset and compute alignment
        sampled_a = X_a(a_idx,:,:);
        sampled_b = X_b(b_idx,:,:);
        %if including multiple samples from each trial, concatenate trials
        %after the trials have been selected
        if size(sampled_a,3) > 1
            sampled_a = reshape(...
                permute(sampled_a,[3,1,2]),...
                [],num_units);
            sampled_b = reshape(...
                permute(sampled_b,[3,1,2]),...
                [],num_units);
        end
        [coeff_a,~,~,~,exp_a,~] = ...
            pca(sampled_a,'NumComponents',num_comp, 'Centered',false);
        [coeff_b,~,~,~,exp_b,~] = ...
            pca(sampled_b,'NumComponents',num_comp, 'Centered',false);
        explained_a(:,iter) = [exp_a; ones(num_units-length(exp_a),1)];
        explained_b(:,iter) = [exp_b; ones(num_units-length(exp_b),1)];
        alignment(iter) = trace(cov(sampled_a*coeff_a*(coeff_a')*coeff_b))/...
            trace(cov(sampled_a*coeff_a));
    end
    
else
    num_trials = min(num_trials_a, num_trials_b);
    for iter = 1:iterations
        a_idx = randperm(num_trials_a, floor(num_trials/2)); %randomly select half the trials rounded down
        b_idx = randperm(num_trials_b, floor(num_trials/2)); %randomly select half the trials rounded down
        sampled_a = X_a(a_idx,:,:);
        sampled_b = X_b(b_idx,:,:);
        %if including multiple samples from each trial, concatenate trials
        %after the trials have been selected
        if size(sampled_a,3) > 1
            sampled_a = reshape(...
                permute(sampled_a,[3,1,2]),...
                [],num_units);
            sampled_b = reshape(...
                permute(sampled_b,[3,1,2]),...
                [],num_units);
        end
        [coeff_a,~,~,~,exp_a,~] = ...
            pca(sampled_a,'NumComponents',num_comp, 'Centered',false);
        [coeff_b,~,~,~,exp_b,~] = ...
            pca(sampled_b,'NumComponents',num_comp, 'Centered',false);
        explained_a(:,iter) = [exp_a; ones(num_units-length(exp_a),1)];
        explained_b(:,iter) = [exp_b; ones(num_units-length(exp_b),1)];
        alignment(iter) = trace(cov(sampled_a*coeff_a*(coeff_a')*coeff_b))/...
            trace(cov(sampled_a*coeff_a));
    end
end



