function [unit_data, fig] = check_repeated_units(unit_data)

%
% Calculates metrics for each unit, including arm preference, modulation,
% etc
%
%
% INPUTS: 
%
% unit_data - struct: 1 x num_units
%             unit-separated data struct containing hemisphere, brain area,
%             firing rate and other metrics
%
% OUTPUTS:
%
% unit_data - struct: 1 x num_units
%             same as unit_data, but with a new field indicating the
%             correlation coefficient with the previous unit's firing rate.          
%
% fig - figure handle for distribution of pairwise correlation coefficients
%


%% add integer labels to identify which session each unit was recorded in
unit_data(1).session_num = 1;
for unit = 2:length(unit_data)
    if strcmp(unit_data(unit-1).session(1:13),...
              unit_data(unit).session(1:13))
        unit_data(unit).session_num = unit_data(unit-1).session_num;
    else
        unit_data(unit).session_num = unit_data(unit-1).session_num + 1;
    end
end


%% compute unit firing rate correlations within each session
num_sessions = unit_data(end).session_num;
corr_mat = cell(num_sessions,1);
rep_unit_full = logical([]);

for session = 1:num_sessions
    session_data = unit_data([unit_data.session_num]==session);
    % log first unit and use to preallocate the firing rate matrix
    fr_mat = cat_unit_fr(session_data(1));
    fr_mat = repmat(fr_mat, [1,length(session_data)]);
    for unit = 2:length(session_data)
        fr_mat(:,unit) = cat_unit_fr(session_data(unit));
    end
    % compute the correlation matrix and log it for plotting later. include
    % only the upper triangle so that each pair of units only has one
    % non-zero correlation value represented in the matrix
    corr_mat{session} = triu(corrcoef(fr_mat), 1);
    % log boolean flags of repeated units, as determined by a threshold 
    % value of correlation with any other unit
    [row, ~] = find(corr_mat{session} > 0.9);
    rep_unit_session = false(length(session_data),1);
    rep_unit_session(row,1) = true;
    rep_unit_full = [rep_unit_full; rep_unit_session];
end
% log the repeated unit flag in the unit_data struct
new_field = num2cell(rep_unit_full);
[unit_data.repeat] = new_field{:};


%% plot the distribution of pairwise correlations
% combine into a single vector
corr_vals = [];
for session = 1:num_sessions
    holder = triu(corr_mat{session},1);
    corr_vals = [corr_vals; holder(holder~=0)];
end
% plot
fig = figure('Name','Pairwise firing rate correlations',...
    'Position',[600, 600, 700, 300]);
subplot(1,2,1)
histogram(corr_vals, -1:0.05:1, 'normalization','Probability')
hold on
xline(0.9, 'r')
xlim([-1,1])
xlabel('r')
ylabel('Proportion of pairwise correlations')
subplot(1,2,2)
histogram(corr_vals, 0.25:0.05:1, 'normalization','Probability')
hold on
xline(0.9, 'r')
xlim([0.25,1])
xlabel('r')
ylabel('Proportion of pairwise correlations')

end


%% helper functions

function unit_fr_cat = cat_unit_fr(unit_data_row)
% concatenate the firing rate data for a single unit. make the proper
% adjustment so that each sample corresponds to the same timepoint for all
% units, in particular accounting for ipsi/contra referring to the same
% hand depending on which hemisphere the unit is in
if unit_data_row.hem
    unit_fr_cat = [cat_cond_fr(unit_data_row.contra.config(1));...
                   cat_cond_fr(unit_data_row.contra.config(2));...
                   cat_cond_fr(unit_data_row.contra.config(3));...
                   cat_cond_fr(unit_data_row.ipsi.config(1));...
                   cat_cond_fr(unit_data_row.ipsi.config(2));...
                   cat_cond_fr(unit_data_row.ipsi.config(3))];
else
    unit_fr_cat = [cat_cond_fr(unit_data_row.ipsi.config(1));...
                   cat_cond_fr(unit_data_row.ipsi.config(2));...
                   cat_cond_fr(unit_data_row.ipsi.config(3));...
                   cat_cond_fr(unit_data_row.contra.config(1));...
                   cat_cond_fr(unit_data_row.contra.config(2));...
                   cat_cond_fr(unit_data_row.contra.config(3))];
end

    function cond_cat = cat_cond_fr(cond_data)
    % concatenate the data from a single condition (hand, config combo)
        cond_cat = [vertcat(cond_data.target.rest),...
                    vertcat(cond_data.target.prep),...
                    vertcat(cond_data.target.move)];
        cond_cat = cond_cat'; 
        cond_cat = cond_cat(:);
    end
end







