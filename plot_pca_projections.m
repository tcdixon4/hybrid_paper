function [] = plot_pca_projections(unit_data, norm_method)


%% Prepare data

% isolate only non-repeated units
unit_data = unit_data([unit_data.repeat]==false);

% prepare data from each sub-population separately
[X_l_hem, ~] = ...
    prep_pca_trial_avg(unit_data([unit_data.hem]==0), norm_method);
[X_r_hem, ~] = ...
    prep_pca_trial_avg(unit_data([unit_data.hem]==1), norm_method);

% create an indexing array for preferred arm
% compute arm preferences with the stationary hand center config
[arm_pref, ~, ~, ~] = calc_limb_dedication(unit_data, 0);
% r_pref_idx.prep = arm_pref.prep>0;
r_pref = arm_pref.move>0;
[X_l_pref, X_r_pref] = extract_pref_pop(X_l_hem, X_r_hem, r_pref);


%% Train models and compute projections

% left hem, ipsi trained
[l_hem_l_train_native.prep,...
    l_hem_l_train_l_hand.prep, l_hem_l_train_r_hand.prep]...
    = compute_projections(X_l_hem, 'prep', 'left');
[l_hem_l_train_native.move,...
    l_hem_l_train_l_hand.move, l_hem_l_train_r_hand.move]...
    = compute_projections(X_l_hem, 'move', 'left');

% right hem, ipsi trained
[r_hem_r_train_native.prep,...
    r_hem_r_train_r_hand.prep, r_hem_r_train_l_hand.prep]...
    = compute_projections(X_r_hem, 'prep', 'right');
[r_hem_r_train_native.move,...
    r_hem_r_train_r_hand.move, r_hem_r_train_l_hand.move]...
    = compute_projections(X_r_hem, 'move', 'right');

% left hem, contra trained
[l_hem_r_train_native.prep,...
    l_hem_r_train_r_hand.prep, l_hem_r_train_l_hand.prep]...
    = compute_projections(X_r_hem, 'prep', 'right');
[l_hem_r_train_native.move,...
    l_hem_r_train_r_hand.move, l_hem_r_train_l_hand.move]...
    = compute_projections(X_r_hem, 'move', 'right');

% right hem, contra trained
[r_hem_l_train_native.prep,...
    r_hem_l_train_l_hand.prep, r_hem_l_train_r_hand.prep]...
    = compute_projections(X_l_hem, 'prep', 'left');
[r_hem_l_train_native.move,...
    r_hem_l_train_l_hand.move, r_hem_l_train_r_hand.move]...
    = compute_projections(X_l_hem, 'move', 'left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% left pref, left hand (pref) trained
[l_pref_l_train_native.prep,...
    l_pref_l_train_l_hand.prep, l_pref_l_train_r_hand.prep]...
    = compute_projections(X_l_pref, 'prep', 'left');
[l_pref_l_train_native.move,...
    l_pref_l_train_l_hand.move, l_pref_l_train_r_hand.move]...
    = compute_projections(X_l_pref, 'move', 'left');

% right pref, right hand (pref) trained
[r_pref_r_train_native.prep,...
    r_pref_r_train_r_hand.prep, r_pref_r_train_l_hand.prep]...
    = compute_projections(X_r_pref, 'prep', 'right');
[r_pref_r_train_native.move,...
    r_pref_r_train_r_hand.move, r_pref_r_train_l_hand.move]...
    = compute_projections(X_r_pref, 'move', 'right');

% left pref, right hand (non-pref) trained
[l_pref_r_train_native.prep,...
    l_pref_r_train_r_hand.prep, l_pref_r_train_l_hand.prep]...
    = compute_projections(X_l_pref, 'prep', 'right');
[l_pref_r_train_native.move,...
    l_pref_r_train_r_hand.move, l_pref_r_train_l_hand.move]...
    = compute_projections(X_l_pref, 'move', 'right');

% right pref, left hand (non-pref) trained
[r_pref_l_train_native.prep,...
    r_pref_l_train_l_hand.prep, r_pref_l_train_r_hand.prep]...
    = compute_projections(X_r_pref, 'prep', 'left');
[r_pref_l_train_native.move,...
    r_pref_l_train_l_hand.move, r_pref_l_train_r_hand.move]...
    = compute_projections(X_r_pref, 'move', 'left');



%% Plot results

% % left hem, ipsi trained
% plot_projections(l_hem_l_train_native, l_hem_l_train_l_hand,...
%                  l_hem_l_train_r_hand, {1:3, 1:3}, 'left hem, ipsi trained')
%              
% % right hem, ipsi trained
% plot_projections(r_hem_r_train_native, r_hem_r_train_r_hand,...
%                  r_hem_r_train_l_hand, {1:3, 1:3}, 'right hem, ipsi trained')
%              
% % left hem, contra trained
% plot_projections(l_hem_r_train_native, l_hem_r_train_r_hand,...
%                  l_hem_l_train_l_hand, {1:3, 1:3}, 'left hem, contra trained')
%              
% % right hem, contra trained
% plot_projections(r_hem_l_train_native, r_hem_l_train_l_hand,...
%                  r_hem_l_train_r_hand, {1:3, 1:3}, 'right hem, contra trained')
%              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% left pref, left hand (pref) trained
plot_projections(l_pref_l_train_native, l_pref_l_train_l_hand,...
                 l_pref_l_train_r_hand, {[2,4,3], [4,3,2]}, 'left pref, left (pref) trained')
             
% right pref, right hand (pref) trained
plot_projections(r_pref_r_train_native, r_pref_r_train_r_hand,...
                 r_pref_r_train_l_hand, {1:3, 1:3}, 'right pref, right (pref) trained')
             
% left pref, right hand (non-pref) trained
plot_projections(l_pref_r_train_native, l_pref_r_train_r_hand,...
                 l_pref_r_train_l_hand, {1:3, 1:3}, 'left pref, right (non-pref) trained')
             
% right pref, left hand (non-pref) trained
plot_projections(r_pref_l_train_native, r_pref_l_train_l_hand,...
                 r_pref_l_train_r_hand, {1:3, 1:3}, 'right pref, left (non-pref) trained')



end


%% Helper functions


function [X_l_pref, X_r_pref] = extract_pref_pop(X_l_hem, X_r_hem, r_pref)

for config = 1:3  
    % recombine the left and right hemisphere data
    l_hand_holder.rest = ...
        [X_l_hem.l_hand.config(config).rest;...
        X_r_hem.l_hand.config(config).rest];
    l_hand_holder.prep = ...
        [X_l_hem.l_hand.config(config).prep;...
        X_r_hem.l_hand.config(config).prep];
    l_hand_holder.move = ...
        [X_l_hem.l_hand.config(config).move;...
        X_r_hem.l_hand.config(config).move];
    r_hand_holder.rest = ...
        [X_l_hem.r_hand.config(config).rest;...
        X_r_hem.r_hand.config(config).rest];
    r_hand_holder.prep = ...
        [X_l_hem.r_hand.config(config).prep;...
        X_r_hem.r_hand.config(config).prep];
    r_hand_holder.move = ...
        [X_l_hem.r_hand.config(config).move;...
        X_r_hem.r_hand.config(config).move];
    
    % separate right hand preferring and left hand preferring populations
    X_l_pref.l_hand.config(config).rest = l_hand_holder.rest(~r_pref,:,:);
    X_l_pref.r_hand.config(config).rest = r_hand_holder.rest(~r_pref,:,:);
    X_r_pref.l_hand.config(config).rest = l_hand_holder.rest(r_pref,:,:);
    X_r_pref.r_hand.config(config).rest = r_hand_holder.rest(r_pref,:,:);
    
    X_l_pref.l_hand.config(config).prep = l_hand_holder.prep(~r_pref,:,:);
    X_l_pref.r_hand.config(config).prep = r_hand_holder.prep(~r_pref,:,:);
    X_r_pref.l_hand.config(config).prep = l_hand_holder.prep(r_pref,:,:);
    X_r_pref.r_hand.config(config).prep = r_hand_holder.prep(r_pref,:,:);
    
    X_l_pref.l_hand.config(config).move = l_hand_holder.move(~r_pref,:,:);
    X_l_pref.r_hand.config(config).move = r_hand_holder.move(~r_pref,:,:);
    X_r_pref.l_hand.config(config).move = l_hand_holder.move(r_pref,:,:);
    X_r_pref.r_hand.config(config).move = r_hand_holder.move(r_pref,:,:);
end


end


function [native_score, cross_score_a, cross_score_b] =...
    compute_projections(X, phase, train_hand)

% select the indicated datasets. Train on the eccentric config, but for
% comparison of the two arms project the center config data
if strcmp(phase, 'prep')
    if strcmp(train_hand, 'left')
        X_train = cat(2, X.l_hand.config(3).rest,...
                       X.l_hand.config(3).prep(:,1:end-1,:));
        X_cross_a = cat(2, X.l_hand.config(2).rest,...
                      X.l_hand.config(2).prep(:,1:end-1,:));
        X_cross_b = cat(2, X.r_hand.config(2).rest,...
                      X.r_hand.config(2).prep(:,1:end-1,:));
    else
        X_train = cat(2, X.r_hand.config(1).rest,...
                       X.r_hand.config(1).prep(:,1:end-1,:));
        X_cross_a = cat(2, X.r_hand.config(2).rest,...
                      X.r_hand.config(2).prep(:,1:end-1,:));
        X_cross_b = cat(2, X.l_hand.config(2).rest,...
                      X.l_hand.config(2).prep(:,1:end-1,:));
    end
elseif strcmp(phase, 'move')
    if strcmp(train_hand, 'left')
        X_train = X.l_hand.config(3).move;
        X_cross_a = X.l_hand.config(2).move;
        X_cross_b = X.r_hand.config(2).move;
    else
        X_train = X.r_hand.config(1).move;
        X_cross_a = X.r_hand.config(2).move;
        X_cross_b = X.l_hand.config(2).move;
    end
else
    error('Improper phase parameter.')
end

% organize for PCA
X_train = [X_train(:,:,1), X_train(:,:,2), X_train(:,:,3),...
            X_train(:,:,4), X_train(:,:,5), X_train(:,:,6)];
X_cross_a = [X_cross_a(:,:,1), X_cross_a(:,:,2), X_cross_a(:,:,3),...
           X_cross_a(:,:,4), X_cross_a(:,:,5), X_cross_a(:,:,6)];
X_cross_b = [X_cross_b(:,:,1), X_cross_b(:,:,2), X_cross_b(:,:,3),...
           X_cross_b(:,:,4), X_cross_b(:,:,5), X_cross_b(:,:,6)];
       
% train model and log the within-sample (native) projection
[coeff, native_score] = pca(X_train', 'Centered',false);

% compute the across-sample (cross) projections
cross_score_a = X_cross_a'*coeff;
cross_score_b = X_cross_b'*coeff;

end


function [] = plot_projections(native_score, cross_a, cross_b, pcs, ttl)

figure('Position',[-1247,565,1000,400], 'Name',['3D projections - ', ttl])
subplot(1,2,1)
plot_3d_projections(native_score.prep, pcs{1})
title('Instruct')
subplot(1,2,2)
plot_3d_projections(native_score.move, pcs{2})
title('Move')

% figure('Position',[-1247,0,400,1000], 'Name',['Instruct - ', ttl])
% plot_projection_comparison(cross_a.prep, cross_b.prep)
% 
% figure('Position',[-1247,0,400,1000], 'Name',['Move - ', ttl])
% plot_projection_comparison(cross_a.move, cross_b.move)

end


function [] = plot_3d_projections(native_score, pcs)
p = cell(6,1);
k = size(native_score,1)/6;
if k==41 % Instruct phase x-axis information
    fade = create_color_gradients(15,k-16,1);
else % Move phase x-axis information
    fade = create_color_gradients(3,k-7,4);
end
for i = 1:6
    p{i} = plot3(native_score(k*(i-1)+1:k*i,pcs(1)),...
                 native_score(k*(i-1)+1:k*i,pcs(2)),...
                 native_score(k*(i-1)+1:k*i,pcs(3)), 'LineWidth',5);
    hold on
    marker{i} = plot3(native_score(k*i,pcs(1)),...
                      native_score(k*i,pcs(2)),...
                      native_score(k*i,pcs(3)), 'ko',...
                      'MarkerFaceColor',uint8(fade{i}(1:3,end)),...
                      'MarkerEdgeColor','none', 'MarkerSize',10);
end
grid on
xlabel(['PC', num2str(pcs(1))]); 
ylabel(['PC', num2str(pcs(2))]); 
zlabel(['PC', num2str(pcs(3))])

drawnow
set(p{1}.Edge, 'ColorBinding','interpolated', 'ColorData',uint8(fade{1}))
set(p{2}.Edge, 'ColorBinding','interpolated', 'ColorData',uint8(fade{2}))
set(p{3}.Edge, 'ColorBinding','interpolated', 'ColorData',uint8(fade{3}))
set(p{4}.Edge, 'ColorBinding','interpolated', 'ColorData',uint8(fade{4}))
set(p{5}.Edge, 'ColorBinding','interpolated', 'ColorData',uint8(fade{5}))
set(p{6}.Edge, 'ColorBinding','interpolated', 'ColorData',uint8(fade{6}))

legend([marker{:}],...
    {'Target 1','Target 2','Target 3','Target 4','Target 5','Target 6'})

view(53, 33)
xticklabels([])
yticklabels([])
zticklabels([])
% xlim([-1,1])
% xticks(-1:1:1)
% ylim([-1,1])
% yticks(-1:1:1)
% zlim([-1,1])
% zticks(-1:1:1)
ax = gca;
ax.Clipping = 'off';
pbaspect([1,1,1])
end


function [] = plot_both_3d_projections(cross_a, cross_b, pcs)
p = cell(6,1);
k = size(cross_a,1)/6;
if k==41 % Instruct phase x-axis information
    fade = create_color_gradients(15,k-16,1);
else % Move phase x-axis information
    fade = create_color_gradients(3,k-7,4);
end
for i = 1:6
    p{i} = plot3(cross_a(k*(i-1)+1:k*i,pcs(1)),...
                 cross_a(k*(i-1)+1:k*i,pcs(2)),...
                 -cross_a(k*(i-1)+1:k*i,pcs(3)),...
                 'LineWidth',5, 'Color',uint8(fade{i}(1:3,end)));
    hold on
    p{i} = plot3(cross_b(k*(i-1)+1:k*i,pcs(1)),...
                 cross_b(k*(i-1)+1:k*i,pcs(2)),...
                 -cross_b(k*(i-1)+1:k*i,pcs(3)),...
                 'LineWidth',4, 'Color',uint8(fade{i}(1:3,end)), 'LineStyle',':');
    marker{i} = plot3(cross_a(k*i,pcs(1)),...
                      cross_a(k*i,pcs(2)),...
                      -cross_a(k*i,pcs(3)), 'ko',...
                      'MarkerFaceColor',uint8(fade{i}(1:3,end)),...
                      'MarkerEdgeColor','none', 'MarkerSize',10);
    marker{i} = plot3(cross_b(k*i,pcs(1)),...
                      cross_b(k*i,pcs(2)),...
                      -cross_b(k*i,pcs(3)), 'ko',...
                      'MarkerFaceColor',uint8(fade{i}(1:3,end)),...
                      'MarkerEdgeColor','none', 'MarkerSize',10);
end
grid on
xlabel(['PC', num2str(pcs(1))]); 
ylabel(['PC', num2str(pcs(2))]); 
zlabel(['PC', num2str(pcs(3))])

legend([marker{:}],...
    {'Target 1','Target 2','Target 3','Target 4','Target 5','Target 6'})

view(53, 33)
xticklabels([])
yticklabels([])
zticklabels([])
% xlim([-1,1])
% xticks(-1:1:1)
% ylim([-1,1])
% yticks(-1:1:1)
% zlim([-1,1])
% zticks(-1:1:1)
ax = gca;
ax.Clipping = 'off';
pbaspect([1,1,1])
end


function [] = plot_projection_comparison(cross_a, cross_b)
k = size(cross_a,1)/6;
if k==41 % Instruct phase x-axis information
    fade = create_color_gradients(6,k-7,1);
    x = -300:20:500;
    xt = [-300,0,500];
    xtl = {'','Instruct','+500ms'};
    xl = [-320,520];
else % Move phase x-axis information
    fade = create_color_gradients(1,k-9,8);
    x = -200:20:500;
    xt = [-200,0,500];
    xtl = {'','Move','+500ms'};
    xl = [-220,620];
end

for pc = 1:4
    subplot(4,1,pc)
    for i = 1:6
        p{i} = plot(x, cross_a(k*(i-1)+1:k*i,pc), 'LineWidth',5,...
            'Color',uint8(fade{i}(1:3,end)));
        hold on
        p{i+6} = plot(x, cross_b(k*(i-1)+1:k*i,pc), 'LineWidth',1,...
            'Color',uint8(fade{i}(1:3,end)));
        xlim(xl)
        xticks(xt)
        xticklabels({'','',''})
        xline(0);
        yticks([])
        box off
    end
    title(['PC ', num2str(pc)])
end
% legend([p{1}, p{7}], {'Native-hand', 'Cross-hand'});
xticklabels(xtl);
end


function [fade] = create_color_gradients(b1, b2, b3)

fade{1} = [...
    [180*ones(1,b1), linspace(180,0,b2),   zeros(1,b3)];...
    [180*ones(1,b1), linspace(180,180,b2), 180*ones(1,b3)];...
    [180*ones(1,b1), linspace(180,0,b2),   zeros(1,b3)];...
    ones(1,b1+b2+b3)];

fade{2} = [...
    [180*ones(1,b1), linspace(180,0,b2),   zeros(1,b3)];...
    [180*ones(1,b1), linspace(180,180,b2), 180*ones(1,b3)];...
    [180*ones(1,b1), linspace(180,180,b2), 180*ones(1,b3)];...
    ones(1,b1+b2+b3)];

fade{3} = [...
    [180*ones(1,b1), linspace(180,0,b2),   zeros(1,b3)];...
    [180*ones(1,b1), linspace(180,0,b2),   zeros(1,b3)];...
    [180*ones(1,b1), linspace(180,180,b2), 180*ones(1,b3)];...
    ones(1,b1+b2+b3)];

fade{4} = [...
    [180*ones(1,b1), linspace(180,180,b2), 180*ones(1,b3)];...
    [180*ones(1,b1), linspace(180,180,b2), 180*ones(1,b3)];...
    [180*ones(1,b1), linspace(180,0,b2),   zeros(1,b3)];...
    ones(1,b1+b2+b3)];

fade{5} = [...
    [180*ones(1,b1), linspace(180,180,b2), 180*ones(1,b3)];...
    [180*ones(1,b1), linspace(180,0,b2),   zeros(1,b3)];...
    [180*ones(1,b1), linspace(180,0,b2),   zeros(1,b3)];...
    ones(1,b1+b2+b3)];

fade{6} = [...
    [180*ones(1,b1), linspace(180,180,b2), 180*ones(1,b3)];...
    [180*ones(1,b1), linspace(180,0,b2),   zeros(1,b3)];...
    [180*ones(1,b1), linspace(180,180,b2), 180*ones(1,b3)];...
    ones(1,b1+b2+b3)];

end






