
% do only pain animals, to see if there are clusters with many aniamls

animal_list = unique(animal_list_used_after_analysis, 'stable');
p_a= endsWith(animal_list, '_F');


cond_inds_F = zeros(1,length(analysisstruct.condition_inds)); % sorting per animal
new_iid = 0;
for iid = 1:length(animal_list)
    
    animal_ID = animal_list{iid};
    idx = ismember(animal_list_used_after_analysis, animal_ID);
    if p_a (iid)
        new_iid =new_iid + iid;
    else
        new_iid = 0;
    end
    % new_iid = p_a (iid) *  double(p_a (iid)) + iid;
    cond_inds_F(idx) = new_iid;
end

cond_inds_S= zeros(1,length(analysisstruct.condition_inds)); % sorting per animal
new_iid = 0;
for iid = 1:length(animal_list)
    
    animal_ID = animal_list{iid};
    idx = ismember(animal_list_used_after_analysis, animal_ID);
    if ~p_a (iid)
        new_iid =new_iid + iid;
    else
        new_iid = 0;
    end
    % new_iid = p_a (iid) *  double(p_a (iid)) + iid;
    cond_inds_S(idx) = new_iid;
end


%% Figure
figure
gscatter(zvals(:,1), zvals(:,2), cond_inds_F)
% plot(zvals(:,1),zvals(:,2),'ob','MarkerFaceColor','b', 'MarkerSize',2)
title({['Granu: ',num2str(analysisparams.tsnegranularity)], ['Perp: ', num2str(perplexity)]})
set(gcf,'Position',([100 100 1100 1100]))
set(gcf, 'color', 'w')

%% Plot only pain
figure
u_c = unique(cond_inds_F);

for i_c = 1:numel(u_c)
    if (u_c(i_c)) == 0, continue, end
    ic_idx = ismember(cond_inds_F, (u_c(i_c)));
    hold on
    scatter(zvals(ic_idx,1), zvals(ic_idx,2), 'Marker','.')
    % plot(zvals(:,1),zvals(:,2),'ob','MarkerFaceColor','b', 'MarkerSize',2)


    title({['Granu: ',num2str(analysisparams.tsnegranularity)], ['Perp: ', num2str(perplexity)]})
end
hold off

set(gcf,'Position',([100 100 1100 1100]))
set(gcf, 'color', 'w')

%% do subplots to plot then in 2 plots

figure
u_c = unique(cond_inds_F);

for i_c = 1:numel(u_c)
    if (u_c(i_c)) == 0
        ax = subplot(1,2,1);
    else
        ax = subplot(1,2,2);
    end
    ic_idx = ismember(cond_inds_F, (u_c(i_c)));
    hold on
    scatter(ax, zvals(ic_idx,1), zvals(ic_idx,2), 'Marker','.')
    % plot(zvals(:,1),zvals(:,2),'ob','MarkerFaceColor','b', 'MarkerSize',2)


    title({['Granu: ',num2str(analysisparams.tsnegranularity)], ['Perp: ', num2str(perplexity)]})
end
hold off

set(gcf,'Position',([100 100 1100 1100]))
set(gcf, 'color', 'w')

%% do sal in gray 
figure
u_c = unique(cond_inds_F);
colors_p = {'black', 'red', 'green', 'cyan', 'blue'};
for i_c = 1:numel(u_c)
    if (u_c(i_c)) == 0
       col = colors_p{1};
       alpha = 0.05;
    else
        col = colors_p{i_c};
        alpha = 0.5;
    end
    ic_idx = ismember(cond_inds_F, (u_c(i_c)));
    hold on
    scatter(zvals(ic_idx,1), zvals(ic_idx,2), 1,'Marker','o', ...
        'MarkerEdgeColor',col, 'MarkerFaceAlpha',alpha, 'MarkerEdgeAlpha',alpha)
    % plot(zvals(:,1),zvals(:,2),'ob','MarkerFaceColor','b', 'MarkerSize',2)


end
hold off
title({['Granu: ',num2str(analysisparams.tsnegranularity)], ['Perp: ', num2str(perplexity)]})


set(gcf,'Position',([100 100 1100 1100]))
set(gcf, 'color', 'w')

%% Fo the same but for the saline
figure
u_c = unique(cond_inds_S);
colors_p = {'black', 'red', 'green', 'cyan', 'blue'};
for i_c = 1:numel(u_c)
    if any(ismember([2     4     6     8], (u_c(i_c)) ))
       col = colors_p{1};
       alpha = 0.05;
    else
        col = colors_p{i_c};
        alpha = 0.5;
    end
    ic_idx = ismember(cond_inds_S, (u_c(i_c)));
    hold on
    scatter(zvals(ic_idx,1), zvals(ic_idx,2), 1,'Marker','o', ...
        'MarkerEdgeColor',col, 'MarkerFaceAlpha',alpha, 'MarkerEdgeAlpha',alpha)
    % plot(zvals(:,1),zvals(:,2),'ob','MarkerFaceColor','b', 'MarkerSize',2)


end
hold off
title({['Granu: ',num2str(analysisparams.tsnegranularity)], ['Perp: ', num2str(perplexity)]})


set(gcf,'Position',([100 100 1100 1100]))
set(gcf, 'color', 'w')