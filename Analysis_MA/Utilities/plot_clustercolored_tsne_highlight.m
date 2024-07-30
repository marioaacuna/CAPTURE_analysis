function  [h2, colors, clustersused] = plot_clustercolored_tsne_highlight(analysisstruct, iter, wshedflag, h2, params, CLUSTER_ID)
% Plot a tsne map
%Inputs.
%    Analysisstruct - a behavioral analysisstruct
%    Iter: an index to change the figure number
%    wshedflag: a flag indicating whether or not to plot the watershed
%    h2: a figure handle input
%    params: a paramters struct with inputs:
%           params.nameplot: whether to plot cluster names
%           params.density_plot: whether to plot as a density map
%           params.watershed: whether to overlay watershed lines
%           params.sorted:  whether to use the sorted cluster anmes
%           params.markersize: size of markers in the tsne plot
%           params.coarseboundary: whether to plot the coarse cluster
%           boundary, if present
%           params.do_coarse: whether to do coarse



% if plotting points by their coarse names
coarse_names = {'Rearing','LGroom','RGroom','LScratch','RScratch','FaceGroom','PostureAdjust','Walk','WetDogShake','Other'};
coarse_colors = {[1 0.65 0],'r','g','r','g','b','k',[0.5 0 0.5],[0.5 0 0.5],'k'};


if nargin<4
    h2=figure(606+iter)
end
colors = othercolor('Mrainbow',analysisstruct.density_objects);
set(h2,'Color','k')
clustersused = analysisstruct.sorted_clust_ind;
markersize_here = 0.2;
do_coarse =0;
jitter = 0;
coarseboundary = 0;

if nargin<5
    nameplot = 1;
    density_plot = 0;
    sorted =1;
    params = [];
    colors = colors;
else
    nameplot=   params.nameplot ;
    density_plot=params.density_plot ;
    sorted = params.sorted;
    if isfield(params,'markersize')
        markersize_here = params.markersize;
    end
    if isfield(params,'do_coarse')
        do_coarse = params.do_coarse;
    end
    if isfield(params,'jitter')
        jitter = params.jitter;
    end
    if isfield(params,'color')
        colors = params.color;
    end
    if isfield(params,'coarseboundary')
        coarseboundary = params.coarseboundary;
    end
end

if isfield(params,'plotthresh')
    plotthresh = params.plotthresh;
else
    plotthresh=1;
end
%plotthresh = 50;

if ~density_plot
    if ~do_coarse

        highlight_color = [1, 0, 0]; % using red to highlight, but can be changed
        highlight_markersize = 2; % Increased marker size for highlighting

        % Plotting logic
        for ll = 1:analysisstruct.density_objects
            current_cluster = analysisstruct.sorted_clust_ind(ll);
            idx = find(analysisstruct.annot_reordered{end, end} == current_cluster);

            if numel(idx) > plotthresh
                if any(ismember(CLUSTER_ID , current_cluster))
                    %current_cluster== CLUSTER_ID
                    % Highlight the specific cluster
                    plot(analysisstruct.zValues(idx, 1), ...
                        analysisstruct.zValues(idx, 2), ...
                        'o', 'MarkerSize', highlight_markersize, ...
                        'Color', colors(ll, :), 'MarkerFaceColor', colors(ll, :), ...
                        'MarkerEdgeColor', highlight_color);
                else
                    plot(analysisstruct.zValues(idx, 1), ...
                        analysisstruct.zValues(idx, 2), ...
                        'o', 'MarkerSize', markersize_here, ...
                        'Color', colors(ll, :), 'MarkerFaceColor', colors(ll, :));
                end
                hold on;
            end
        end
        hold off;

        xlabel('Tsne 1')
        ylabel('Tsne 2')
        set(gca,'fontsize',18)
        box off
    else
        for ll = (1:analysisstruct.density_objects)
            if numel(find(analysisstruct.annot_reordered{end,end}==analysisstruct.sorted_clust_ind(ll)))>plotthresh
                plot(analysisstruct.zValues(find(analysisstruct.annot_reordered{end,end}==analysisstruct.sorted_clust_ind(ll)),1),...
                    analysisstruct.zValues(find(analysisstruct.annot_reordered{end,end}==analysisstruct.sorted_clust_ind(ll)),2),...
                    'o','MarkerSize',markersize_here,'Color',analysisstruct.coarse_cluster_color{analysisstruct.sorted_clust_ind(ll)},...
                    'MarkerFaceColor',analysisstruct.coarse_cluster_color{analysisstruct.sorted_clust_ind(ll)})
                hold on
            end
        end

    end




    if coarseboundary
        figure(607)

        hold on
        B = analysisstruct.coarse_borders
        for kk = 1:numel(B)
            if numel(B{kk})
                plot(analysisstruct.xx((B{kk}{1}(:,2))),...
                    analysisstruct.yy((B{kk}{1}(:,1))),...
                    'Color',[1 1 1],'linewidth',0.01)
            end
        end
        hold off
    end

    if (wshedflag)
        % If you want to display the cluster number
        if wshedflag
            if length(CLUSTER_ID) > 1
                s2 = regionprops(analysisstruct.sorted_watershed, 'Centroid');
                for cli = 1: length(CLUSTER_ID)
                    
                    c = s2(CLUSTER_ID(cli)).Centroid;
                    text(analysisstruct.xx(floor(c(1))), analysisstruct.yy(floor(c(2))), num2str(CLUSTER_ID(cli)), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', 'Color', 'y', 'FontWeight', 'Bold', 'FontSize', 25);
                end
            else

    
                s2 = regionprops(analysisstruct.sorted_watershed, 'Centroid');
                c = s2(CLUSTER_ID).Centroid;
                text(analysisstruct.xx(floor(c(1))), analysisstruct.yy(floor(c(2))), num2str(CLUSTER_ID), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'Color', 'y', 'FontWeight', 'Bold', 'FontSize', 50);
            end
        end


        %% plot names as well
        if nameplot
            hold on
            s2 = regionprops(analysisstruct.sorted_watershed, 'Centroid');
            for k = 1:numel(analysisstruct.sorted_clust_ind')
                % get the ind of the sorted cluster
                if sorted
                    khere = k;%analysisstruct.sorted_clust_ind(k);%analysisstruct.sorted_clust_ind(k);%(analysisstruct.sorted_clust_ind(k));%k;
                else
                    khere = find(analysisstruct.sorted_clust_ind==k);%
                end
                c = s2((khere)).Centroid;

                if numel(find(analysisstruct.annot_reordered{end,end}==khere))>plotthresh
                    text(analysisstruct.xx(floor(c(1))), analysisstruct.yy(floor(c(2))), num2str((k)), ... %sprintf('%d', integer(ind)),
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle','Color','w','FontWeight','Bold');
                end
            end
            hold off
        end
        
        nnn = analysisstruct.sorted_watershed;
        nnn(nnn>0) = 1;
        B = bwboundaries(((nnn)));

        % Define the target cluster ID to highlight
        % CLUSTER_ID = ...;  % Replace ... with the actual cluster ID you want to highlight

        %figure(333)
        hold on
        for kk = 1:numel(B)
            if numel(find(ismember(analysisstruct.sorted_clust_ind,kk)))
                if numel(find(analysisstruct.annot_reordered{end,end}==find(analysisstruct.sorted_clust_ind==kk)))>plotthresh
                    % Default boundary color
                    boundary_color = 'w';
                    boundary_thickness = 1;

                    % Check if the current boundary is for the target CLUSTER_ID
                    if find(analysisstruct.sorted_clust_ind==kk) == CLUSTER_ID
                        boundary_color = 'r';  % Set the highlight color
                        boundary_thickness = 5;
                    end

                    % Plot the boundary with the selected color
                    plot(analysisstruct.xx(B{kk}(:,2)),analysisstruct.yy(B{kk}(:,1)),boundary_color,'LineWidth' ,boundary_thickness)
                end
            end
        end
        hold off
        % nnn = analysisstruct.sorted_watershed;
        % nnn(nnn>0) = 1;
        % B = bwboundaries(((nnn)));
        % %figure(333)
        % hold on
        % for kk = 1:numel(B)
        %     if numel(find(ismember(analysisstruct.sorted_clust_ind,kk)))
        %         if numel(find(analysisstruct.annot_reordered{end,end}==find(analysisstruct.sorted_clust_ind==kk)))>plotthresh
        %             plot(analysisstruct.xx(B{kk}(:,2)),analysisstruct.yy(B{kk}(:,1)),'w')
        %         end
        %     end
        % end
        % hold off

    end
end
%% also plot colored density
% zvals_cell_array = cell(1,analysisstruct.density_objects);
% for ll = 1:analysisstruct.density_objects
%     zvals_cell_array{ll} = analysisstruct.zValues(find(analysisstruct.annotation_vec{end,end}==analysisstruct.sorted_clust_ind(ll)),:);
% end

if ~jitter
    zvals_cell_array = cell(1,max(unique(analysisstruct.condition_inds)));
    for ll = unique(analysisstruct.condition_inds)
        zvals_cell_array{ll} = analysisstruct.zValues(find(analysisstruct.condition_inds==ll),:);
    end
else
    zvals_cell_array = cell(1,max(unique(analysisstruct.condition_inds)));
    for ll = unique(analysisstruct.condition_inds)'
        zvals_cell_array{ll} = analysisstruct.zValues_jitter(find(analysisstruct.condition_inds==ll),:);
    end
end


if density_plot
    badclust = find(cellfun(@numel,strfind(analysisstruct.clusternames,'BadTracking')));
    badframes = find(ismember(analysisstruct.annot_reordered{end,end},badclust));
    goodframes = setxor(1:numel(analysisstruct.annot_reordered{end,end}),badframes);
    if numel(zvals_cell_array) == 1 && numel(badclust)>1
        zvals_cell_array{1} = zvals_cell_array{1}(goodframes,:);
    end
    %h3 = subplot(1,1,1)
    fighand_in = h2%figure(444)
    set(h2,'Color','w')
    plotdensitymaps({cat(1,zvals_cell_array{:})},1,fighand_in,analysisstruct.params.density_width,...
        max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor,analysisstruct.params.density_res)

    figure(608)
    if (wshedflag)
        nnn = analysisstruct.sorted_watershed;
        nnn(nnn>0) = 1;
        B = bwboundaries((flipud(nnn)));
        %figure(333)
        hold on
        for kk = 1:numel(B)
            if numel(find(analysisstruct.annot_reordered{end,end}==find(analysisstruct.sorted_clust_ind==kk)))>plotthresh
                plot(B{kk}(:,2),B{kk}(:,1),'Color',[0.95 0.95 0.95])
            end
        end
        hold off
    end


    if coarseboundary
        figure(607)

        hold on
        B = analysisstruct.coarse_borders
        for kk = 1:numel(B)
            %  plot((B{kk}{1}(:,2)),analysisstruct.params.density_res-(B{kk}{1}(:,1)),'Color',[0.95 0.95 0.95],'linewidth',0.01)
            plot((B{kk}{1}(:,2)),analysisstruct.params.density_res-(B{kk}{1}(:,1)),'Color',[0 0 0],'linewidth',0.01)

        end
        hold off

        figure(608)

        hold on
        B = analysisstruct.coarse_borders
        for kk = 1:numel(B)
            %  plot((B{kk}{1}(:,2)),analysisstruct.params.density_res-(B{kk}{1}(:,1)),'Color',[0.95 0.95 0.95],'linewidth',0.01)
            plot((B{kk}{1}(:,2)),analysisstruct.params.density_res-(B{kk}{1}(:,1)),'Color',[0 0 0],'linewidth',0.01)

        end
        hold off
    end

    %% plot names as well
    if (wshedflag)

        hold on
        s2 = regionprops(analysisstruct.sorted_watershed, 'Centroid');
        for k = 1:numel(analysisstruct.sorted_clust_ind')
            %     % get the ind of the sorted cluster
            %     c = s2((k)).Centroid;
            %     text(analysisstruct.xx(floor(c(1))), analysisstruct.yy(floor(c(2))), num2str((k)), ... %sprintf('%d', integer(ind)),
            %         'HorizontalAlignment', 'center', ...
            %         'VerticalAlignment', 'middle','Color','Red','FontWeight','Bold');
        end
        hold off
    end
end
ax = gca;
ax.Color = [0, 0, 0];
axis off
%if numel(legendinput)
%legend(legendinput)
%end
end