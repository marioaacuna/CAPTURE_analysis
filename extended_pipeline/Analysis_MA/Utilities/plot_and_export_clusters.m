function plot_and_export_clusters(analysisstruct_temp, rootpath, ratname, plot_poses)
%PLOT_AND_EXPORT_CLUSTERS Plot t-SNE clusters and cluster poses; export to files.
%   plot_and_export_clusters(analysisstruct_temp, rootpath, ratname, plot_poses)
%
% Ensures minimal fields for downstream plotting and saves:
%   - Tsne_clusters.pdf
%   - Poses_clusters.pdf (if plot_poses true)

arguments
    analysisstruct_temp struct
    rootpath (1,:) char
    ratname (1,:) char = 'rat'
    plot_poses (1,1) logical = true
end

% behavior plots and movies meta
try
    analysisstruct_temp.conditionnames = {ratname};
    analysisstruct_temp.ratnames = {ratname};
catch
    % ignore
end
try
    if isfield(analysisstruct_temp,'mocapstruct_reduced_agg') && ~isempty(analysisstruct_temp.mocapstruct_reduced_agg)
        filesz = size(analysisstruct_temp.mocapstruct_reduced_agg{1}.aligned_mean_position,1);
    else
        filesz = size(analysisstruct_temp.zValues,1);
    end
    analysisstruct_temp.filesizes = {filesz};
catch
    analysisstruct_temp.filesizes = {size(analysisstruct_temp.zValues,1)};
end

% Plot a tsne map
h1 = figure(609); clf;
params = struct();
params.nameplot = 1;
params.density_plot = 0;
params.watershed = 1;
params.sorted = 1;
params.markersize = 1;
params.coarseboundary = 0;
params.do_coarse = 0;
plot_clustercolored_tsne(analysisstruct_temp, 1, params.watershed, h1, params);
set(h1,'Position',([100 100 1100 1100]));

% Export clusters map
if ~exist(rootpath, 'dir'); mkdir(rootpath); end
cluster_figure_filename = fullfile(rootpath,'Tsne_clusters.pdf');
try
    exportgraphics(h1, cluster_figure_filename);
catch
    warning('Failed to export t-SNE clusters figure to %s', cluster_figure_filename);
end

% Plot cluster poses
try
    [cls, ~, ~] = unique(analysisstruct_temp.annot_reordered{end}, 'stable');
catch
    % If annot_reordered not present, skip poses
    cls = [];
end

if plot_poses && ~isempty(cls)
    fig_poses = figure('pos', [10,10,1500,1200]);
    nclus = numel(cls);
    n_rows = ceil(sqrt(nclus));
    n_cols = ceil(sqrt(nclus));
    for ic = 1:numel(cls)
        subplot(n_rows, n_cols, ic);
        this_cls = cls(ic);
        try
            plot_mean_cluster_aligned(analysisstruct_temp.mocapstruct_reduced_agg{1}, ...
                find(analysisstruct_temp.annot_reordered{end}==this_cls), ['cl nr :  ', num2str(this_cls)]);
        catch ME
            title(sprintf('cl %d (plot error: %s)', this_cls, ME.message), 'Interpreter','none');
        end
        title(num2str(this_cls));
    end
    % Export poses
    cluster_poses_figure_filename = fullfile(rootpath, 'Poses_clusters.pdf');
    try
        exportgraphics(fig_poses, cluster_poses_figure_filename);
    catch
        warning('Failed to export poses figure to %s', cluster_poses_figure_filename);
    end
end

end
