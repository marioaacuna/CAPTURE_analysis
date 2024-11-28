function plot_identified_clusters(clusters, analysisstruct)
    
    to_take = clusters;
    fig_predominant = figure('pos', [10,300,1500,1900]);
    n_rows = ceil(sqrt(numel(to_take)));
    n_cols = ceil(sqrt(numel(to_take)));
    
    for ic = 1:numel(to_take)
        subplot(n_rows, n_cols, ic)
        this_cls = to_take(ic);
        fprintf('ic = %i - \n', this_cls)
        plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1},...
            find(analysisstruct.annot_reordered{end}==this_cls),['cl nr :  ', num2str(this_cls)]);
        title(this_cls)
    end
    
end