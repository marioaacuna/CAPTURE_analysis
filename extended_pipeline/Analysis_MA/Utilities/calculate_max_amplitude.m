function [max_amplitudes, unique_clusters] = calculate_max_amplitude(traces_interpolated, cluster_vector_ds)
    %% Collect max amplitude of all cells in this session for all clusters
    unique_clusters = unique(cluster_vector_ds);
    n_clusters = length(unique_clusters);
    n_rois = size(traces_interpolated, 1);
    max_amplitudes = zeros(n_rois, n_clusters);

    for i = 1:n_clusters
        current_cluster = unique_clusters(i);
        idx = cluster_vector_ds == current_cluster; % indices where current cluster is present
        for j = 1:n_rois
            roi_trace = traces_interpolated(j, idx);
            max_roi = max(roi_trace);
            if max_roi<0
                max_roi = 0;
            end
            max_amplitudes(j, i) =max_roi;
        end
    end
end