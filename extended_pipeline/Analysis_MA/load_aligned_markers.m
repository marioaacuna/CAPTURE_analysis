
function markers_aligned_ds = load_aligned_markers(markers_aligned_preproc, repfactor, smooth_factor)
    markers = fieldnames(markers_aligned_preproc);
    markers_aligned_ds = struct();
    for i = 1:length(markers)
        marker_name = markers{i};
        marker_data = markers_aligned_preproc.(marker_name);
        % Downsample
        marker_data_ds = marker_data(1:repfactor:end, :);
        % Smooth data to remove high-frequency noise
        marker_data_smooth = smoothdata(marker_data_ds, 'sgolay', smooth_factor);
        markers_aligned_ds.(marker_name) = marker_data_smooth;
    end
end