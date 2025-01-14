function cluster_vector_ds = downsample_vector(cluster_vector, traces)
    %% Downsample cluster vector to match fluorescence sampling rate
    % Considering missing frames in fluorescence, let's linearly interpolate
    %missing_frames_ratio = length(cluster_vector)/length_traces;
    %ds_factor_adjusted = ds_factor * missing_frames_ratio;
    %cluster_vector_ds = downsample(cluster_vector, round(ds_factor));

    % %%
    % Linearly interpolate the cluster_vector to match the fluorescence traces length
    cluster_vector_interpolated = interp1(1:length(cluster_vector), cluster_vector, linspace(1, length(cluster_vector), length(traces)), 'nearest');


    % cluster_vector_ds = downsample(cluster_vector, ds_factor);
    % Ensure lengths are consistent after downsampling
    min_length = min(length(cluster_vector_interpolated), length(traces));
    cluster_vector_ds = cluster_vector_interpolated(1:min_length);
end
