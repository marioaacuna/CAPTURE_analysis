%% Create GIF for Significant Clusters
% Identify significantly higher clusters in F condition
sig_clusters = to_take;

% Set up the figure for the GIF
n_sig_clusters = length(sig_clusters);
n_rows = ceil(sqrt(n_sig_clusters));
n_cols = ceil(n_sig_clusters / n_rows);

fig = figure('Position', [100, 100, 1200, 900], 'Visible','off');

% Create the GIF
filename = fullfile(rootpath, 'significant_clusters_animation.gif');
frame_duration = 0.1; % Duration of each frame in seconds

for frame = 1:1000 % Adjust the number of frames as needed
    clf; % Clear the figure for each new frame
    
    for ic = 1:n_sig_clusters
        subplot(n_rows, n_cols, ic);
        this_cls = sig_clusters(ic);
        
        % Find frames for this cluster
        cluster_frames = find(analysisstruct.annot_reordered{end} == this_cls);
        
        % Randomly select one frame from this cluster for this animation frame
        if ~isempty(cluster_frames)
            random_frame = cluster_frames(randi(length(cluster_frames)));
            
            % Animate the markers for this frame
            animate_markers_aligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1}, random_frame, gcf, ['Cluster ', num2str(this_cls)]);
        end
        
        title(['Cluster ', num2str(this_cls)]);
    end
    
    % Capture the plot as an image
    frame_img = getframe(fig);
    im = frame2im(frame_img);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to the GIF File
    if frame == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', frame_duration);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_duration);
    end
end

close(fig);
disp(['GIF saved as: ', filename]);
%%