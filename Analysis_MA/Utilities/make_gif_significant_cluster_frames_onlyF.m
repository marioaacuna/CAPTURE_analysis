
%% Create GIF for Significant Clusters (F condition only)
% Identify significantly higher clusters in F condition
% sig_clusters = [];
% for i = 1:length(clusters)
%     cluster_name = ['id_', num2str(clusters(i))];
%     if p_values.(cluster_name) < 0.05 && mean(cluster_proportions_F(:,i)) > mean(cluster_proportions_S(:,i))
%         sig_clusters = [sig_clusters, clusters(i)];
%     end
% end

sig_clusters =predominantS;

% Set up the figure for the GIF
n_sig_clusters = length(sig_clusters);
n_rows = ceil(sqrt(n_sig_clusters));
n_cols = ceil(n_sig_clusters / n_rows);

fig = figure('Position', [100, 100, 1200, 900], Visible='on');

% Create the GIF
filename = fullfile(rootpath, 'predominantS_clusters_animation.gif');
frame_duration = 0.1; % Duration of each frame in seconds

% Get F condition frames
F_condition_frames = contains(upsampled_identifiers(good_frames), '_F');

for frame = 1:250 % Adjust the number of frames as needed
    clf; % Clear the figure for each new frame
    
    for ic = 1:n_sig_clusters
        subplot(n_rows, n_cols, ic);
        this_cls = sig_clusters(ic);
        
        % Find frames for this cluster that are also in F condition
        cluster_frames = (find(analysisstruct.annot_reordered{end} == this_cls & logical(F_condition_frames)'));
        
        % Randomly select one frame from this cluster for this animation frame
        if ~isempty(cluster_frames)
            random_frame = cluster_frames(randi(length(cluster_frames)));
            
            % Get the actual frame number from good_frames
            actual_frame = random_frame;
            
            % Animate the markers for this frame
            animate_markers_aligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1}, actual_frame, gcf, ['Cluster ', num2str(this_cls)]);
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