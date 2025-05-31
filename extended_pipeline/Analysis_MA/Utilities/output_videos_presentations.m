%% Preamble
% This script outputs representative videos of predictions and ms video for
% a particular animal and a particular session. So far, this is only for
% showing-off reasons, maybe putting them in a presentation and so. 
% 
% UPDATED VERSION: Now supports multiple conditions (baseline and Formalin_injection)
% and uses flexible path structure for data from multiple sources.
%
% Variables:
%   - animal_ID: Animal identifier (e.g., 'ID_1386')
%   - condition: 'baseline' or 'Formalin_injection'
%   - Data sources:
%     * Behavior videos: K:\Mario\BioMed_students_2023\Anna\exp_6cam_miniscope_data\[condition]\6cam_data\[animal]\[date]\videos\Camera1\0.mp4
%     * FOV recordings: K:\Ca_imaging_pain\2_motion_corrected_movies\[animal]\[animal]_[date]_SP.h5
%     * Raw miniscope: K:\Mario\BioMed_students_2023\Anna\exp_6cam_miniscope_data\[condition]\miniscope_data\[animal]\[date_miniscope]\msCam0.avi
%     * Predictions: [behavior_path]\DANNCE\predict_results\predictions.mat
%     * Clusters: From preprocessing pipeline clusters_struct.mat
%     * Analysis struct: From general_configs filename_analysis
%
% The basis of this analysis is the use of cluster ids and output only
% those frames that belong to a cluster. We need to output the 3D predictions,
% and the neuronal activity (FOV of active neurons)
% They all have different frame rates, however, they all have the same duration
% (in this case 30 min). 
% Maybe it would be good to downsample everything to the lowest (FOV recording)
% Specific cases:
%   FOV: this is a h5 file, that could be too large to fit in memory at once, so
%       it needs to be read in chunks. Or particularly, we need to read only those
%       important frames, ie. where we see a cluster id (>0). Frame rate 5 Hz
%   Cluster id vector: this is a matfile (or csv) containing the vector of cluster
%       ids. Frame rate: 120Hz
%   Predictions: this is a matfile structure that contains cells representing each
%       body part. Frame rate 120. This is not centered, just raw.
%   behavior: This is a mp4 file from a single camera. Frame rate 120 Hz.
%   Raw miniscope: This is an avi file with raw calcium imaging data. Frame rate ~20 Hz.


%%
clear
clc
close all
%%
TEST_clusters = 0; % if debugging to test cluster ids
%%
global FR_to_downsample

% Add general configs
GC = general_configs();

% Configuration
animal_to_take = 'ID_1386';
condition = 'Formalin_injection'; % 'baseline' or 'Formalin_injection'

% Define condition-specific parameters
condition_info = struct();
condition_info.baseline.folder = 'baseline';
condition_info.baseline.date = '240918';
condition_info.baseline.date_miniscope = '20240918';
condition_info.baseline.cluster_suffix = 'B';

condition_info.Formalin_injection.folder = 'Formalin_injection';
condition_info.Formalin_injection.date = '241030';
condition_info.Formalin_injection.date_miniscope = '241030';
condition_info.Formalin_injection.cluster_suffix = 'F';

% Root paths
behavior_rootpath = 'K:\Mario\BioMed_students_2023\Anna\exp_6cam_miniscope_data';
FOV_rootpath = 'K:\Ca_imaging_pain\2_motion_corrected_movies';
miniscope_raw_rootpath = 'K:\Mario\BioMed_students_2023\Anna\exp_6cam_miniscope_data';

% Get condition-specific info
curr_condition = condition_info.(condition);

% Build paths based on condition
sixcam_path = fullfile(behavior_rootpath, curr_condition.folder, '6cam_data', animal_to_take, curr_condition.date);
behavior_file = fullfile(sixcam_path, 'videos', 'Camera1', '0.mp4');
predictions_file = fullfile(sixcam_path, 'DANNCE', 'predict_results', 'predictions.mat');

% FOV file path
FOV_file = fullfile(FOV_rootpath, animal_to_take, [animal_to_take, '_', curr_condition.date, '_SP.h5']);

% Alternative raw miniscope file
raw_miniscope_file = fullfile(miniscope_raw_rootpath, curr_condition.folder, 'miniscope_data', animal_to_take, curr_condition.date_miniscope, 'msCam0.avi');

% Cluster file from preprocessing
clusters_struct_file = fullfile(GC.preprocessing_rootpath, 'clusters_struct.mat');
cluster_field_name = [animal_to_take, '_', curr_condition.cluster_suffix];

% Analysis struct file
analysisstruct_filename = GC.filename_analysis;

% Load files
h5_info = h5info(FOV_file);
dataset = '/1'; %dataset of h5 file
%read here h5 file
h5_data = h5read(FOV_file, dataset);

% Load alternative raw miniscope file
if exist(raw_miniscope_file, 'file')
    raw_miniscope_video = VideoReader(raw_miniscope_file);
    raw_miniscope_FR = raw_miniscope_video.FrameRate;
    fprintf('Raw miniscope video found with frame rate: %.2f Hz\n', raw_miniscope_FR);
else
    warning('Raw miniscope file not found: %s', raw_miniscope_file);
    raw_miniscope_video = [];
end

% read clusters vector from clusters_struct
if exist(clusters_struct_file, 'file')
    clusters_data = load(clusters_struct_file);
    if isfield(clusters_data.clusters_struct, cluster_field_name)
        clusters = clusters_data.clusters_struct.(cluster_field_name);
        fprintf('Loaded clusters for %s\n', cluster_field_name);
    else
        error('Cluster field %s not found in clusters_struct', cluster_field_name);
    end
else
    error('Clusters struct file not found: %s', clusters_struct_file);
end

% read predictions
pred = load(predictions_file);
% pred = pred.predictions;
predictions = pred.predictions;

% read behavioral movie
beh_movie = VideoReader(behavior_file);

% Load analysis struct
load(analysisstruct_filename)

% upsample FOV h5 file 120 hz OR downsample everything to 5Hz
% to make coincidence between all the variables
% n_frames = h5_info.Datasets.Dataspace.Size(3);
% FR_FOV = 1/((30*60)/n_frames)
% FR_FOV = 5;
% FR_pred = 120;

% If you decide to downsample everything to 5Hz:
% clusters = resample(clusters, 1, FR_pred/FR_FOV);
% predictions = structfun(@(x) resample(x, 1, FR_pred/FR_FOV), predictions, 'UniformOutput', false);

% Calculate the downsampling rate
% downsampling_rate = 120/5; % which is 24

% Identify frames where cluster IDs > 0
% Linearly interpolate the cluster_vector to match the fluorescence traces length

%% NEW
% frames_with_clusters = find(cluster_vector_interpolated > 0);
% 
% 
% % Preallocate the upsampled frames array with zeros
% frames_with_clusters_upsampled = zeros(1, beh_movie.NumFrames);
% 
% % Populate the upsampled frames array based on the frames_with_clusters
% for frame = frames_with_clusters
%     start_idx = (frame-1) * downsampling_rate + 1;
%     end_idx = frame * downsampling_rate;
%     frames_with_clusters_upsampled(start_idx:end_idx) = 1;
% end
% 
% 
% 
% 
% % Downsample those frames: upsample to cover the total frames of the movie
% frames_to_read = find(frames_with_clusters_upsampled);
% frames_to_read = downsample(frames_to_read, downsampling_rate);
% % Clusters downsampled
% % clusters_frames = clusters(frames_with_clusters);
% clusters_ds = cluster_vector_interpolated(cluster_vector_interpolated >0);
% 
% % Calculate equivalent frame indices for 5 Hz data
% frames_ro_read_FOV = frames_with_clusters;




%% Original thought
% cluster_vector_interpolated = interp1(1:length(clusters), clusters, linspace(1, length(clusters), size(h5_data, 3)), 'nearest');
% 
% frames_with_clusters = find(clusters > 0);
% 
% % Downsample those frames
% frames_to_read = frames_with_clusters(1:downsampling_rate:end);
% 
% % Clusters downsampled
% clusters_frames = clusters(frames_with_clusters);
% clusters_ds = clusters_frames(1:downsampling_rate:end);
% 
% % Calculate equivalent frame indices for 5 Hz data
% frames_ro_read_FOV = find(cluster_vector_interpolated>0);
% % frames_ro_read_FOV = ceil(frames_to_read / downsampling_rate);

%% NEW 2
FR_video =   round(beh_movie.FrameRate); 
FR_FOV = 5;
length_FOV = size(h5_data, 3);
length_movie = length(pred.predictions.AnkleL); % this is same as predictions
% Calculate the downsampling rate
downsampling_rate = FR_video/FR_FOV; % which is 24

% Align clusters to the length of behavior video and predictions
cluster_vector_interpolated = interp1(1:length(clusters), clusters, linspace(1, length(clusters), length_movie), 'nearest');

% Identify frames where cluster IDs > 0
frames_with_clusters = find(cluster_vector_interpolated > 0);

% Downsample those frames to be used for beh movie and predictions
frames_to_read = frames_with_clusters(1:downsampling_rate:end);

% Get the downsampled clusters
clusters_ds = cluster_vector_interpolated(frames_to_read);

% Calculate equivalent frame indices for 5 Hz data (FOV)
frames_to_read_FOV = ceil(frames_to_read / downsampling_rate);

% Check to ensure we don't request more frames than available in FOV
frames_to_read_FOV = frames_to_read_FOV(frames_to_read_FOV <= length_FOV);

% Assume frames_to_read_FOV is our reference length
reference_length = length(frames_to_read_FOV);

% Interpolate frames_to_read to match reference length
frames_to_read_interp = round(interp1(1:length(frames_to_read), frames_to_read, linspace(1, length(frames_to_read), reference_length)));
frames_to_read = frames_to_read_interp;

% Interpolate cluster_vector_interpolated to match reference length
clusters_ds_interp = interp1(1:length(clusters_ds), clusters_ds, linspace(1, length(clusters_ds), reference_length), 'nearest');
clusters_ds = clusters_ds_interp;

% At this point, the lengths should match
assert(length(frames_to_read_interp) == length(frames_to_read_FOV));
assert(length(clusters_ds_interp) == length(frames_to_read_FOV));
% % Adjust all vectors to be of the same length (by trimming to the shortest length)
% min_length = min([length(frames_to_read), length(frames_to_read_FOV), length(clusters_ds)]);
% 
% frames_to_read = frames_to_read(1:min_length);
% frames_to_read_FOV = frames_to_read_FOV(1:min_length);
% clusters_ds = clusters_ds(1:min_length);

%%
% Read ROIs - commented out for now as we're not using the mask
% TODO, load ROI info from animal
% ROIs
% Transform ROIS for mask
% Assuming your image size is MxN
FOV_size = size(h5_data);
FOV_size = FOV_size([1,2]);
M = FOV_size(1); 
N = FOV_size(2);

% Initialize the mask to zeros (no mask for now)
mask = ones(M, N); % Use ones to not mask anything

% % Define the expansion size (in pixels)
% expand_size = 5;
% 
% % Create a disk-shaped structuring element for dilation
% SE = strel('disk', expand_size);
% 
% % Loop over each ROI
% for i = 1:size(ROI_info, 1)
%     x = ROI_info{i, 5}(:, 1);
%     y = ROI_info{i, 5}(:, 2);
%     
%     % Method 1: Expand the polygon and then create the mask
%     roi_mask = poly2mask(x, y, M, N);
%     roi_mask_dilated = imdilate(roi_mask, SE);
%     mask = mask | roi_mask_dilated;
% 
%     % Method 2: Create a round shape from the centroid
%     % centroid = mean([x, y]);
%     % circle_mask = insertShape(zeros(M, N, 'uint8'), 'FilledCircle', [centroid, expand_size*2], 'Color', 'white', 'Opacity', 1);
%     % mask = mask | circle_mask;
% end
% process h5 file
processed_frames = processFOVFrames(h5_data);


%% TEST specific cluster 
if TEST_clusters
    this_camera = pred.cameras.(['Camera', num2str(1)]);
    K = this_camera.IntrinsicMatrix;
    R =  this_camera.rotationMatrix;
    t = this_camera.translationVector;
    M = [R; t] * K;
    K1 = K';

    vid1 = beh_movie;
    camparams = cameraParameters(K = K1, RadialDistortion = this_camera.RadialDistortion,...
        TangentialDistortion= this_camera.TangentialDistortion);


    Fig_vids = figure('Position',[20 20 1800 850], 'Visible','on', 'color', 'k'); % TODO: check position and size of the figure
    unique_cls = unique(clusters_ds, 'stable');
    % write video
    output = fullfile(output_base_dir, 'test_rearing.mp4'); % cls 23, 45, 47
 
    writerObj = VideoWriter(fullfile(output), 'MPEG-4');
    writerObj.Quality = 100;
    writerObj.FrameRate = FR_to_downsample; % TODO: this needs to be done
    open(writerObj);
    this_cl = [28];
    % for ic = 1:length(unique_cls)
        % this_cl = unique_cls(ic);
        frames_to_choose = frames_to_read(ismember(clusters_ds, this_cl));
        f = frames_to_choose;
        for iframe = 1:length(f)
            this_frame_idx = f(iframe) ;
            this_frame = read(vid1, this_frame_idx);

            h1 = imagesc(undistortImage(this_frame, camparams));
            title(num2str(this_cl))
            writeVideo(writerObj, getframe(Fig_vids));
            drawnow
        end
    % end
    close(writerObj);
end

%% Run outputs
% Create output directory structure
output_base_dir = fullfile(GC.temp_root, 'output_videos', animal_to_take, condition);
if ~exist(output_base_dir, 'dir')
    mkdir(output_base_dir);
end

output_vid_beh = fullfile(output_base_dir, 'beh.mp4');
output_vid_FVO = fullfile(output_base_dir, 'FOV.mp4');
output_vid_pred = fullfile(output_base_dir, 'pred.mp4');
output_vid_clusters = fullfile(output_base_dir, 'clusters.mp4');

% Also create outputs for raw miniscope if available
if ~isempty(raw_miniscope_video)
    output_vid_raw_miniscope = fullfile(output_base_dir, 'raw_miniscope.mp4');
end

FR_to_downsample = 5; % 5Hz

% pred.predictions = predictions;
run_beh_video(beh_movie, frames_to_read, output_vid_beh, pred)
run_predictions(pred, frames_to_read, output_vid_pred, clusters_ds)
run_FOV(processed_frames, frames_to_read_FOV, output_vid_FVO, mask)
run_cluster_highlights(analysisstruct, frames_to_read, output_vid_clusters, clusters_ds)

% Process raw miniscope video if available
if ~isempty(raw_miniscope_video)
    % Calculate frames for raw miniscope video
    raw_miniscope_length = raw_miniscope_video.NumFrames;
    FOV_length = size(h5_data, 3);
    
    % Map FOV frames to raw miniscope frames
    frames_to_read_raw = round(frames_to_read_FOV * (raw_miniscope_length / FOV_length));
    frames_to_read_raw = frames_to_read_raw(frames_to_read_raw <= raw_miniscope_length & frames_to_read_raw > 0);
    
    run_raw_miniscope_video(raw_miniscope_video, frames_to_read_raw, output_vid_raw_miniscope, clusters_ds(1:length(frames_to_read_raw)));
end

fprintf('Video processing completed for %s condition %s\n', animal_to_take, condition);
fprintf('Output videos saved to: %s\n', output_base_dir);



%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions
% These functions will make the videos, all videos need to have the same size and the same frame rate.
% these videos are later going to be used for making a professional video in Adobe Premier Pro

% function run_FOV(m, f, output)
    
    % this function will get the frames and then output a mp4 video of the FOV according to the clusters at 5 Hz
    % n_frames = h5_info.Datasets.Dataspace.Size(3);
    % TODO: fill this up
%    Fig_vids = figure('Position',[20 20 1800 850], 'Visible',visibility); % TODO: check position and size of the figure. Needs to match all figures(videos)
%    writerObj = VideoWriter(fullfile(output), 'MPEG-4');
%    writerObj.Quality = 100;
%    writerObj.FrameRate = FR_to_downsample; % TODO: this needs to be done
    
    % I'll add just a structure, but it's not clear what needs to be done inside:
%    open(writerObj);
%    for frameIdx = f
        % extract frame data from h5_data using frameIdx
%        frameData = h5_data(:,:,frameIdx); 
        % visualize or process frameData
        % write to the video
        
%        writeVideo(writerObj, getframe(Fig_vids));
%    end
%    close(writerObj);



% end

function run_FOV(m, f, output, mask)
    global FR_to_downsample
    writerObj = VideoWriter(fullfile(output), 'MPEG-4');
    writerObj.Quality = 100;
    writerObj.FrameRate = FR_to_downsample; % TODO: this needs to be done

    % Initialize the video writer
    open(writerObj);
    
    % Create the figure for visualization
    Fig_vids = figure('Position',[20 20 1800 850], 'Visible','off', 'color', 'k'); % 'off' ensures the figure doesn't pop up during the loop
    
    % Iterate through the frames
    for frameIdx = f
        % Extract frame data from h5_data
        frameData = m(:,:,frameIdx); 
        segmented_image = frameData .* (mask);
        % control for over saturation
        mean_seg = mean(segmented_image(segmented_image>0));
        min_lim = 0.45;
        if mean_seg > 0.4
            % keyboard
            % segmented_image=  segmented_image .* 0.75;
            segmented_image = mask .* 0;
            segmented_image = segmented_image - mean_seg;
            min_lim = mean_seg + mean_seg*0.15;
            if min_lim >0.85
                min_lim = 0.85;
            end
        end
        % Visualize the data (assuming it's 2D image data)
        imagesc(segmented_image);
        colormap('gray'); % choose an appropriate colormap
        clim([min_lim, 1])
        axis off; % to turn off axis numbering and ticks
        
        % Optional: You can add a title or any annotations if needed
        % title(['Frame: ' num2str(frameIdx)]);
        
        % Save this visualization to the video
        writeVideo(writerObj, getframe(Fig_vids));
        % drawnow
    end
    
    % Close the video writer
    close(writerObj);
    
    % Close the visualization figure
    close(Fig_vids);
end


function run_predictions(preds, frames, output, cls)
    global FR_to_downsample
    % this function will get the frames and then output a mp4 video of the clusters at 5 Hz
    this_camera = preds.cameras.(['Camera', num2str(1)]);
    K = this_camera.IntrinsicMatrix;
    R =  this_camera.rotationMatrix;
    t = this_camera.translationVector;
    M = [R; t] * K;
    K1 = K';
    % TODO: init video, with output variable
 
    Fig_vids = figure('Position',[20 20 1800 850], 'Visible','off', 'color', 'k'); % TODO: check position and size of the figure. Needs to match all figures(videos)
    % ax = axes;
    writerObj = VideoWriter(fullfile(output), 'MPEG-4');
    writerObj.Quality = 100;
    writerObj.FrameRate = FR_to_downsample; % TODO: this needs to be done
    open(writerObj);
    sklton = fullfile("H:\Cam_cal_and_label3D_Nevian\Label3D\skeletons\mouse22.mat");
    skeleton = load(sklton);
    
    names = fieldnames(preds.predictions);
    names(ismember(names, 'sampleID')) = [];
    handles_here = cell(1,numel(names));
    n_frames = length(frames);
    cls_idx = 1;
    for iframe = frames% hasFrame(vid) && hasFrame(ani)
        % iid = iid +1;
        % plot real
        clf
            
        ind_to_plot = iframe;
        pts_this_frame = NaN(numel(names),3);
        %     title(ax, str_title, 'Color','w')
        for jj = 1:numel(names)
            % don't plot markers that drop out
            % if ~isnan(sum(preds.predictions.(names{jj})(ind_to_plot,:),2))
            if (~sum(preds.predictions.(names{jj})(ind_to_plot,:),2) == 0)
                xx = squeeze(preds.predictions.(names{jj})(ind_to_plot,1));
                yy = squeeze(preds.predictions.(names{jj})(ind_to_plot,2));
                zz = squeeze(preds.predictions.(names{jj})(ind_to_plot,3));
                % handles_here{jj} = line(xx,yy,zz,'Marker','o','Color',skeleton.color(jj,:),'MarkerFaceColor',skeleton.color(jj,:),'MarkerSize',5);
                pts_this_frame(jj,:) = [xx,yy,zz];            
                hold on
                marker_plot(jj) = 1;
            else
                marker_plot(jj) = 0;
            end

            % end
        end

        pts = pts_this_frame;
        projPts = [pts, ones(size(pts, 1), 1)] * M;
        projPts(:, 1:2) = projPts(:, 1:2) ./ projPts(:, 3);
        scatter(projPts(:,1), projPts(:,2), 'ro', 'filled', 'Marker','o')
        hold on
        %% plot the links between markers
        links = skeleton.joints_idx;
        colors = skeleton.color;
        n_links = length(links);
        for mm = 1:(n_links)
           
            xx = [squeeze(preds.predictions.(names{links(mm,1)})(ind_to_plot,1)) ...
                squeeze(preds.predictions.(names{links(mm,2)})(ind_to_plot,1)) ];
            yy = [squeeze(preds.predictions.(names{links(mm,1)})(ind_to_plot,2)) ...
                squeeze(preds.predictions.(names{links(mm,2)})(ind_to_plot,2)) ];
            zz = [squeeze(preds.predictions.(names{links(mm,1)})(ind_to_plot,3)) ...
                squeeze(preds.predictions.(names{links(mm,2)})(ind_to_plot,3)) ];

            % x
            pts = [xx(1), yy(1), zz(1)];      

            projPts = [pts, ones(size(pts, 1), 1)] * M;
            projPts(:, 1:2) = projPts(:, 1:2) ./ projPts(:, 3);

            pts2  =[xx(2), yy(2), zz(2)];
            projPts2 = [pts2, ones(size(pts2, 1), 1)] * M;
            projPts2(:, 1:2) = projPts2(:, 1:2) ./ projPts2(:, 3);
            
            xx = [projPts(1), projPts2(1)];
            yy = [projPts(2), projPts2(2)];
            %zz = [projPts(3), projPts2(3)];
            this_color = colors(mm, 1:3);
            line(xx,yy,'Color',this_color,'LineWidth',1);
            
        
        end
        drawnow
        xlim([1,1800])
        ylim([1,1000])
        ax = gca;
        ax.YDir = "reverse";
        % ax = axes;
        ax.Color = [0, 0, 0];
        axis off
        % this_frame_idx = frames(iframe);
        title(['Cls - ', num2str(cls(cls_idx))], 'Color','w')
        hold off
        writeVideo(writerObj,getframe(Fig_vids))
        cla
        cls_idx = cls_idx + 1;
    end
    close(writerObj);
    

end

function run_beh_video(m, f, output, preds)
    global FR_to_downsample
    % this function will get the frames and then output a mp4 video of the behavioral video according to the clusters at 5 Hz
    % Fill this up
    %% Cemera parameters to do undistorted images
    this_camera = preds.cameras.(['Camera', num2str(1)]);
    K = this_camera.IntrinsicMatrix;
    R =  this_camera.rotationMatrix;
    t = this_camera.translationVector;
    M = [R; t] * K;
    K1 = K';

    % video
    vid1 = (m);
    % Init writer
    writerObj = VideoWriter(fullfile(output), 'MPEG-4');
    writerObj.Quality = 100;
    writerObj.FrameRate = FR_to_downsample;
    open(writerObj);

    camparams = cameraParameters(K = K1, RadialDistortion = this_camera.RadialDistortion,...
        TangentialDistortion= this_camera.TangentialDistortion);


    Fig_vids = figure('Position',[20 20 1800 850], 'Visible','off', 'color', 'k'); % TODO: check position and size of the figure

    for iframe = 1:length(f)
        this_frame_idx = f(iframe);
        this_frame = read(vid1, this_frame_idx);

        h1 = imagesc(undistortImage(this_frame, camparams));
        writeVideo(writerObj,getframe(Fig_vids))

    end
    close(writerObj);
    
end


function processed_frames = processFOVFrames(raw_frames)

    % Background Subtraction
    median_frame = median(raw_frames, 3); % Compute the median across frames (3rd dimension)
    subtracted_frames = bsxfun(@minus, raw_frames, median_frame);

    % Temporal Median Filtering
    temporal_window_size = 5; % Adjust this based on your data
    padded_frames = padarray(subtracted_frames, [0 0 floor(temporal_window_size/2)], 'both', 'replicate');
    filtered_frames = zeros(size(subtracted_frames));

    for i = 1:size(raw_frames, 3)
        filtered_frames(:,:,i) = median(padded_frames(:,:,i:i+temporal_window_size-1), 3);
    end

    % Spatial Gaussian Filtering
    spatial_sigma = 1; % Adjust this based on your data
    h = fspecial('gaussian', [5 5], spatial_sigma);
    smoothed_frames = imfilter(filtered_frames, h);

    % Normalization or Contrast Stretching
    processed_frames = zeros(size(smoothed_frames));
    for i = 1:size(smoothed_frames, 3)
        frame = smoothed_frames(:,:,i);
        processed_frames(:,:,i) = (frame - min(frame(:))) / (max(frame(:)) - min(frame(:)));
    end

    % % Zscore it
    % processed_frames_z = zeros(size(smoothed_frames));
    % mean_frames = mean(processed_frames, 3);
    % std_frames = std(processed_frames,[], 3);
    % for i = 1:size(processed_frames, 3)
    %     frame = processed_frames(:,:,i);
    %     processed_frames_z(:,:,i) = (frame - mean(mean_frames(:))) / mean(std_frames(:));
    % end


end

function  run_cluster_highlights(analysisstruct, frames, output, cls)
    global FR_to_downsample
    % idx = 1;
    params.nameplot=1;
    params.density_plot =0;
    params.watershed = 1;
    params.sorted = 1;
    params.markersize = 1;
    params.coarseboundary =0;
    params.do_coarse = 0;

    Fig_vids = figure('Position',[20 20 1800 850], 'Visible','off', 'color', 'k');
    writerObj = VideoWriter(fullfile(output), 'MPEG-4');
    writerObj.Quality = 100;
    writerObj.FrameRate = FR_to_downsample; % TODO: this needs to be done
    open(writerObj);


    for idx = 1:length(cls)
        % TODO: check position and size of the figure. Needs to match all figures(videos)
        % ax = axes;
       
        CLUSTER_ID = cls(idx);
        clf;
        plot_clustercolored_tsne_highlight(analysisstruct,1,params.watershed,Fig_vids,params, CLUSTER_ID); % create this function, mnodifying the original
        % idx = idx +1;
    
        % ... the rest of the code 
        writeVideo(writerObj,getframe(Fig_vids))
        
        % drawnow
    
    end
    close(writerObj);
    % Close the visualization figure

end

function run_raw_miniscope_video(raw_video, frames, output, cls)
    global FR_to_downsample
    writerObj = VideoWriter(fullfile(output), 'MPEG-4');
    writerObj.Quality = 100;
    writerObj.FrameRate = FR_to_downsample;

    % Initialize the video writer
    open(writerObj);
    
    % Create the figure for visualization
    Fig_vids = figure('Position',[20 20 1800 850], 'Visible','off', 'color', 'k');
    
    cls_idx = 1;
    % Iterate through the frames
    for frameIdx = frames
        % Extract frame data from raw miniscope video
        if frameIdx <= raw_video.NumFrames
            frameData = read(raw_video, frameIdx);
            
            % Convert to grayscale if needed
            if size(frameData, 3) == 3
                frameData = rgb2gray(frameData);
            end
            
            % Normalize the frame
            frameData = double(frameData);
            frameData = (frameData - min(frameData(:))) / (max(frameData(:)) - min(frameData(:)));
            
            % Display the frame
            imagesc(frameData);
            colormap('gray');
            clim([0, 1]);
            axis off;
            
            % Add cluster information if available
            if cls_idx <= length(cls)
                title(['Raw Miniscope - Cls: ' num2str(cls(cls_idx))], 'Color','w');
            else
                title('Raw Miniscope', 'Color','w');
            end
            
            % Save this visualization to the video
            writeVideo(writerObj, getframe(Fig_vids));
            cls_idx = cls_idx + 1;
        end
    end
    
    % Close the video writer
    close(writerObj);
    
    % Close the visualization figure
    close(Fig_vids);
end
