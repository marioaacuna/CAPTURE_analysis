%% Peamble
% this script is going to output the video for a specific cluster

%% Inputs
animal_ID = 'JH_326';
camera='Camera1';
animal_number = strsplit(animal_ID, '_');
animal_number = animal_number{2};
date = []; % TODO
this_cl = [41:43]; % this can be also a group of clusters;
write_video_as_mp4 = 0;

animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
conditions = {'F', 'F', 'S', 'S', 'N', 'S', 'F', 'F', 'F'};

% Loop throug animals
%TODO

animal_condition = conditions{ismember(animal_list, animal_number)};

% switch condition to match the folder name
switch animal_condition
    case 'S'
        animal_condition = 'saline';
    case 'F'
        animal_condition = 'PFA';        
end

% Paths
dannce_root = fullfile ("H:\DANNCE\6cam_behavior", animal_condition, animal_number,'DANNCE_ready');
cluster_file = fullfile("H:\Mario\Results_Capture\clusters\clusters_struct.mat");
predictions_file = fullfile(dannce_root,'DANNCE\predict_results_net_8', 'predictions.mat' );
behavior_file = fullfile(dannce_root, 'videos', camera, '0.mp4');
analysisstruct_filename = fullfile("D:\test_CAPTURE\CAPTURE\analysisstruct_clusters.mat");


%% Load files
% read clusters vector
clusters = load(cluster_file);
clusters = clusters.clusters_struct;
cluster_this_animal = clusters.(animal_ID);
% read predictions
pred = load(predictions_file);
% pred = pred.predictions;
predictions = pred.predictions;
% read behavioral movie
beh_movie = VideoReader(behavior_file);
% Load analysis struct
load(analysisstruct_filename)

%% Start the main script
done = make_video_of_cluster(beh_movie, pred, this_cl, cluster_this_animal, animal_number, write_video_as_mp4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% Helper fucntions
function done = make_video_of_cluster(beh_movie, pred, this_cl, clusters_ds, animal_number, write_video_as_mp4)
    % Camera and Video Initialization
    this_camera = pred.cameras.Camera1; % Direct access if only one camera
    K = this_camera.IntrinsicMatrix;
    R = this_camera.rotationMatrix;
    t = this_camera.translationVector;
    M = [R; t] * K;
    K1 = K';

    vid1 = beh_movie;
    camparams = cameraParameters(K = K1, RadialDistortion = this_camera.RadialDistortion,...
        TangentialDistortion= this_camera.TangentialDistortion);

    % Cluster String Formation
    cluster_to_make = num2str(this_cl, '%d_'); % Vectorized approach
    cluster_to_make = cluster_to_make(1:end-1); % Remove last underscore

    % Video Output Preparation
    output_dir = fullfile("D:\test_CAPTURE\out_videos_example_clusters\", animal_number);
    output = fullfile(output_dir, strcat(cluster_to_make, '.mp4')); 
    if ~exist(output_dir, 'dir')
        mkdir(output_dir)
    end
    writerObj = VideoWriter(output, 'MPEG-4');
    writerObj.Quality = 100;
    writerObj.FrameRate = 120;

    % Figure Preparation
    Fig_vids = figure('Position',[20 20 1800 850], 'Visible', 'on', 'color', 'k');

    % Frame Processing
    frames_to_choose = find(ismember(clusters_ds, this_cl));
    if write_video_as_mp4
        open(writerObj);
    end
    for iframe = frames_to_choose
        % Calculate the time to seek to
        seekTime = (iframe - 1) / vid1.FrameRate;
        vid1.CurrentTime = seekTime;

        % Read the frame at the specified time
        this_frame = readFrame(vid1);
        imagesc(undistortImage(this_frame, camparams)); % Removed unnecessary handle
        title(num2str(clusters_ds(iframe)), 'Color','w');
        
        if write_video_as_mp4
            writeVideo(writerObj, getframe(Fig_vids));
        else
            drawnow;
        end
    end

    % Closing Operations
    if write_video_as_mp4
        close(writerObj);
    end
    done = 1;
end

% % function done = make_video_of_cluster(beh_movie, pred, this_cl, clusters_ds, animal_number, write_video_as_mp4);
% %     this_camera = pred.cameras.(['Camera', num2str(1)]);
% %     K = this_camera.IntrinsicMatrix;
% %     R = this_camera.rotationMatrix;
% %     t = this_camera.translationVector;
% %     M = [R; t] * K;
% %     K1 = K';
% % 
% %     vid1 = beh_movie;
% %     camparams = cameraParameters(K = K1, RadialDistortion = this_camera.RadialDistortion,...
% %         TangentialDistortion= this_camera.TangentialDistortion);
% % 
% % 
% % 
% %     % make a string out of the cluster numbers
% %     if length(this_cl) > 1
% %         cluster_to_make = strjoin(string(this_cl), '_');
% %     else
% %         cluster_to_make = num2str(this_cl);
% %     end
% %     visibility = 'on';
% %     if write_video_as_mp4
% %         % write video
% %         output_dir = fullfile("D:\test_CAPTURE\out_videos_example_clusters\", animal_number);
% %         output = fullfile("D:\test_CAPTURE\out_videos_example_clusters\", animal_number, num2str(cluster_to_make),'.mp4'); 
% %         if ~exist("output_dir", 'dir')
% %             mkdir(output_dir)
% %         end
% %         writerObj = VideoWriter(fullfile(output), 'MPEG-4');
% %         writerObj.Quality = 100;
% %         writerObj.FrameRate = 120; 
% %         open(writerObj);
% %         visibility = 'off';
% %     end
% % 
% %     Fig_vids = figure('Position',[20 20 1800 850], 'Visible', visibility, 'color', 'k'); % TODO: check position and size of the figure
% %     % run video of clusters
% %     frames_to_choose = find(ismember(clusters_ds, this_cl));
% %     f = frames_to_choose;
% %     for iframe = 1:length(f)
% %         this_frame_idx = f(iframe) ;
% %         this_frame = read(vid1, this_frame_idx);
% % 
% %         h1 = imagesc(undistortImage(this_frame, camparams));
% %         title(num2str(clusters_ds(this_frame_idx)), 'Color','w')
% % 
% %         if write_video_as_mp4
% %             writeVideo(writerObj, getframe(Fig_vids));
% %         else
% %             drawnow
% %         end
% %     end
% %     % close video
% %     if write_video_as_mp4
% %         close(writerObj);
% %     end
% %     done = 1;
% % end
% % 
% % 
% % 
% % 
