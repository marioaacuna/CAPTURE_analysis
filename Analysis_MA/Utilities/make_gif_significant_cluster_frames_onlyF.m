
%% Create GIF for Significant Clusters (F condition only)
% Stand alone script. 
% TODO: re-check variable names.

clear;
close all;
clc;

%% set output vars
filename = fullfile(rootpath, 'predominantF_clusters_animation_1_5width.gif');
frame_duration = 0.1; % Duration of each frame in seconds

%% INIT
% Set root path based on operating system
if ispc
    rootpath = 'D:\CAPTURE';
else
    rootpath = '/Users/mario/Library/Mobile Documents/com~apple~CloudDocs/_todo/SpontaneousPainProject/240704_data'; % I have updated where the files are to gather them more quickly
end

% File paths
analysis_filename = fullfile(rootpath, 'raw_concat_analysis.mat');
predictions_filename = fullfile(rootpath, 'agg_predictions.mat');
ratception_filename = fullfile(rootpath, 'ratception_prediction.mat');

%% Load Data
% Load analysis structure
load(analysis_filename, 'analysisstruct');

% Load predictions
load(predictions_filename, 'predictions', 'animal_condition_identifier');

% Load ratception structure
load(ratception_filename, 'ratception_struct');

%% 1. Cluster Composition Analysis
disp('Performing Cluster Composition Analysis...');

%% Prepare Data
% Upsample animal_condition_identifier
input_params.repfactor = 3;
% input_params = analysisstruct.input_params;
upsampled_identifiers = repelem(animal_condition_identifier, input_params.repfactor);

% Get the frames with good tracking
good_frames = analysisstruct.frames_with_good_tracking{1, 1};

% Extract cluster IDs and corresponding identifiers for good frames
cluster_ids = analysisstruct.annot_reordered{end};
frame_identifiers = upsampled_identifiers(good_frames);

% Extract unique animal IDs and conditions
[animal_ids, ~, animal_indices] = unique(cellfun(@(x) x(1:end-2), frame_identifiers, 'UniformOutput', false));
conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);

unique_clusters = unique(cluster_ids);
num_clusters = length(unique_clusters);
num_animals = length(animal_ids);

%%
% Initialize storage for the density of each cluster per condition
cluster_density_condition_1 = zeros(size(unique_clusters));
cluster_density_condition_0 = zeros(size(unique_clusters));

% Initialize matrices to store results
clusterComposition = zeros(num_clusters, 2);  % [S_count, F_count]

% animal_frames_ids = frame_identifiers(good_frames);
    
for c = 1:length(unique_clusters)
    cluster_index = unique_clusters(c);
    
    % Find all frames belonging to the current cluster
    frames_in_cluster = find(cluster_ids == cluster_index);
    
    % Identify the animal-condition combinations for these frames
    animal_conditions_in_cluster = frame_identifiers(frames_in_cluster);
    
    % Count frames for each condition
    cluster_density_condition_1(c) = sum(cellfun(@(x) endsWith(x, '_F'), animal_conditions_in_cluster));
    cluster_density_condition_0(c) = sum(cellfun(@(x) endsWith(x, '_S'), animal_conditions_in_cluster));
    
    
    clusterComposition(c, 1) =  sum(cellfun(@(x) endsWith(x, '_F'), animal_conditions_in_cluster));
    clusterComposition(c, 2) =  sum(cellfun(@(x) endsWith(x, '_S'), animal_conditions_in_cluster));
end

totalFrames = sum(clusterComposition, 2);
clusterProportions = clusterComposition ./ totalFrames;

% Identify predominantly associated clusters
threshold = 0.99;  % Define threshold for "predominant" association
predominantF = find(clusterProportions(:, 1) >= threshold);
predominantS = find(clusterProportions(:, 2) >= threshold);
%%
sig_clusters = predominantF;

% Set up the figure for the GIF
n_sig_clusters = length(sig_clusters);
n_rows = ceil(sqrt(n_sig_clusters));
n_cols = ceil(n_sig_clusters / n_rows);

fig = figure('Position', [100, 100, 1200, 900], Visible='on');

%% Create the GIF


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