% 1. Initialization
logger('Starting walking analysis script', 'INFO');
% clear;
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Load Data
logger('Loading data', 'INFO');
load(GC.filename_analysis, 'analysisstruct')
load(GC.filename_ratception, 'ratception_struct');
load(GC.filename_predictions, 'animal_condition_identifier');
input_params.repfactor = GC.repfactor;

% Preprocess data
logger('Preprocessing data', 'INFO');
markers_aligned_ds = load_aligned_markers(ratception_struct.markers_aligned_preproc, input_params.repfactor, 10);
markers_not_aligned_ds = load_aligned_markers(ratception_struct.markers_preproc, input_params.repfactor, 15);
% Extract conditions
frame_identifiers = animal_condition_identifier;
conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);
unique_conditions = unique(conditions);

%% 2. Detect pain phenotypes
% TODO: add short descripton of what they are
% With custom parameters
params = struct();
params.velocity_percentile = 87; % More stringent threshold
params.min_bout_duration = 0.1; % Longer minimum bout
params.sampling_rate = 100; % Hz
params.smoothing_window = 5; % frames
params.z_smoothing_window = 5;
params.z_threshold_percentile = 99.3;  % Threshold for Z displacement
params.direction_threshold = 75;  % Max angle deviation from heading (degrees)
params.cycle_window = [-0.5 0.5]; % Time window for cycle normalization (seconds)
params.n_normalized_points = 200; % Number of points after interpolation
params.interp_method = 'pchip'; % Interpolation method (pchip/linear/spline)
params.markers_to_study = {'WristL', 'WristR', 'KneeL', 'KneeR','AnkleL', 'AnkleR', 'HindpawL', 'HindpawR'}; % Markers to study


[pain_frames, metrics]  = detectPainPhenotypes(markers_not_aligned_ds,markers_aligned_ds, params);

% Output a small video with a walking example
logger('Creating walking example video', 'INFO');
s = struct();
s.markers_aligned_preproc = markers_aligned_ds;
s.markernames = ratception_struct.markernames;
s.markercolor =ratception_struct.markercolor; 
s.markers_preproc = markers_not_aligned_ds;
s.links = ratception_struct.links;
h = figure;
frames_to_plot = find(pain_frames);
% g = animate_markers_aligned_fullmovie_demo(s,frames_to_plot, h , '');

% Create GIF
logger('Creating GIF', 'INFO');
gif_filename = 'pain_frames.gif';
h = figure;
for frame = 1:length(frames_to_plot)
    % Plot 
    animate_markers_aligned_fullmovie_demo(s,frames_to_plot(frame), h , '');

    
    % Capture the plot as an image
    frame_img = getframe(h);
    img = frame2im(frame_img);
    [imind, cm] = rgb2ind(img, 256);
    
    % Write to the GIF File
    if frame == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end