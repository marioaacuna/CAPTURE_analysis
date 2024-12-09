% Initialization
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;
logger('Starting leg movement plot script', 'INFO');

% Load Data
logger('Loading data', 'INFO');
%load(GC.filename_analysis, 'analysisstruct')
load(GC.filename_ratception, 'ratception_struct');
load(GC.filename_predictions, 'animal_condition_identifier');
input_params.repfactor = GC.repfactor;

% Preprocess data
logger('Preprocessing data', 'INFO');
markers_aligned_ds = load_aligned_markers(ratception_struct.markers_aligned_preproc, input_params.repfactor, 15);
markers_not_aligned_ds = load_aligned_markers(ratception_struct.markers_preproc, input_params.repfactor, 15);

% Extract conditions
frame_identifiers = animal_condition_identifier;
conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);
unique_conditions = unique(conditions);

% Load parameters
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



% Output a small video with a walking example
logger('Creating walking example video', 'INFO');
s = struct();
s.markers_aligned_preproc = markers_aligned_ds;
s.markernames = ratception_struct.markernames;
s.markercolor = ratception_struct.markercolor; 
s.markers_preproc = markers_not_aligned_ds;
s.links = ratception_struct.links;
free_SpineM = markers_not_aligned_ds.SpineM;
free_Snout = markers_not_aligned_ds.Snout;
[walking_bouts, metrics] = detectWalkingFrom2D(free_SpineM,free_Snout, params);

bout_starts = find(diff([0; walking_bouts]) == 1);
bout_ends = find(diff([walking_bouts; 0]) == -1);

frames_to_plot = [];
min_frames = 100;
num_bouts = 10;

for i = 1:length(bout_starts)
    if (bout_ends(i) - bout_starts(i) + 1) >= min_frames
        frames_to_plot = [frames_to_plot, bout_starts(i):bout_ends(i)];
        if length(frames_to_plot) >= num_bouts * min_frames
            frames_to_plot = frames_to_plot(1:num_bouts * min_frames);
            break;
        end
    end
end

% Define markers to use
markers_to_use = {'SpineF', 'SpineM', 'KneeR', 'AnkleR', 'HindpawR'};

% Create GIF
logger('Creating GIF', 'INFO');
gif_filename = 'leg_movement.gif';
h = figure;

for frame = 1:length(frames_to_plot)
    % Plot 
    animate_leg_markers_aligned(s, frames_to_plot(frame), h, '', markers_to_use);
    
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

logger('Leg movement plot script completed', 'INFO');
