% 1. Initialization
clear;
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Load Data
load(GC.filename_analysis, 'analysisstruct')
load(GC.filename_ratception, 'ratception_struct');
load(GC.filename_predictions, 'animal_condition_identifier');
input_params.repfactor = GC.repfactor;

% Preprocess data
markers_aligned_ds = load_aligned_markers(ratception_struct.markers_aligned_preproc, input_params.repfactor, 50);
markers_not_aligned_ds = load_aligned_markers(ratception_struct.markers_preproc, input_params.repfactor, 50);
% Extract conditions
frame_identifiers = animal_condition_identifier;
conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);
unique_conditions = unique(conditions);

%% 2.1 calculate walking based on 2D positions
% Basic usage
% walking_bouts = detectWalkingFrom2D(markers_not_aligned_ds.SpineM);

% With custom parameters
params = struct();
params.velocity_percentile = 87; % More stringent threshold
params.min_bout_duration = 0.1; % Longer minimum bout
params.sampling_rate = 100; % Hz
params.smoothing_window = 5; % frames
params.z_smoothing_window = 5;
params.z_threshold_percentile = 99.3;  % Threshold for Z displacement
params.direction_threshold = 75;  % Max angle deviation from heading (degrees)

free_SpineM = markers_not_aligned_ds.SpineM;
free_Snout = markers_not_aligned_ds.Snout;
[walking_bouts, metrics] = detectWalkingFrom2D(free_SpineM,free_Snout, params);

%% 3. Analysis of angles at walking in egocentric reference
markers = markers_aligned_ds;

% Calculate angles and perform analysis
analyze_angles(markers, walking_bouts, conditions, unique_conditions, frame_identifiers, animal_condition_identifier);

disp('done')

%% 4. Kinematic analysis
plot_gait_trajectories(markers, walking_bouts, params)
analyze_gait_cycles_3d(markers,walking_bouts, params)

% TODO

%% --- END OF SCRIPT ---

%% Next step
% TODO: analyse trajectories of the main markers for left and right paw while walking

