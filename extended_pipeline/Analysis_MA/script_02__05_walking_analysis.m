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
markers_aligned_ds = load_aligned_markers(ratception_struct.markers_aligned_preproc, input_params.repfactor, 10);
markers_not_aligned_ds = load_aligned_markers(ratception_struct.markers_preproc, input_params.repfactor, 15);
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
cycles = analyze_gait_cycles_3d(markers,walking_bouts, params);

% TODO: Separate the analysis into walking bouts for each condition

CYCLES = struct();
% Separate data by conditions
for c = 1:length(unique_conditions)
    condition = unique_conditions{c};
    condition_mask = strcmp(conditions, condition) ;
    
    % Perform gait analysis for each condition
    fprintf('Gait analysis for condition: %s\n', condition);
    
    % Restrict markers to current condition
    condition_markers = structfun(@(x) x(condition_mask, :), markers, 'UniformOutput', false);
    
    % Analyze gait cycles for the current condition
    cycles = analyze_gait_cycles_3d(condition_markers, walking_bouts(strcmp(conditions, condition)), params);
    
    CYCLES.(condition) = cycles;
end

disp('done')

%% --- END OF SCRIPT ---

