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

%% 2.1 calculate walking bouts based on 2D positions (x,y)
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
params.cycle_window = [-0.5 0.5]; % Time window for cycle normalization (seconds)
params.n_normalized_points = 200; % Number of points after interpolation
params.interp_method = 'pchip'; % Interpolation method (pchip/linear/spline)
params.markers_to_study = {'WristL', 'WristR', 'AnkleL', 'AnkleR', 'HindpawL', 'HindpawR'}; % Markers to study

free_SpineM = markers_not_aligned_ds.SpineM;
free_Snout = markers_not_aligned_ds.Snout;
[walking_bouts, metrics] = detectWalkingFrom2D(free_SpineM,free_Snout, params);

%% 3. Analysis of angles at walking in egocentric reference
markers = markers_aligned_ds;

% Calculate angles and perform analysis
analyze_angles(markers, walking_bouts, conditions, unique_conditions, frame_identifiers, animal_condition_identifier);

disp('done')

%% 4. Kinematic analysis of gait cycles
cycles_all = analyze_gait_cycles_3d(markers, walking_bouts, params);

%% 4.1 Kinematics of gait cycles per condition
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

%% 5. Plotting cycles per marker for each condition
marker_names = params.markers_to_study;
n_markers = length(marker_names);
n_conditions = length(unique_conditions);

% Determine subplot grid size
n_cols = ceil(sqrt(n_markers));
n_rows = ceil(n_markers / n_cols);

% Define colors for each condition
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
          [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
legendHandles = [];  % Store handles for legend

figure;
for m = 1:n_markers
    subplot(n_rows, n_cols, m);
    hold on;
    for c = 1:n_conditions
        condition = unique_conditions{c};
        cycles = CYCLES.(condition);
        
        % Extract data for the current marker
        marker_data = arrayfun(@(x) x.(marker_names{m}), cycles, 'UniformOutput', false);
        marker_data = cat(1, marker_data{:});
        
        % Calculate mean and SEM
        marker_mean = mean(marker_data, 1);
        marker_sem = std(marker_data, 0, 1) / sqrt(size(marker_data, 1));
        
        % Plot mean and SEM as shaded area
        t = linspace(params.cycle_window(1), params.cycle_window(2), params.n_normalized_points);
        fill([t, fliplr(t)], [marker_mean + marker_sem, fliplr(marker_mean - marker_sem)], ...
            colors{c}, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        h_line = plot(t, marker_mean, 'Color', colors{c}, 'LineWidth', 1.5);
        
        if m == 1  % Only store handles from first subplot for legend
            legendHandles = [legendHandles h_line];
        end
    end
    title(marker_names{m});
    xlabel('Time (s)');
    ylabel('Position (AU)');
    grid on;
end
legend(legendHandles, unique_conditions, 'Location', 'best');
sgtitle('Mean Gait Cycle per Marker and Condition');
hold off;

%% --- END OF SCRIPT ---

