% 1. Initialization
clc
logger('Starting walking analysis script', 'INFO');
% clear;
close all;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Load Data
logger('Loading data', 'INFO');

data_path = '/media/mario/no backup/CAPTURE/data/0_preprocessing_BSFC_300hz/';

load(fullfile(data_path, "raw_concat_analysis.mat"), 'analysisstruct')
load(fullfile(data_path, "ratception_prediction.mat"), 'ratception_struct');
load(fullfile(data_path,"agg_predictions.mat"), 'animal_condition_identifier');
input_params.repfactor = GC.repfactor;

% Preprocess data
logger('Preprocessing data', 'INFO');
markers_aligned_ds = load_aligned_markers(ratception_struct.markers_aligned_preproc, input_params.repfactor, 10);
markers_not_aligned_ds = load_aligned_markers(ratception_struct.markers_preproc, input_params.repfactor, 15);
% Extract conditions
frame_identifiers = animal_condition_identifier;
conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);
unique_conditions = unique(conditions, 'stable');
animals_condition = unique(frame_identifiers, 'stable');
animals_comditions_cell = cellfun(@(x) strsplit(x, '_'), animals_condition, 'UniformOutput', false);

%% 2. Detect pain phenotypes
% TODO: add short descripton of what they are
% With custom parameters
params = struct();
params.velocity_percentile = 87; % More stringent threshold
params.min_bout_duration = 0.1; % Longer minimum bout
params.sampling_rate = 100; % Hz
params.smoothing_window = 200; % frames
params.z_smoothing_window = 5;
params.z_threshold_percentile = 99.3;  % Threshold for Z displacement
params.direction_threshold = 75;  % Max angle deviation from heading (degrees)
params.cycle_window = [-0.5 0.5]; % Time window for cycle normalization (seconds)
params.n_normalized_points = 200; % Number of points after interpolation
params.interp_method = 'pchip'; % Interpolation method (pchip/linear/spline)
params.markers_to_study = {'WristL', 'WristR', 'KneeL', 'KneeR','AnkleL', 'AnkleR', 'HindpawL', 'HindpawR'}; % Markers to study

%% Detect licking/attention
%[pain_frames, metrics]  = detectPainPhenotypes(markers_not_aligned_ds,markers_aligned_ds, params);
params.pain_threshold = 0.59;
params.temporal_window = 5; % 0.1 seconds at 100 Hz, example.

params.velocity_percentile = 5; % For movement threshold, default: 15
params.feature_smoothing_window = 5; % For temporal smoothing
params.max_smoothing_window = 20; % Maximum smoothing window
[pain_frames, metrics]  = detectPainPhenotypes_v2(markers_not_aligned_ds,markers_aligned_ds, params);


% TEst
% pain_frames = metrics.paw_licking_detected;
% %% Output a small video with a walking example
% logger('Creating walking example video', 'INFO');
% s = struct();
% s.markers_aligned_preproc = markers_aligned_ds;
% s.markernames = ratception_struct.markernames;
% s.markercolor =ratception_struct.markercolor; 
% s.markers_preproc = markers_not_aligned_ds;
% s.links = ratception_struct.links;
% frames_to_plot = find(pain_frames);
% % g = animate_markers_aligned_fullmovie_demo(s,frames_to_plot, h , '');
% 
% % Create GIF
% logger('Creating GIF', 'INFO');
% gif_filename = 'pain_frames_BSFC.gif';
% h = figure;
% for frame = 1:length(frames_to_plot)
%     % Plot 
%     animate_markers_aligned_fullmovie_demo(s,frames_to_plot(frame), h , '');
% 
% 
%     % Capture the plot as an image
%     frame_img = getframe(h);
%     img = frame2im(frame_img);
%     [imind, cm] = rgb2ind(img, 256);
% 
%     % Write to the GIF File
%     if frame == 1
%         imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
%     else
%         imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
%     end
% end

%% Calculate proportion of pain frames per animal and condition
logger('Calculating proportion of pain frames per animal and condition', 'INFO');
animal_ids = unique(frame_identifiers);
n_animals = length(animal_ids);
pain_proportions = nan(n_animals, length(unique_conditions));

for a = 1:n_animals
    animal_id = animal_ids{a};
    animal_mask = strcmp(frame_identifiers, animal_id);
    total_frames = sum(animal_mask);

    c = (ismember(unique_conditions, animal_id(end)));

    pain_frames_animal_condition = pain_frames & animal_mask;
    pain_proportions(a, c) = sum(pain_frames_animal_condition);

end

% Convert proportions to time (seconds)
pain_times = pain_proportions / params.sampling_rate;

% %% Perform ANOVA and posthoc comparisons
% logger('Performing ANOVA and posthoc comparisons', 'INFO');
% [p, tbl, stats] = anova1(pain_times);
% % c = multcompare(stats, 'CType', 'bonferroni');
% 
% c = multcompare(stats);
% 
% 
% % Extract significant comparisons
% significant_comparisons = c(c(:,6) < 0.05, :);

%% Plot bar plot with errorbars for each condition
logger('Plotting bar plot with errorbars for each condition', 'INFO');
mean_pain_times = nanmean(pain_times, 1);
sem_pain_times = zeros(1, length(unique_conditions));

% Perform t-test comparisons
logger('Performing t-test comparisons', 'INFO');
for c = 1:length(unique_conditions)
    condition_mask = strcmp(conditions, unique_conditions{c});
    n_animals_condition = sum(condition_mask);
    data = pain_times(~isnan(pain_times(:, c)),c);
    sem_pain_times(c) = std(data) / sqrt(length(data));
    
    % Perform t-test comparisons
    if strcmp(unique_conditions{c}, 'F')
        F_data = data;
    elseif strcmp(unique_conditions{c}, 'S')
        S_data = data;
    elseif strcmp(unique_conditions{c}, 'N')
        N_data = data;
    elseif strcmp(unique_conditions{c}, 'H')
        H_data = data;
    end
end

% T-test for F vs S
if exist('F_data', 'var') && exist('S_data', 'var')
    [~, p_F_vs_S] = ttest(F_data, S_data);
    disp(['p-value for F vs S: ', num2str(p_F_vs_S)]);
end

% T-test2 for N vs H
if exist('N_data', 'var') && exist('H_data', 'var')
    [~, p_N_vs_H] = ttest2(N_data, H_data);
    disp(['p-value for N vs H: ', num2str(p_N_vs_H)]);
end

figure;
bar(mean_pain_times, 'FaceColor', 'flat');
hold on;
errorbar(1:length(unique_conditions), mean_pain_times, sem_pain_times, 'k.', 'LineWidth', 1.5);
set(gca, 'XTickLabel', unique_conditions);
xlabel('Condition');
ylabel('Time of Pain Frames (s)');
title('Time of Pain Frames per Condition');
grid on;

% Add stars and horizontal bars for significant comparisons
hold on;
y_max = max(mean_pain_times + sem_pain_times) * 1.1;
for i = 1:size(significant_comparisons, 1)
    group1 = significant_comparisons(i, 1);
    group2 = significant_comparisons(i, 2);
    p_value = significant_comparisons(i, 6);
    
    % Determine star symbol based on p-value
    if p_value < 0.001
        star = '***';
    elseif p_value < 0.01
        star = '**';
    elseif p_value < 0.05
        star = '*';
    else
        star = '';
    end
    
    % Plot horizontal bar and star
    plot([group1, group2], [y_max, y_max], 'k-', 'LineWidth', 1.5);
    text(mean([group1, group2]), y_max, star, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
    y_max = y_max * 1.05; % Increment y_max for next comparison
end
hold off;