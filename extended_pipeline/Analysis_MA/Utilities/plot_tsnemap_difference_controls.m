%% Initialization
global logger GC
%logger('test tsne map difference controls vs baseline', 'INFO'); % Modify
clear;
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Configuration for visualization and export
debugging = true;  % set to true for debugging mode
if debugging
    visualize = 'on';  % Show figures during debugging
    do_export = false;  % Don't export during debugging
else
    visualize = 'off';  % Don't show figures in production mode
    do_export = true;   % Export figures in production mode
end

% Export folder
export_folder = fullfile(GC.temp_root, 'figs_presentation_painAI_controls');
if ~exist(export_folder, 'dir')
    mkdir(export_folder);
end

%% Load Data
logger('Loading data', 'INFO');
% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');

% Load predictions
load(GC.filename_predictions, 'predictions', 'animal_condition_identifier');

% Load ratception structure
load(GC.filename_ratception, 'ratception_struct');

% Extract animal list and conditions
upsamplig_factor = GC.repfactor;
long_animal_frames_identifier = repelem(animal_condition_identifier,upsamplig_factor);
animal_list_used_after_analysis =  long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});

%%
zvals = analysisstruct.zValues;
cls_all = analysisstruct.annot_reordered{2}; % list of all clusters on the t-SNE map

%% Calculate difference between baseline and control conditions
% Control conditions to compare with baseline ('B')
conditions_to_compare = {{'B'}, {'S'}; {'B'}, {'H'}}; % B vs S and B vs H
% Initialize a structure to hold the differences
diff_struct = struct();

%% Extract cluster information for each comparison with animal ID matching
logger('Analyzing cluster differences between baseline and control conditions', 'INFO');

% Get cluster assignments for all frames
cluster_assignments = analysisstruct.annot_reordered{end,end}; % Final cluster assignments
unique_clusters = unique(cluster_assignments);
unique_clusters = unique_clusters(unique_clusters > 0); % Remove background/noise clusters

% Map cluster assignments to animal conditions and extract animal IDs
frame_conditions = cell(length(cluster_assignments), 1);
frame_animal_ids = cell(length(cluster_assignments), 1);
for i = 1:length(cluster_assignments)
    if i <= length(animal_list_used_after_analysis)
        this_a = animal_list_used_after_analysis{i};
        this_cond = this_a(end); % Get last character as condition
        % Extract animal ID (everything before the last '_')
        underscore_pos = find(this_a == '_', 1, 'last');
        if ~isempty(underscore_pos)
            this_animal_id = this_a(1:underscore_pos-1);
        else
            this_animal_id = this_a; % fallback if no underscore found
        end
        frame_conditions{i} = this_cond;
        frame_animal_ids{i} = this_animal_id;
    else
        frame_conditions{i} = 'Unknown';
        frame_animal_ids{i} = 'Unknown';
    end
end

%% Process each condition comparison (baseline vs controls) with animal matching
for comp_idx = 1:size(conditions_to_compare, 1)
    baseline_conditions = conditions_to_compare{comp_idx, 1}; % Always 'B'
    control_conditions = conditions_to_compare{comp_idx, 2}; % 'S' or 'H'
    
    comparison_name = sprintf('%s_vs_%s', strjoin(control_conditions, ''), strjoin(baseline_conditions, ''));
    logger(sprintf('Processing comparison: %s (with animal ID matching)', comparison_name), 'INFO');
    
    % Find animals that have both baseline and control conditions
    baseline_animal_ids = {};
    control_animal_ids = {};
    
    % Collect animal IDs for baseline frames
    for i = 1:length(frame_conditions)
        if any(strcmp(frame_conditions{i}, baseline_conditions))
            baseline_animal_ids{end+1} = frame_animal_ids{i};
        end
    end
    
    % Collect animal IDs for control frames
    for i = 1:length(frame_conditions)
        if any(strcmp(frame_conditions{i}, control_conditions))
            control_animal_ids{end+1} = frame_animal_ids{i};
        end
    end
    
    % Find common animal IDs (animals that have both baseline and control data)
    common_animal_ids = intersect(unique(baseline_animal_ids), unique(control_animal_ids));
    logger(sprintf('Found %d animals with both baseline and %s conditions', ...
        length(common_animal_ids), strjoin(control_conditions, ',')), 'INFO');
    
    if isempty(common_animal_ids)
        logger(sprintf('Warning: No common animals found for comparison %s', comparison_name), 'WARNING');
        continue;
    end
    
    % Find frames belonging to each condition for matched animals only
    baseline_frames = false(length(frame_conditions), 1);
    control_frames = false(length(frame_conditions), 1);
    
    for i = 1:length(frame_conditions)
        animal_id = frame_animal_ids{i};
        if any(strcmp(animal_id, common_animal_ids))
            if any(strcmp(frame_conditions{i}, baseline_conditions))
                baseline_frames(i) = true;
            elseif any(strcmp(frame_conditions{i}, control_conditions))
                control_frames(i) = true;
            end
        end
    end
    
    % Calculate cluster prevalence for each condition
    baseline_cluster_counts = zeros(length(unique_clusters), 1);
    control_cluster_counts = zeros(length(unique_clusters), 1);
    total_cluster_frames = zeros(length(unique_clusters), 1);
    
    for c_idx = 1:length(unique_clusters)
        cluster_id = unique_clusters(c_idx);
        cluster_frames = cluster_assignments == cluster_id;
        
        baseline_cluster_counts(c_idx) = sum(cluster_frames' & baseline_frames);
        control_cluster_counts(c_idx) = sum(cluster_frames' & control_frames);
        total_cluster_frames(c_idx) = sum(cluster_frames); % Total frames in this cluster
    end
    
    % Calculate total frames for normalization (for reference)
    total_baseline_frames = sum(baseline_frames);
    total_control_frames = sum(control_frames);
    
    % Calculate cluster-normalized prevalence (percentage within each cluster)
    baseline_prevalence = zeros(length(unique_clusters), 1);
    control_prevalence = zeros(length(unique_clusters), 1);
    
    for c_idx = 1:length(unique_clusters)
        if total_cluster_frames(c_idx) > 0
            baseline_prevalence(c_idx) = baseline_cluster_counts(c_idx) / total_cluster_frames(c_idx);
            control_prevalence(c_idx) = control_cluster_counts(c_idx) / total_cluster_frames(c_idx);
        end
    end
    
    % Calculate difference (control - baseline) normalized by cluster size
    cluster_difference = control_prevalence - baseline_prevalence;
    
    % Identify clusters that are gained (higher in control) or lost (higher in baseline)
    gained_clusters = unique_clusters(cluster_difference > 0);
    lost_clusters = unique_clusters(cluster_difference < 0);
    
    % Store results
    diff_struct.(comparison_name).baseline_conditions = baseline_conditions;
    diff_struct.(comparison_name).control_conditions = control_conditions;
    diff_struct.(comparison_name).common_animal_ids = common_animal_ids;
    diff_struct.(comparison_name).cluster_ids = unique_clusters;
    diff_struct.(comparison_name).baseline_prevalence = baseline_prevalence;
    diff_struct.(comparison_name).control_prevalence = control_prevalence;
    diff_struct.(comparison_name).cluster_difference = cluster_difference;
    diff_struct.(comparison_name).gained_clusters = gained_clusters;
    diff_struct.(comparison_name).lost_clusters = lost_clusters;
    diff_struct.(comparison_name).total_baseline_frames = total_baseline_frames;
    diff_struct.(comparison_name).total_control_frames = total_control_frames;
    
    logger(sprintf('Found %d gained clusters and %d lost clusters for %s (matched animals: %d)', ...
        length(gained_clusters), length(lost_clusters), comparison_name, length(common_animal_ids)), 'INFO');
end

%% Create visualization function for t-SNE difference maps (baseline vs controls)
function plot_tsne_difference_map_controls(analysisstruct, diff_data, comparison_name, visualize)
    % Create figure
    fig = figure('Visible', visualize);
    set(fig, 'Position', [100, 100, 1200, 800]);
    set(fig, 'Color', 'w');
    
    % Plot base t-SNE map (all points in light gray)
    subplot(1, 3, 1);
    plot(analysisstruct.zValues(:,1), analysisstruct.zValues(:,2), '.', ...
        'Color', [0.8, 0.8, 0.8], 'MarkerSize', 2);
    hold on;
    
    % Highlight gained clusters in red
    for i = 1:length(diff_data.gained_clusters)
        cluster_id = diff_data.gained_clusters(i);
        cluster_frames = analysisstruct.annot_reordered{end,end} == cluster_id;
        if sum(cluster_frames) > 0
            plot(analysisstruct.zValues(cluster_frames,1), ...
                 analysisstruct.zValues(cluster_frames,2), '.', ...
                 'Color', 'red', 'MarkerSize', 4);
        end
    end
    
    title(sprintf('Gained Clusters (%s > %s)', ...
        strjoin(diff_data.control_conditions, ','), ...
        strjoin(diff_data.baseline_conditions, ',')));
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    axis equal;
    grid on;
    
    % Plot lost clusters
    subplot(1, 3, 2);
    plot(analysisstruct.zValues(:,1), analysisstruct.zValues(:,2), '.', ...
        'Color', [0.8, 0.8, 0.8], 'MarkerSize', 2);
    hold on;
    
    % Highlight lost clusters in blue
    for i = 1:length(diff_data.lost_clusters)
        cluster_id = diff_data.lost_clusters(i);
        cluster_frames = analysisstruct.annot_reordered{end,end} == cluster_id;
        if sum(cluster_frames) > 0
            plot(analysisstruct.zValues(cluster_frames,1), ...
                 analysisstruct.zValues(cluster_frames,2), '.', ...
                 'Color', 'blue', 'MarkerSize', 4);
        end
    end
    
    title(sprintf('Lost Clusters (%s < %s)', ...
        strjoin(diff_data.control_conditions, ','), ...
        strjoin(diff_data.baseline_conditions, ',')));
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    axis equal;
    grid on;
    
    % Plot combined difference map
    subplot(1, 3, 3);
    plot(analysisstruct.zValues(:,1), analysisstruct.zValues(:,2), '.', ...
        'Color', [0.8, 0.8, 0.8], 'MarkerSize', 2);
    hold on;
    
    % Gained clusters in red
    for i = 1:length(diff_data.gained_clusters)
        cluster_id = diff_data.gained_clusters(i);
        cluster_frames = analysisstruct.annot_reordered{end,end} == cluster_id;
        if sum(cluster_frames) > 0
            plot(analysisstruct.zValues(cluster_frames,1), ...
                 analysisstruct.zValues(cluster_frames,2), '.', ...
                 'Color', 'red', 'MarkerSize', 4);
        end
    end
    
    % Lost clusters in blue
    for i = 1:length(diff_data.lost_clusters)
        cluster_id = diff_data.lost_clusters(i);
        cluster_frames = analysisstruct.annot_reordered{end,end} == cluster_id;
        if sum(cluster_frames) > 0
            plot(analysisstruct.zValues(cluster_frames,1), ...
                 analysisstruct.zValues(cluster_frames,2), '.', ...
                 'Color', 'blue', 'MarkerSize', 4);
        end
    end
    
    % Add legend
    legend({'Background', 'Gained (Control > Baseline)', 'Lost (Control < Baseline)'}, ...
        'Location', 'best');
    
    title(sprintf('Combined Difference Map: %s', comparison_name));
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    axis equal;
    grid on;
    
    % Add main title with animal count
    sgtitle(sprintf('t-SNE Cluster Differences: %s (n=%d animals)', ...
        strrep(comparison_name, '_', ' vs '), length(diff_data.common_animal_ids)), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

%% Create enhanced watershed-based visualization for controls
function plot_tsne_watershed_difference_controls(analysisstruct, diff_data, comparison_name, visualize)
    % Create figure with watershed boundaries
    fig = figure('Visible', visualize);
    set(fig, 'Position', [100, 100, 1400, 600]);
    set(fig, 'Color', 'w');
    
    % Define colors for different cluster states
    gained_color = [1, 0, 0]; % Red for gained
    lost_color = [0, 0, 1];   % Blue for lost
    neutral_color = [0.7, 0.7, 0.7]; % Gray for unchanged
    
    % Create masks for different cluster types
    gained_mask = false(size(analysisstruct.zValues, 1), 1);
    lost_mask = false(size(analysisstruct.zValues, 1), 1);
    
    for i = 1:length(diff_data.gained_clusters)
        cluster_id = diff_data.gained_clusters(i);
        gained_mask = gained_mask | (analysisstruct.annot_reordered{end,end} == cluster_id)';
    end
    
    for i = 1:length(diff_data.lost_clusters)
        cluster_id = diff_data.lost_clusters(i);
        lost_mask = lost_mask | (analysisstruct.annot_reordered{end,end} == cluster_id)';
    end
    
    % Plot with watershed boundaries if available
    subplot(1, 2, 1);
    
    % Plot all points first
    plot(analysisstruct.zValues(:,1), analysisstruct.zValues(:,2), '.', ...
        'Color', neutral_color, 'MarkerSize', 1);
    hold on;
    
    % Overlay gained clusters
    if sum(gained_mask) > 0
        plot(analysisstruct.zValues(gained_mask,1), analysisstruct.zValues(gained_mask,2), '.', ...
            'Color', gained_color, 'MarkerSize', 3);
    end
    
    % Overlay lost clusters
    if sum(lost_mask) > 0
        plot(analysisstruct.zValues(lost_mask,1), analysisstruct.zValues(lost_mask,2), '.', ...
            'Color', lost_color, 'MarkerSize', 3);
    end
    
    % Add watershed boundaries if available
    if isfield(analysisstruct, 'sorted_watershed') && isfield(analysisstruct, 'xx') && isfield(analysisstruct, 'yy')
        nnn = analysisstruct.sorted_watershed;
        nnn(nnn > 0) = 1;
        B = bwboundaries(nnn);
        
        for kk = 1:numel(B)
            if size(B{kk}, 1) > 0
                plot(analysisstruct.xx(B{kk}(:,2)), analysisstruct.yy(B{kk}(:,1)), ...
                    'k-', 'LineWidth', 0.5);
            end
        end
    end
    
    title('Control vs Baseline Difference Map');
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    axis equal;
    legend({'Unchanged', 'Gained', 'Lost', 'Cluster Boundaries'}, 'Location', 'best');
    
    % Create difference intensity plot
    subplot(1, 2, 2);
    
    % Create a colormap based on difference magnitude
    difference_values = zeros(size(analysisstruct.zValues, 1), 1);
    
    for i = 1:length(diff_data.cluster_ids)
        cluster_id = diff_data.cluster_ids(i);
        cluster_mask = analysisstruct.annot_reordered{end,end} == cluster_id;
        difference_values(cluster_mask) = diff_data.cluster_difference(i);
    end
    
    % Create scatter plot colored by difference
    scatter(analysisstruct.zValues(:,1), analysisstruct.zValues(:,2), 10, difference_values, 'filled');
    
    % Use a blue-white-red colormap if available, otherwise use default
    try
        colormap(bluewhitered(256)); % Blue-white-red colormap
    catch
        colormap(jet(256)); % Alternative colormap
    end
    
    colorbar;
    caxis([-max(abs(difference_values)), max(abs(difference_values))]);
    
    title('Cluster Difference Intensity (Control vs Baseline)');
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    axis equal;
    
    % Add main title
    sgtitle(sprintf('Enhanced Control Analysis: %s (n=%d animals)', ...
        strrep(comparison_name, '_', ' vs '), length(diff_data.common_animal_ids)), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

%% Generate visualizations for each comparison
for comp_idx = 1:size(conditions_to_compare, 1)
    baseline_conditions = conditions_to_compare{comp_idx, 1};
    control_conditions = conditions_to_compare{comp_idx, 2};
    comparison_name = sprintf('%s_vs_%s', strjoin(control_conditions, ''), strjoin(baseline_conditions, ''));
    
    if isfield(diff_struct, comparison_name)
        plot_tsne_difference_map_controls(analysisstruct, diff_struct.(comparison_name), ...
            comparison_name, visualize);
        
        plot_tsne_watershed_difference_controls(analysisstruct, diff_struct.(comparison_name), ...
            comparison_name, visualize);
    end
end

%% Create summary statistics table for control comparisons
logger('Creating summary statistics for control vs baseline comparisons', 'INFO');

% Create a summary table of differences
summary_table = table();
comparison_names = {};
num_animals = [];
num_gained = [];
num_lost = [];
total_clusters = [];
max_gain = [];
max_loss = [];
top_20_increased_clusters = {};
bottom_20_decreased_clusters = {};
num_top_20_increased = [];
num_bottom_20_decreased = [];
baseline_frames_count = [];
control_frames_count = [];

comp_idx = 1;
field_names = fieldnames(diff_struct);
for i = 1:length(field_names)
    comparison_name = field_names{i};
    data = diff_struct.(comparison_name);
    
    comparison_names{comp_idx} = strrep(comparison_name, '_', ' vs ');
    num_animals(comp_idx) = length(data.common_animal_ids);
    num_gained(comp_idx) = length(data.gained_clusters);
    num_lost(comp_idx) = length(data.lost_clusters);
    total_clusters(comp_idx) = length(data.cluster_ids);
    max_gain(comp_idx) = max(data.cluster_difference);
    max_loss(comp_idx) = min(data.cluster_difference);
    baseline_frames_count(comp_idx) = data.total_baseline_frames;
    control_frames_count(comp_idx) = data.total_control_frames;
    
    % Calculate z-scores for cluster differences
    cluster_diff_zscore = zscore(data.cluster_difference);
    
    % Find top 20% increased clusters (highest z-scores)
    n_clusters = length(data.cluster_ids);
    top_20_threshold = ceil(0.2 * n_clusters);
    bottom_20_threshold = ceil(0.2 * n_clusters);
    
    % Sort by z-score to get top and bottom clusters
    [sorted_zscores, sort_idx] = sort(cluster_diff_zscore, 'descend');
    
    % Top 20% increased (highest positive z-scores)
    top_20_idx = sort_idx(1:top_20_threshold);
    top_20_increased_cluster_ids = data.cluster_ids(top_20_idx)';
    top_20_increased_zscores = sorted_zscores(1:top_20_threshold);
    
    % Bottom 20% decreased (lowest negative z-scores)
    bottom_20_idx = sort_idx(end-bottom_20_threshold+1:end);
    bottom_20_decreased_cluster_ids = data.cluster_ids(bottom_20_idx);
    bottom_20_decreased_zscores = sorted_zscores(end-bottom_20_threshold+1:end);
    
    % Store results
    top_20_increased_clusters{comp_idx} = top_20_increased_cluster_ids;
    bottom_20_decreased_clusters{comp_idx} = bottom_20_decreased_cluster_ids;
    num_top_20_increased(comp_idx) = length(top_20_increased_cluster_ids);
    num_bottom_20_decreased(comp_idx) = length(bottom_20_decreased_cluster_ids);
    
    % Store z-score results in diff_struct for later use
    diff_struct.(comparison_name).cluster_diff_zscore = cluster_diff_zscore;
    diff_struct.(comparison_name).top_20_increased_clusters = top_20_increased_cluster_ids;
    diff_struct.(comparison_name).bottom_20_decreased_clusters = bottom_20_decreased_cluster_ids;
    diff_struct.(comparison_name).top_20_increased_zscores = top_20_increased_zscores;
    diff_struct.(comparison_name).bottom_20_decreased_zscores = bottom_20_decreased_zscores;
    
    comp_idx = comp_idx + 1;
end

summary_table.Comparison = comparison_names';
summary_table.Matched_Animals = num_animals';
summary_table.Baseline_Frames = baseline_frames_count';
summary_table.Control_Frames = control_frames_count';
summary_table.Total_Clusters = total_clusters';
summary_table.Gained_Clusters = num_gained';
summary_table.Lost_Clusters = num_lost';
summary_table.Max_Gain = max_gain';
summary_table.Max_Loss = max_loss';
summary_table.Top_20pct_Increased_Count = num_top_20_increased';
summary_table.Bottom_20pct_Decreased_Count = num_bottom_20_decreased';
summary_table.Top_20pct_Increased_Clusters = top_20_increased_clusters';
summary_table.Bottom_20pct_Decreased_Clusters = bottom_20_decreased_clusters';

% Display summary
disp('=== Control vs Baseline t-SNE Cluster Difference Summary ===');
disp(summary_table);

% Save summary table if export is enabled
if do_export
    summary_filename = fullfile(export_folder, 'control_vs_baseline_difference_summary.csv');
    writetable(summary_table, summary_filename);
    logger(sprintf('Saved control vs baseline summary table: %s', summary_filename), 'INFO');
end

%% Create per-animal analysis (optional detailed breakdown)
logger('Creating per-animal detailed analysis', 'INFO');

% For each comparison, analyze individual animal contributions
field_names = fieldnames(diff_struct);
per_animal_analysis = struct();

for i = 1:length(field_names)
    comparison_name = field_names{i};
    data = diff_struct.(comparison_name);
    
    per_animal_analysis.(comparison_name) = struct();
    per_animal_analysis.(comparison_name).animal_ids = data.common_animal_ids;
    
    % For each animal, calculate their individual contribution to cluster differences
    animal_contributions = zeros(length(data.common_animal_ids), length(data.cluster_ids));
    
    for animal_idx = 1:length(data.common_animal_ids)
        animal_id = data.common_animal_ids{animal_idx};
        
        % Find frames for this specific animal in baseline and control conditions
        animal_baseline_frames = false(length(frame_conditions), 1);
        animal_control_frames = false(length(frame_conditions), 1);
        
        for j = 1:length(frame_conditions)
            if strcmp(frame_animal_ids{j}, animal_id)
                if any(strcmp(frame_conditions{j}, data.baseline_conditions))
                    animal_baseline_frames(j) = true;
                elseif any(strcmp(frame_conditions{j}, data.control_conditions))
                    animal_control_frames(j) = true;
                end
            end
        end
        
        % Calculate cluster prevalence for this animal
        for c_idx = 1:length(data.cluster_ids)
            cluster_id = data.cluster_ids(c_idx);
            cluster_frames = cluster_assignments == cluster_id;
            
            animal_baseline_count = sum(cluster_frames' & animal_baseline_frames);
            animal_control_count = sum(cluster_frames' & animal_control_frames);
            
            total_animal_baseline = sum(animal_baseline_frames);
            total_animal_control = sum(animal_control_frames);
            
            if total_animal_baseline > 0 && total_animal_control > 0
                baseline_prev = animal_baseline_count / total_animal_baseline;
                control_prev = animal_control_count / total_animal_control;
                animal_contributions(animal_idx, c_idx) = control_prev - baseline_prev;
            end
        end
    end
    
    per_animal_analysis.(comparison_name).contributions = animal_contributions;
    per_animal_analysis.(comparison_name).cluster_ids = data.cluster_ids;
end

%% Analyze cluster proportions with 80% threshold
logger('Analyzing cluster proportions with 80% threshold', 'INFO');

% Define threshold for significant cluster changes
threshold = 0.8; % 80% threshold

% Initialize arrays for proportion analysis
proportion_analysis = table();
comparison_names_prop = {};
total_clusters_prop = [];
num_gained_80 = [];
num_lost_80 = [];
proportion_gained = [];
proportion_lost = [];
proportion_unchanged = [];

comp_idx = 1;
field_names = fieldnames(diff_struct);
for i = 1:length(field_names)
    comparison_name = field_names{i};
    data = diff_struct.(comparison_name);
    
    % Count clusters with > 80% increase (gained)
    gained_80_clusters = sum(data.cluster_difference > threshold);
    
    % Count clusters with < -80% decrease (lost)
    lost_80_clusters = sum(data.cluster_difference < -threshold);
    
    % Count unchanged clusters (between -80% and +80%)
    unchanged_clusters = sum(abs(data.cluster_difference) <= threshold);
    
    % Calculate proportions
    total_clusters = length(data.cluster_ids);
    prop_gained = gained_80_clusters / total_clusters;
    prop_lost = lost_80_clusters / total_clusters;
    prop_unchanged = unchanged_clusters / total_clusters;
    
    % Store results
    comparison_names_prop{comp_idx} = strrep(comparison_name, '_', ' vs ');
    total_clusters_prop(comp_idx) = total_clusters;
    num_gained_80(comp_idx) = gained_80_clusters;
    num_lost_80(comp_idx) = lost_80_clusters;
    proportion_gained(comp_idx) = prop_gained;
    proportion_lost(comp_idx) = prop_lost;
    proportion_unchanged(comp_idx) = prop_unchanged;
    
    % Store in diff_struct for later use
    diff_struct.(comparison_name).gained_80_clusters = gained_80_clusters;
    diff_struct.(comparison_name).lost_80_clusters = lost_80_clusters;
    diff_struct.(comparison_name).unchanged_clusters = unchanged_clusters;
    diff_struct.(comparison_name).proportion_gained = prop_gained;
    diff_struct.(comparison_name).proportion_lost = prop_lost;
    diff_struct.(comparison_name).proportion_unchanged = prop_unchanged;
    diff_struct.(comparison_name).threshold_used = threshold;
    
    logger(sprintf('For %s: %d gained (%.1f%%), %d lost (%.1f%%), %d unchanged (%.1f%%)', ...
        comparison_names_prop{comp_idx}, gained_80_clusters, prop_gained*100, ...
        lost_80_clusters, prop_lost*100, unchanged_clusters, prop_unchanged*100), 'INFO');
    
    comp_idx = comp_idx + 1;
end

% Create proportion analysis table
proportion_analysis.Comparison = comparison_names_prop';
proportion_analysis.Total_Clusters = total_clusters_prop';
proportion_analysis.Gained_80pct_Count = num_gained_80';
proportion_analysis.Lost_80pct_Count = num_lost_80';
proportion_analysis.Proportion_Gained = proportion_gained';
proportion_analysis.Proportion_Lost = proportion_lost';
proportion_analysis.Proportion_Unchanged = proportion_unchanged';

% Display proportion analysis
disp('=== Cluster Proportion Analysis (80% Threshold) ===');
disp(proportion_analysis);

%% Create pie charts for cluster proportions
function create_cluster_proportion_pie_charts(diff_struct, threshold, visualize, export_folder, do_export)
    field_names = fieldnames(diff_struct);
    
    for i = 1:length(field_names)
        comparison_name = field_names{i};
        data = diff_struct.(comparison_name);
        
        % Create figure for pie chart
        fig = figure('Visible', visualize);
        set(fig, 'Position', [100, 100, 800, 600]);
        set(fig, 'Color', 'w');
        
        % Prepare data for pie chart
        gained_count = data.gained_80_clusters;
        lost_count = data.lost_80_clusters;
        unchanged_count = data.unchanged_clusters;
        
        % Only include non-zero categories
        pie_data = [];
        pie_labels = {};
        pie_colors = [];
        
        if gained_count > 0
            pie_data = [pie_data, gained_count];
            pie_labels{end+1} = sprintf('Gained (>%d%%) - %d clusters', threshold*100, gained_count);
            pie_colors = [pie_colors; 1, 0, 0]; % Red
        end
        
        if lost_count > 0
            pie_data = [pie_data, lost_count];
            pie_labels{end+1} = sprintf('Lost (<%d%%) - %d clusters', -threshold*100, lost_count);
            pie_colors = [pie_colors; 0, 0, 1]; % Blue
        end
        
        if unchanged_count > 0
            pie_data = [pie_data, unchanged_count];
            pie_labels{end+1} = sprintf('Unchanged (±%d%%) - %d clusters', threshold*100, unchanged_count);
            pie_colors = [pie_colors; 0.7, 0.7, 0.7]; % Gray
        end
        
        % Create pie chart
        if ~isempty(pie_data)
            pie_handle = pie(pie_data);
            
            % Customize colors
            for j = 1:2:length(pie_handle)
                color_idx = ceil(j/2);
                if color_idx <= size(pie_colors, 1)
                    set(pie_handle(j), 'FaceColor', pie_colors(color_idx, :));
                end
            end
            
            % Add legend
            legend(pie_labels, 'Location', 'eastoutside', 'FontSize', 10);
            
            % Add title with statistics
            title_str = sprintf('Cluster Changes: %s (n=%d animals)\nTotal Clusters: %d | Threshold: ±%d%%', ...
                strrep(comparison_name, '_', ' vs '), length(data.common_animal_ids), ...
                length(data.cluster_ids), threshold*100);
            title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
            
            % Export if requested
            if do_export
                pie_filename = fullfile(export_folder, sprintf('cluster_proportions_pie_%s.png', comparison_name));
                saveas(fig, pie_filename);
                logger(sprintf('Saved pie chart: %s', pie_filename), 'INFO');
            end
        else
            logger(sprintf('Warning: No data to plot for %s', comparison_name), 'WARNING');
        end
    end
end

%% Create bar chart for comparative analysis
function create_cluster_proportion_bar_chart(diff_struct, threshold, visualize, export_folder, do_export)
    field_names = fieldnames(diff_struct);
    
    % Prepare data for bar chart
    comparison_names = {};
    gained_props = [];
    lost_props = [];
    unchanged_props = [];
    
    for i = 1:length(field_names)
        comparison_name = field_names{i};
        data = diff_struct.(comparison_name);
        
        comparison_names{i} = strrep(comparison_name, '_', ' vs ');
        gained_props(i) = data.proportion_gained * 100; % Convert to percentage
        lost_props(i) = data.proportion_lost * 100;
        unchanged_props(i) = data.proportion_unchanged * 100;
    end
    
    % Create stacked bar chart
    fig = figure('Visible', visualize);
    set(fig, 'Position', [100, 100, 1000, 600]);
    set(fig, 'Color', 'w');
    
    % Create stacked bar chart
    bar_data = [gained_props', lost_props', unchanged_props'];
    bar_handle = bar(bar_data, 'stacked');
    
    % Customize colors
    set(bar_handle(1), 'FaceColor', [1, 0, 0]); % Red for gained
    set(bar_handle(2), 'FaceColor', [0, 0, 1]); % Blue for lost
    set(bar_handle(3), 'FaceColor', [0.7, 0.7, 0.7]); % Gray for unchanged
    
    % Customize axes
    set(gca, 'XTickLabel', comparison_names);
    ylabel('Percentage of Clusters', 'FontSize', 12);
    xlabel('Comparison', 'FontSize', 12);
    title(sprintf('Cluster Proportion Analysis (Threshold: ±%d%%)', threshold*100), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Add legend
    legend({'Gained', 'Lost', 'Unchanged'}, 'Location', 'best');
    
    % Add percentage labels on bars
    for i = 1:length(comparison_names)
        if gained_props(i) > 5 % Only show label if segment is large enough
            text(i, gained_props(i)/2, sprintf('%.1f%%', gained_props(i)), ...
                'HorizontalAlignment', 'center', 'Color', 'white', 'FontWeight', 'bold');
        end
        if lost_props(i) > 5
            text(i, gained_props(i) + lost_props(i)/2, sprintf('%.1f%%', lost_props(i)), ...
                'HorizontalAlignment', 'center', 'Color', 'white', 'FontWeight', 'bold');
        end
        if unchanged_props(i) > 5
            text(i, gained_props(i) + lost_props(i) + unchanged_props(i)/2, sprintf('%.1f%%', unchanged_props(i)), ...
                'HorizontalAlignment', 'center', 'Color', 'black', 'FontWeight', 'bold');
        end
    end
    
    grid on;
    ylim([0, 100]);
    
    % Export if requested
    if do_export
        bar_filename = fullfile(export_folder, 'cluster_proportions_comparison_bar.png');
        saveas(fig, bar_filename);
        logger(sprintf('Saved bar chart: %s', bar_filename), 'INFO');
    end
end

%% Generate proportion visualizations
logger('Creating cluster proportion visualizations', 'INFO');

% Create pie charts for each comparison
create_cluster_proportion_pie_charts(diff_struct, threshold, visualize, export_folder, do_export);

% Create comparative bar chart
create_cluster_proportion_bar_chart(diff_struct, threshold, visualize, export_folder, do_export);

% Save proportion analysis table if export is enabled
if do_export
    proportion_filename = fullfile(export_folder, 'cluster_proportion_analysis_80pct.csv');
    writetable(proportion_analysis, proportion_filename);
    logger(sprintf('Saved proportion analysis table: %s', proportion_filename), 'INFO');
end

%% Individuality Analysis - Analyze cluster dominance by individual animals
logger('Starting individuality analysis - evaluating pose uniqueness per animal', 'INFO');

% Define individuality threshold (80% of frames in a cluster belong to one animal)
individuality_threshold = 0.8;

% Get all unique animal IDs from baseline condition
all_baseline_animal_ids = {};
for i = 1:length(frame_conditions)
    if strcmp(frame_conditions{i}, 'B') && ~strcmp(frame_animal_ids{i}, 'Unknown')
        all_baseline_animal_ids{end+1} = frame_animal_ids{i};
    end
end
unique_baseline_animals = unique(all_baseline_animal_ids);

logger(sprintf('Found %d unique animals with baseline condition for individuality analysis', ...
    length(unique_baseline_animals)), 'INFO');

% Initialize individuality analysis structure
individuality_analysis = struct();
individuality_analysis.threshold = individuality_threshold;
individuality_analysis.unique_animals = unique_baseline_animals;
individuality_analysis.cluster_ids = unique_clusters;

% For each cluster, calculate animal composition (only baseline frames)
cluster_animal_composition = cell(length(unique_clusters), 1);
cluster_dominant_animal = cell(length(unique_clusters), 1);
cluster_dominance_percentage = zeros(length(unique_clusters), 1);
cluster_is_individual = false(length(unique_clusters), 1);

logger('Analyzing cluster composition for each animal...', 'INFO');

for c_idx = 1:length(unique_clusters)
    cluster_id = unique_clusters(c_idx);
    cluster_frames = cluster_assignments == cluster_id;
    
    % Count frames per animal in this cluster (baseline only)
    animal_frame_counts = containers.Map();
    total_baseline_frames_in_cluster = 0;
    
    for i = 1:length(frame_conditions)
        if cluster_frames(i) && strcmp(frame_conditions{i}, 'B') && ~strcmp(frame_animal_ids{i}, 'Unknown')
            animal_id = frame_animal_ids{i};
            if isKey(animal_frame_counts, animal_id)
                animal_frame_counts(animal_id) = animal_frame_counts(animal_id) + 1;
            else
                animal_frame_counts(animal_id) = 1;
            end
            total_baseline_frames_in_cluster = total_baseline_frames_in_cluster + 1;
        end
    end
    
    % Find dominant animal and calculate percentage
    if total_baseline_frames_in_cluster > 0
        animals = keys(animal_frame_counts);
        counts = values(animal_frame_counts);
        counts = cell2mat(counts);
        
        [max_count, max_idx] = max(counts);
        dominant_animal = animals{max_idx};
        dominance_percentage = max_count / total_baseline_frames_in_cluster;
        
        cluster_dominant_animal{c_idx} = dominant_animal;
        cluster_dominance_percentage(c_idx) = dominance_percentage;
        cluster_is_individual(c_idx) = dominance_percentage >= individuality_threshold;
        
        % Store full composition
        composition = struct();
        for j = 1:length(animals)
            composition.(animals{j}) = counts(j) / total_baseline_frames_in_cluster;
        end
        cluster_animal_composition{c_idx} = composition;
    else
        cluster_dominant_animal{c_idx} = 'None';
        cluster_dominance_percentage(c_idx) = 0;
        cluster_is_individual(c_idx) = false;
        cluster_animal_composition{c_idx} = struct();
    end
end

% Store results in individuality analysis structure
individuality_analysis.cluster_animal_composition = cluster_animal_composition;
individuality_analysis.cluster_dominant_animal = cluster_dominant_animal;
individuality_analysis.cluster_dominance_percentage = cluster_dominance_percentage;
individuality_analysis.cluster_is_individual = cluster_is_individual;

% Calculate summary statistics
total_clusters_analyzed = sum(cluster_dominance_percentage > 0);
individual_clusters = sum(cluster_is_individual);
individual_percentage = (individual_clusters / total_clusters_analyzed) * 100;

logger(sprintf('Individuality Summary: %d/%d clusters (%.1f%%) show individual dominance (>%.0f%%)', ...
    individual_clusters, total_clusters_analyzed, individual_percentage, individuality_threshold*100), 'INFO');

% Per-animal individuality statistics
animal_individual_clusters = containers.Map();
animal_total_dominant_clusters = containers.Map();

for c_idx = 1:length(unique_clusters)
    if ~strcmp(cluster_dominant_animal{c_idx}, 'None')
        animal_id = cluster_dominant_animal{c_idx};
        
        % Count total dominant clusters per animal
        if isKey(animal_total_dominant_clusters, animal_id)
            animal_total_dominant_clusters(animal_id) = animal_total_dominant_clusters(animal_id) + 1;
        else
            animal_total_dominant_clusters(animal_id) = 1;
        end
        
        % Count individual clusters per animal
        if cluster_is_individual(c_idx)
            if isKey(animal_individual_clusters, animal_id)
                animal_individual_clusters(animal_id) = animal_individual_clusters(animal_id) + 1;
            else
                animal_individual_clusters(animal_id) = 1;
            end
        end
    end
end

% Store per-animal statistics
individuality_analysis.animal_individual_clusters = animal_individual_clusters;
individuality_analysis.animal_total_dominant_clusters = animal_total_dominant_clusters;

%% Create individuality visualization functions

function plot_individuality_tsne_map(analysisstruct, individuality_data, visualize)
    % Create figure showing individual clusters on t-SNE map
    fig = figure('Visible', visualize);
    set(fig, 'Position', [100, 100, 1400, 800]);
    set(fig, 'Color', 'w');
    
    % Subplot 1: Individual vs non-individual clusters
    subplot(1, 2, 1);
    plot(analysisstruct.zValues(:,1), analysisstruct.zValues(:,2), '.', ...
        'Color', [0.9, 0.9, 0.9], 'MarkerSize', 1);
    hold on;
    
    % Highlight individual clusters
    for c_idx = 1:length(individuality_data.cluster_ids)
        if individuality_data.cluster_is_individual(c_idx)
            cluster_id = individuality_data.cluster_ids(c_idx);
            cluster_frames = analysisstruct.annot_reordered{end,end} == cluster_id;
            if sum(cluster_frames) > 0
                plot(analysisstruct.zValues(cluster_frames,1), ...
                     analysisstruct.zValues(cluster_frames,2), '.', ...
                     'Color', [1, 0, 0], 'MarkerSize', 4);
            end
        end
    end
    
    title(sprintf('Individual Clusters (>%d%% dominance)', individuality_data.threshold*100));
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    axis equal;
    legend({'All clusters', 'Individual clusters'}, 'Location', 'best');
    grid on;
    
    % Subplot 2: Dominance percentage heatmap
    subplot(1, 2, 2);
    
    % Create a colormap based on dominance percentage
    dominance_values = zeros(size(analysisstruct.zValues, 1), 1);
    
    for c_idx = 1:length(individuality_data.cluster_ids)
        cluster_id = individuality_data.cluster_ids(c_idx);
        cluster_frames = analysisstruct.annot_reordered{end,end} == cluster_id;
        dominance_values(cluster_frames) = individuality_data.cluster_dominance_percentage(c_idx);
    end
    
    % Create scatter plot colored by dominance
    scatter(analysisstruct.zValues(:,1), analysisstruct.zValues(:,2), 8, dominance_values, 'filled');
    
    % Set colormap and colorbar
    colormap(hot);
    colorbar;
    caxis([0, 1]);
    
    title('Cluster Dominance Percentage');
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    axis equal;
    
    % Add main title
    sgtitle(sprintf('Pose Individuality Analysis (n=%d animals, %d clusters)', ...
        length(individuality_data.unique_animals), length(individuality_data.cluster_ids)), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

function plot_individuality_per_animal_tsne(analysisstruct, individuality_data, visualize)
    % Create figure showing individual clusters for each animal separately
    unique_animals = individuality_data.unique_animals;
    n_animals = length(unique_animals);
    
    % Create subplot grid
    n_cols = ceil(sqrt(n_animals));
    n_rows = ceil(n_animals / n_cols);
    
    fig = figure('Visible', visualize);
    set(fig, 'Position', [100, 100, 1600, 1200]);
    set(fig, 'Color', 'w');
    
    for a_idx = 1:n_animals
        animal_id = unique_animals{a_idx};
        
        subplot(n_rows, n_cols, a_idx);
        
        % Plot all points in gray
        plot(analysisstruct.zValues(:,1), analysisstruct.zValues(:,2), '.', ...
            'Color', [0.9, 0.9, 0.9], 'MarkerSize', 1);
        hold on;
        
        % Highlight clusters dominated by this animal
        for c_idx = 1:length(individuality_data.cluster_ids)
            if strcmp(individuality_data.cluster_dominant_animal{c_idx}, animal_id) && ...
               individuality_data.cluster_is_individual(c_idx)
                cluster_id = individuality_data.cluster_ids(c_idx);
                cluster_frames = analysisstruct.annot_reordered{end,end} == cluster_id;
                if sum(cluster_frames) > 0
                    plot(analysisstruct.zValues(cluster_frames,1), ...
                         analysisstruct.zValues(cluster_frames,2), '.', ...
                         'Color', [1, 0, 0], 'MarkerSize', 3);
                end
            end
        end
        
        % Calculate individual clusters for this animal
        individual_count = 0;
        if isKey(individuality_data.animal_individual_clusters, animal_id)
            individual_count = individuality_data.animal_individual_clusters(animal_id);
        end
        
        title(sprintf('%s (%d individual clusters)', animal_id, individual_count));
        xlabel('t-SNE 1');
        ylabel('t-SNE 2');
        axis equal;
        axis tight;
    end
    
    sgtitle('Individual Clusters per Animal', 'FontSize', 16, 'FontWeight', 'bold');
end

function plot_individuality_statistics(individuality_data, visualize, export_folder, do_export)
    % Create comprehensive statistical plots
    
    % Figure 1: Overall statistics
    fig1 = figure('Visible', visualize);
    set(fig1, 'Position', [100, 100, 1200, 800]);
    set(fig1, 'Color', 'w');
    
    % Subplot 1: Pie chart of individual vs non-individual clusters
    subplot(2, 2, 1);
    individual_count = sum(individuality_data.cluster_is_individual);
    total_count = length(individuality_data.cluster_ids);
    non_individual_count = total_count - individual_count;
    
    pie_data = [individual_count, non_individual_count];
    pie_labels = {sprintf('Individual (%.1f%%)', (individual_count/total_count)*100), ...
                  sprintf('Shared (%.1f%%)', (non_individual_count/total_count)*100)};
    pie_handle = pie(pie_data, pie_labels);
    
    % Color the pie slices
    set(pie_handle(1), 'FaceColor', [1, 0.2, 0.2]); % Red for individual
    set(pie_handle(3), 'FaceColor', [0.7, 0.7, 0.7]); % Gray for shared
    
    title(sprintf('Cluster Individuality (Threshold: %d%%)', individuality_data.threshold*100));
    
    % Subplot 2: Histogram of dominance percentages
    subplot(2, 2, 2);
    valid_dominance = individuality_data.cluster_dominance_percentage(individuality_data.cluster_dominance_percentage > 0);
    histogram(valid_dominance, 20, 'FaceColor', [0.3, 0.6, 1], 'EdgeColor', 'black');
    xlabel('Dominance Percentage');
    ylabel('Number of Clusters');
    title('Distribution of Cluster Dominance');
    xlim([0, 1]);
    
    % Add threshold line
    hold on;
    line([individuality_data.threshold, individuality_data.threshold], ylim, ...
        'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
    legend('Clusters', sprintf('%d%% Threshold', individuality_data.threshold*100), 'Location', 'best');
    
    % Subplot 3: Per-animal individual cluster counts
    subplot(2, 2, [3, 4]);
    
    % Prepare data for bar chart
    animals = individuality_data.unique_animals;
    individual_counts = zeros(length(animals), 1);
    total_dominant_counts = zeros(length(animals), 1);
    
    for i = 1:length(animals)
        animal_id = animals{i};
        if isKey(individuality_data.animal_individual_clusters, animal_id)
            individual_counts(i) = individuality_data.animal_individual_clusters(animal_id);
        end
        if isKey(individuality_data.animal_total_dominant_clusters, animal_id)
            total_dominant_counts(i) = individuality_data.animal_total_dominant_clusters(animal_id);
        end
    end
    
    % Create grouped bar chart
    bar_data = [individual_counts, total_dominant_counts - individual_counts];
    bar_handle = bar(bar_data, 'stacked');
    set(bar_handle(1), 'FaceColor', [1, 0.2, 0.2]); % Red for individual
    set(bar_handle(2), 'FaceColor', [0.7, 0.7, 0.7]); % Gray for shared
    
    xlabel('Animal ID');
    ylabel('Number of Clusters');
    title('Individual vs Shared Dominant Clusters per Animal');
    legend('Individual Clusters', 'Shared Dominant Clusters', 'Location', 'best');
    
    % Set x-axis labels
    set(gca, 'XTickLabel', animals);
    xtickangle(45);
    
    sgtitle('Pose Individuality Statistical Analysis', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Export if requested
    if do_export
        stats_filename = fullfile(export_folder, 'individuality_statistics.png');
        saveas(fig1, stats_filename);
        logger(sprintf('Saved individuality statistics: %s', stats_filename), 'INFO');
    end
end

%% Generate individuality visualizations
logger('Creating individuality visualizations', 'INFO');

% Create t-SNE individuality maps
plot_individuality_tsne_map(analysisstruct, individuality_analysis, visualize);

% Create per-animal individuality maps
plot_individuality_per_animal_tsne(analysisstruct, individuality_analysis, visualize);

% Create statistical analysis plots
plot_individuality_statistics(individuality_analysis, visualize, export_folder, do_export);

%% Create individuality summary table
logger('Creating individuality summary table', 'INFO');

% Create detailed summary table
individuality_table = table();
animal_ids = individuality_analysis.unique_animals;
individual_cluster_counts = zeros(length(animal_ids), 1);
total_dominant_cluster_counts = zeros(length(animal_ids), 1);
individuality_percentages = zeros(length(animal_ids), 1);

for i = 1:length(animal_ids)
    animal_id = animal_ids{i};
    
    % Get individual cluster count
    if isKey(individuality_analysis.animal_individual_clusters, animal_id)
        individual_cluster_counts(i) = individuality_analysis.animal_individual_clusters(animal_id);
    end
    
    % Get total dominant cluster count
    if isKey(individuality_analysis.animal_total_dominant_clusters, animal_id)
        total_dominant_cluster_counts(i) = individuality_analysis.animal_total_dominant_clusters(animal_id);
    end
    
    % Calculate individuality percentage
    if total_dominant_cluster_counts(i) > 0
        individuality_percentages(i) = (individual_cluster_counts(i) / total_dominant_cluster_counts(i)) * 100;
    end
end

individuality_table.Animal_ID = animal_ids;
individuality_table.Individual_Clusters = individual_cluster_counts;
individuality_table.Total_Dominant_Clusters = total_dominant_cluster_counts;
individuality_table.Individuality_Percentage = individuality_percentages;

% Add overall statistics
total_clusters_analyzed = length(individuality_analysis.cluster_ids);
total_individual_clusters = sum(individuality_analysis.cluster_is_individual);
overall_individuality_percentage = (total_individual_clusters / total_clusters_analyzed) * 100;

% Display results
disp('=== Pose Individuality Analysis Results ===');
fprintf('Threshold: %.0f%% dominance\n', individuality_analysis.threshold * 100);
fprintf('Total clusters analyzed: %d\n', total_clusters_analyzed);
fprintf('Individual clusters: %d (%.1f%%)\n', total_individual_clusters, overall_individuality_percentage);
fprintf('Shared clusters: %d (%.1f%%)\n', total_clusters_analyzed - total_individual_clusters, ...
    100 - overall_individuality_percentage);
disp(' ');
disp('Per-animal individuality:');
disp(individuality_table);

% Save individuality table if export is enabled
if do_export
    individuality_filename = fullfile(export_folder, 'pose_individuality_analysis.csv');
    writetable(individuality_table, individuality_filename);
    logger(sprintf('Saved individuality analysis table: %s', individuality_filename), 'INFO');
end

%% Save analysis results
analysis_filename = fullfile(export_folder, 'tsne_control_vs_baseline_analysis.mat');
%save(analysis_filename, 'diff_struct', 'summary_table', 'per_animal_analysis', 'proportion_analysis', 'conditions_to_compare', 'individuality_analysis');
%logger(sprintf('Saved control vs baseline analysis results: %s', analysis_filename), 'INFO');

logger('Control vs baseline t-SNE difference analysis completed', 'INFO');
