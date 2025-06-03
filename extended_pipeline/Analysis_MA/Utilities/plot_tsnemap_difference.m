%% Initialization
global logger GC
%logger('test tsne map difference', 'INFO'); % Modify
clear;
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Configuration for visualization and export
debugging = false;  % Set to true for debugging mode
if debugging
    visualize = 'on';  % Show figures during debugging
    do_export = false;  % Don't export during debugging
else
    visualize = 'off';  % Don't show figures in production mode
    do_export = true;   % Export figures in production mode
end

% Export folder
export_folder = fullfile(GC.temp_root, 'figs_presentation_painAI');
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

%% Calculate difference between conditions
% select conditions
conditions_to_compare = {{'S'}, {'F'}; {'H'}, {'N'}}; % we start with these combinations(S vs F and H vs N), but it could be expanded later to 'B', 'S', 'F' ; and  'B', 'H', 'N';
% Initialize a structure to hold the differences
diff_struct = struct();

%% Extract cluster information for each comparison
logger('Analyzing cluster differences between conditions', 'INFO');

% Get cluster assignments for all frames
cluster_assignments = analysisstruct.annot_reordered{end,end}; % Final cluster assignments
unique_clusters = unique(cluster_assignments);
unique_clusters = unique_clusters(unique_clusters > 0); % Remove background/noise clusters

% Map cluster assignments to animal conditions
frame_conditions = cell(length(cluster_assignments), 1);
for i = 1:length(cluster_assignments)
    if i <= length(animal_list_used_after_analysis)
        this_a = animal_list_used_after_analysis{i};
        this_cond = this_a(end); % Get last character as condition
        frame_conditions{i} = this_cond;
    else
        frame_conditions{i} = 'Unknown';
    end
end

%% Process each condition comparison
for comp_idx = 1:size(conditions_to_compare, 1)
    control_conditions = conditions_to_compare{comp_idx, 1};
    experimental_conditions = conditions_to_compare{comp_idx, 2};
    
    comparison_name = sprintf('%s_vs_%s', strjoin(experimental_conditions, ''), strjoin(control_conditions, ''));
    logger(sprintf('Processing comparison: %s', comparison_name), 'INFO');
    
    % Find frames belonging to each condition
    control_frames = false(length(frame_conditions), 1);
    experimental_frames = false(length(frame_conditions), 1);
    
    for i = 1:length(frame_conditions)
        if any(strcmp(frame_conditions{i}, control_conditions))
            control_frames(i) = true;
        elseif any(strcmp(frame_conditions{i}, experimental_conditions))
            experimental_frames(i) = true;
        end
    end
    
    % Calculate cluster prevalence for each condition
    control_cluster_counts = zeros(length(unique_clusters), 1);
    experimental_cluster_counts = zeros(length(unique_clusters), 1);
    
    for c_idx = 1:length(unique_clusters)
        cluster_id = unique_clusters(c_idx);
        cluster_frames = cluster_assignments == cluster_id;
        
        control_cluster_counts(c_idx) = sum(cluster_frames' & control_frames);
        experimental_cluster_counts(c_idx) = sum(cluster_frames' & experimental_frames);
    end
    
    % Calculate total frames for normalization
    total_control_frames = sum(control_frames);
    total_experimental_frames = sum(experimental_frames);
    
    % Calculate normalized prevalence (percentage of frames in each condition)
    control_prevalence = control_cluster_counts / max(total_control_frames, 1);
    experimental_prevalence = experimental_cluster_counts / max(total_experimental_frames, 1);
    
    % Calculate difference (experimental - control)
    cluster_difference = experimental_prevalence - control_prevalence;
    
    % Identify clusters that are gained (higher in experimental) or lost (higher in control)
    gained_clusters = unique_clusters(cluster_difference > 0);
    lost_clusters = unique_clusters(cluster_difference < 0);
    
    % Store results
    diff_struct.(comparison_name).control_conditions = control_conditions;
    diff_struct.(comparison_name).experimental_conditions = experimental_conditions;
    diff_struct.(comparison_name).cluster_ids = unique_clusters;
    diff_struct.(comparison_name).control_prevalence = control_prevalence;
    diff_struct.(comparison_name).experimental_prevalence = experimental_prevalence;
    diff_struct.(comparison_name).cluster_difference = cluster_difference;
    diff_struct.(comparison_name).gained_clusters = gained_clusters;
    diff_struct.(comparison_name).lost_clusters = lost_clusters;
    
    logger(sprintf('Found %d gained clusters and %d lost clusters for %s', ...
        length(gained_clusters), length(lost_clusters), comparison_name), 'INFO');
end

%% Create visualization function for t-SNE difference maps
function plot_tsne_difference_map(analysisstruct, diff_data, comparison_name, visualize)
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
        strjoin(diff_data.experimental_conditions, ','), ...
        strjoin(diff_data.control_conditions, ',')));
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
        strjoin(diff_data.experimental_conditions, ','), ...
        strjoin(diff_data.control_conditions, ',')));
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
    legend({'Background', 'Gained (Experimental > Control)', 'Lost (Experimental < Control)'}, ...
        'Location', 'best');
    
    title(sprintf('Combined Difference Map: %s', comparison_name));
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    axis equal;
    grid on;
    
    % Add main title
    sgtitle(sprintf('t-SNE Cluster Differences: %s', strrep(comparison_name, '_', ' vs ')), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

%% Create enhanced watershed-based visualization
function plot_tsne_watershed_difference(analysisstruct, diff_data, comparison_name, visualize)
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
    
    title('Difference Map with Cluster Boundaries');
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
    
    title('Cluster Difference Intensity');
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    axis equal;
    
    % Add main title
    sgtitle(sprintf('Enhanced t-SNE Analysis: %s', strrep(comparison_name, '_', ' vs ')), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

%% Generate visualizations for each comparison
for comp_idx = 1:size(conditions_to_compare, 1)
    control_conditions = conditions_to_compare{comp_idx, 1};
    experimental_conditions = conditions_to_compare{comp_idx, 2};
    comparison_name = sprintf('%s_vs_%s', strjoin(experimental_conditions, ''), strjoin(control_conditions, ''));
    
    plot_tsne_difference_map(analysisstruct, diff_struct.(comparison_name), ...
        comparison_name, visualize);
    
    plot_tsne_watershed_difference(analysisstruct, diff_struct.(comparison_name), ...
        comparison_name, visualize);
end

%% Create summary statistics table
logger('Creating summary statistics', 'INFO');

% Create a summary table of differences
summary_table = table();
comparison_names = {};
num_gained = [];
num_lost = [];
total_clusters = [];
max_gain = [];
max_loss = [];
top_20_increased_clusters = {};
bottom_20_decreased_clusters = {};
num_top_20_increased = [];
num_bottom_20_decreased = [];

comp_idx = 1;
field_names = fieldnames(diff_struct);
for i = 1:length(field_names)
    comparison_name = field_names{i};
    data = diff_struct.(comparison_name);
    
    comparison_names{comp_idx} = strrep(comparison_name, '_', ' vs ');
    num_gained(comp_idx) = length(data.gained_clusters);
    num_lost(comp_idx) = length(data.lost_clusters);
    total_clusters(comp_idx) = length(data.cluster_ids);
    max_gain(comp_idx) = max(data.cluster_difference);
    max_loss(comp_idx) = min(data.cluster_difference);
    
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
      % Log the statistically significant clusters with formatted output
    % cluster_z_pairs_increased = sprintf('Cluster %d (z=%.2f) ', ...
        % [top_20_increased_cluster_ids'; top_20_increased_zscores']);
    % logger(sprintf('For %s: Top 20%% increased clusters (n=%d): %s', ...
        % comparison_names{comp_idx}, length(top_20_increased_cluster_ids), ...
        % cluster_z_pairs_increased), 'INFO');
    
    % cluster_z_pairs_decreased = sprintf('Cluster %d (z=%.2f) ', ...
    %     [bottom_20_decreased_cluster_ids'; bottom_20_decreased_zscores']);
    % logger(sprintf('For %s: Bottom 20%% decreased clusters (n=%d): %s', ...
    %     comparison_names{comp_idx}, length(bottom_20_decreased_cluster_ids), ...
    %     cluster_z_pairs_decreased), 'INFO');
    
    comp_idx = comp_idx + 1;
end

summary_table.Comparison = comparison_names';
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
disp('=== t-SNE Cluster Difference Summary ===');
disp(summary_table);

% Save summary table if export is enabled
if do_export
    summary_filename = fullfile(export_folder, 'cluster_difference_summary.csv');
    writetable(summary_table, summary_filename);
    logger(sprintf('Saved summary table: %s', summary_filename), 'INFO');
end

%% Save analysis results
analysis_filename = fullfile(export_folder, 'tsne_difference_analysis.mat');
save(analysis_filename, 'diff_struct', 'summary_table', 'conditions_to_compare');
logger(sprintf('Saved analysis results: %s', analysis_filename), 'INFO');

logger('t-SNE difference analysis completed', 'INFO');