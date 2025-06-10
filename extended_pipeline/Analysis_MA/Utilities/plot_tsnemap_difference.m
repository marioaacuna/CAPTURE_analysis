%% Initialization
global logger GC
%logger('test tsne map difference', 'INFO'); % Modify
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
    total_cluster_frames = zeros(length(unique_clusters), 1);
    
    for c_idx = 1:length(unique_clusters)
        cluster_id = unique_clusters(c_idx);
        cluster_frames = cluster_assignments == cluster_id;
        
        control_cluster_counts(c_idx) = sum(cluster_frames' & control_frames);
        experimental_cluster_counts(c_idx) = sum(cluster_frames' & experimental_frames);
        total_cluster_frames(c_idx) = sum(cluster_frames); % Total frames in this cluster
    end
    
    % Calculate total frames for normalization (for reference)
    total_control_frames = sum(control_frames);
    total_experimental_frames = sum(experimental_frames);
    
    % Calculate cluster-normalized prevalence (percentage within each cluster)
    control_prevalence = zeros(length(unique_clusters), 1);
    experimental_prevalence = zeros(length(unique_clusters), 1);
    
    for c_idx = 1:length(unique_clusters)
        if total_cluster_frames(c_idx) > 0
            control_prevalence(c_idx) = control_cluster_counts(c_idx) / total_cluster_frames(c_idx);
            experimental_prevalence(c_idx) = experimental_cluster_counts(c_idx) / total_cluster_frames(c_idx);
        end
    end
    
    % Calculate difference (experimental - control) normalized by cluster size
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
            title_str = sprintf('Cluster Changes: %s\nTotal Clusters: %d | Threshold: ±%d%%', ...
                strrep(comparison_name, '_', ' vs '), ...
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

%% Save analysis results
analysis_filename = fullfile(export_folder, 'tsne_difference_analysis.mat');
save(analysis_filename, 'diff_struct', 'summary_table', 'proportion_analysis', 'conditions_to_compare');
logger(sprintf('Saved analysis results: %s', analysis_filename), 'INFO');

logger('t-SNE difference analysis completed', 'INFO');