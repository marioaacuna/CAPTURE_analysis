%% ========================================================================
%                   t-SNE POSE INDIVIDUALITY ANALYSIS (ALL CONDITIONS)
% ========================================================================
%
% DESCRIPTION:
% This script performs comprehensive analysis of behavioral data from t-SNE 
% clustered pose sequences across different experimental conditions, with a 
% focus on evaluating pose individuality across animals. It generates multiple 
% visualizations and quantitative metrics to understand how behavioral clusters 
% are distributed among individual animals and identifies unique vs. shared 
% behavioral patterns.
%
% MAIN OBJECTIVES:
% 1. Visualize t-SNE behavioral landscapes for each animal in selected condition
% 2. Analyze cluster composition and dominance patterns across animals
% 3. Quantify pose individuality - behaviors unique to specific animals
% 4. Generate statistical summaries and export results for further analysis
%
% KEY CONCEPTS:
% - Individual Clusters: Behavioral clusters where ≥80% of frames belong to one animal
% - Shared Clusters: Behavioral clusters with mixed animal representation
% - Dominance Percentage: The proportion of frames in a cluster belonging to the most frequent animal
% - Pose Individuality: The extent to which an animal exhibits unique behavioral patterns
%
% INPUT REQUIREMENTS:
% - analysisstruct: CAPTURE analysis structure with t-SNE embeddings and cluster assignments
% - animal_condition_identifier: Cell array identifying animal and condition for each frame
% - GC (general_configs): Configuration structure with file paths and parameters
%
% ANALYSIS WORKFLOW:
% 1. Data Loading & Preprocessing
%    - Load t-SNE coordinates, cluster assignments, and animal identifiers
%    - User selects condition to analyze (B, S, F, H, N) via dialog
%    - Filter for selected condition frames only
%    - Extract unique animal IDs and generate color mappings
%    - IMPORTANT: Uses GLOBAL t-SNE coordinates as reference space for consistency
%
% 2. Basic Visualization Generation
%    - Density maps: Kernel density estimation for each animal's behavioral space
%    - Scatter plots: Combined and individual animal t-SNE projections
%    - Color-coded visualization showing spatial distribution of behaviors
%    - All visualizations use global t-SNE space with condition data highlighted
%
% 3. Individuality Analysis
%    - Cluster Composition Analysis: Calculate animal+condition representation in each cluster
%    - Dominance Calculation: Identify dominant animal+condition per cluster (highest frame count)
%    - Individual Threshold Application: Mark clusters with ≥80% single animal+condition dominance
%    - Per-Animal+Condition Statistics: Count individual vs shared clusters for each combination
%    - CRITICAL: Uses animal+condition identifiers (e.g., "Animal1_S" vs "Animal1_F") 
%    - This prevents treating same animals in different conditions as identical entities
%    - IMPORTANT: Uses GLOBAL data (all animals, all conditions) for individuality calculation
%    - This prevents inflation of individuality scores in sparse conditions
%
% 4. Advanced Visualizations
%    - Individual Cluster Maps: Highlight unique behaviors per animal on global t-SNE space
%    - Dominance Heatmaps: Color-code points by dominance percentage
%    - Statistical Plots: Pie charts, histograms, and bar charts of individuality metrics
%    - Global context maintained for consistent interpretation across conditions
%
% 5. Quantitative Output
%    - Summary tables with per-animal individuality percentages
%    - Overall dataset statistics (total individual vs. shared clusters)
%    - Export results to CSV files for external analysis
%
% OUTPUT METRICS:
% Cluster-Level Metrics:
% - Individual_Clusters_Global: Number of clusters where animal+condition shows ≥80% dominance (globally)
% - Total_Clusters_Present_Global: Total clusters where animal+condition appears globally (any percentage)
% - Total_Dominant_Clusters_Global: Clusters where animal+condition is the most frequent globally
% - Individuality_Percentage_Global: (Individual_Clusters_Global / Total_Clusters_Present_Global) × 100
%
% Frame-Level Metrics:
% - Individual_Frames_Global: Number of frames spent in individual behaviors (globally)
% - Total_Frames_Global: Total number of frames for this animal+condition combination (globally)
% - Frame_Individuality_Percentage_Global: (Individual_Frames_Global / Total_Frames_Global) × 100
%
% BIOLOGICAL INTERPRETATION:
% High individuality percentage indicates an animal+condition combination has many unique 
% behavioral patterns not commonly exhibited by other animal+condition combinations. 
% Low values suggest more shared/common behavioral repertoires. 
%
% Cluster vs Frame Individuality:
% - Cluster Individuality: Measures how many distinct behavioral types are unique to an animal
% - Frame Individuality: Measures how much time an animal spends in unique behaviors
% - An animal might have few unique behavioral types but spend most of their time in them (high frame, low cluster)
% - Or have many unique behaviors but use them rarely (high cluster, low frame)
%
% This dual analysis helps identify:
% - Animal+condition combinations with distinctive behavioral signatures
% - Whether uniqueness comes from rare specialized behaviors or frequent unique patterns
% - Common vs. rare behavioral patterns in the population
% - Behavioral diversity and specialization within and across conditions
% - How the same animal's behavior changes across different experimental conditions
%
% IMPORTANT: Individuality is calculated using the GLOBAL dataset (all animals, all 
% conditions) as reference, ensuring unbiased interpretation across different 
% experimental conditions. This prevents inflation of individuality scores in 
% conditions with fewer animals or limited behavioral repertoires. The analysis 
% treats each animal+condition combination as a unique entity (e.g., "Animal1_S" 
% vs "Animal1_F"), which is critical for understanding how the same animal's 
% behavior differs across experimental conditions. While visualizations focus on 
% the selected condition, the individuality metrics reflect true behavioral 
% uniqueness in the context of the complete experiment.
%
% CONFIGURATION OPTIONS:
% - individuality_threshold: Dominance threshold for individual classification (default: 0.8)
% - debugging: Toggle between development (visualize='on') and production mode
% - do_export: Enable/disable figure and data export functionality
%
% GENERATED OUTPUTS:
% Visual Outputs:
% - Density maps per animal (PDF)
% - Combined scatter plot of all animals (PDF)  
% - Individual scatter plots per animal (PDF)
% - t-SNE individuality analysis maps (PNG)
% - Per-animal individual cluster visualizations (PNG)
% - Statistical analysis plots (PNG)
%
% Data Outputs:
% - [condition]_pose_individuality_analysis.csv: Detailed per-animal metrics
% - Console logging: Real-time analysis progress and summary statistics
% - MATLAB workspace: Complete individuality_analysis structure
%
% DEPENDENCIES:
% - CAPTURE analysis framework
% - MotionMapper t-SNE implementation
% - plotdensitymaps function for density visualization
% - distinguishable_colors function (optional, falls back to HSV)
% - Statistics and Machine Learning Toolbox (for clustering functions)
%
% USAGE EXAMPLE:
% % Ensure general_configs is properly set up
% GC = general_configs;
% % Run the complete analysis (user will be prompted to select condition)
% run('plot_pose_individuality_analysis.m');
%
% AUTHORS: CAPTURE Analysis Team
% CREATED: 2025
% LAST MODIFIED: June 2025 - Added comprehensive pose individuality analysis
%
% RELATED SCRIPTS:
% - plot_tsnemap_difference_controls.m: Comparative analysis between conditions
% - compute_tsne_features.m: Feature extraction for t-SNE analysis
% - cluster_tsne_maps.m: Clustering of behavioral sequences
%
% ========================================================================

%% Initialization
clc;
logger('Starting pose individuality analysis script', 'INFO');
clear;
close all;

GC = general_configs;
rootpath = GC.preprocessing_rootpath;

%% User Input Dialog for Condition and Debug Mode Selection
% Available conditions
available_conditions = {'B', 'S', 'F', 'H', 'N'};
condition_names = {'Baseline', 'Saline', 'Formalin', 'Sham', 'Neuropathic'};

% Create condition selection dialog
[condition_idx, ok] = listdlg('PromptString', 'Select condition to analyze:', ...
                              'SelectionMode', 'single', ...
                              'ListString', strcat(available_conditions, ' - ', condition_names), ...
                              'Name', 'Condition Selection', ...
                              'ListSize', [300, 150]);

if ~ok
    logger('Analysis cancelled by user', 'WARNING');
    return;
end

selected_condition = available_conditions{condition_idx};
selected_condition_name = condition_names{condition_idx};

% Debug mode selection
debug_choice = questdlg('Select analysis mode:', 'Debug Mode Selection', ...
                       'Debug Mode (Show figures, no export)', ...
                       'Production Mode (Export figures, no display)', ...
                       'Debug Mode (Show figures, no export)');

% Configuration for visualization and export based on user choice
if contains(debug_choice, 'Debug')
    debugging = true;
    visualize = 'on';   % Show figures during debugging
    do_export = false;  % Don't export during debugging
    logger(sprintf('Running in DEBUG MODE for condition %s (%s)', selected_condition, selected_condition_name), 'INFO');
else
    debugging = false;
    visualize = 'off';  % Don't show figures in production mode
    do_export = true;   % Export figures in production mode
    logger(sprintf('Running in PRODUCTION MODE for condition %s (%s)', selected_condition, selected_condition_name), 'INFO');
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

%%
% Extract animal list and conditions
upsamplig_factor = GC.repfactor;
long_animal_frames_identifier = repelem(animal_condition_identifier,upsamplig_factor);
animal_list_used_after_analysis = long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});

% Plot t-SNE maps for selected condition per animal
logger(sprintf('Plotting t-SNE maps for %s condition per animal', selected_condition_name), 'INFO');

%% Extract condition identifiers and animal names
% Extract last character of each identifier to determine condition
conditions = cellfun(@(x) x(end), animal_list_used_after_analysis, 'UniformOutput', false);

% Find frames for the selected condition
condition_mask = strcmp(conditions, selected_condition);
logger(sprintf('Found %d frames for condition %s (%s)', sum(condition_mask), selected_condition, selected_condition_name), 'INFO');

% Extract animal names (everything before the underscore)
animal_names = cellfun(@(x) x(1:find(x=='_',1)-1), animal_list_used_after_analysis, 'UniformOutput', false);

% Filter for selected condition only
idx_condition = strcmp(conditions, selected_condition);
condition_animal_names = animal_names(idx_condition);
condition_zvals = analysisstruct.zValues(idx_condition, :);

% Get unique animal names for selected condition
unique_animals = unique(condition_animal_names);
num_animals = length(unique_animals);

logger(sprintf('Found %d unique animals in %s condition', num_animals, selected_condition_name), 'INFO');

% Generate distinct colors for each animal
% Add path to distinguishable_colors if needed
if ~exist('distinguishable_colors', 'file')
    addpath(fullfile(GC.preprocessing_rootpath, '3rd party toolboxes', '_plotting'));
end

try
    colors = distinguishable_colors(num_animals);
catch
    % Fallback: use hsv colormap if distinguishable_colors is not available
    logger('Using HSV colormap as fallback for animal colors', 'WARNING');
    colors = hsv(num_animals);
end

% Create color map for animals
animal_color_map = containers.Map();
for i = 1:num_animals
    animal_color_map(unique_animals{i}) = colors(i, :);
end

% Store all z-values for easier access
zvals = analysisstruct.zValues;

%% Create density maps for each animal in selected condition
logger(sprintf('Creating density maps for each animal in %s condition', selected_condition_name), 'INFO');
logger('NOTE: Using global t-SNE coordinate space for consistent density estimation', 'INFO');

% Calculate subplot layout
n_cols = ceil(sqrt(num_animals));
n_rows = ceil(num_animals / n_cols);

fig_density = figure('Name', sprintf('Density Maps: %s per Animal (Global Reference)', selected_condition_name), 'Color', 'w', ...
    'Position', [100, 100, 300*n_cols, 300*n_rows], 'Visible', visualize);

% Use global t-SNE parameters for consistent density estimation
global_zvals = analysisstruct.zValues; % Global coordinates for reference

for i = 1:num_animals
    animal = unique_animals{i};
    idx_animal = strcmp(condition_animal_names, animal);
    
    subplot(n_rows, n_cols, i);
    h_ax = gca;
    set(h_ax, 'Color', 'w');
    
    if sum(idx_animal) > 0  % Check if animal has data points
        % Use global coordinate range for density map limits
        plotdensitymaps({condition_zvals(idx_animal,:)}, 1, h_ax, analysisstruct.params.density_width, ...
            max(global_zvals(:))*analysisstruct.params.expansion_factor, analysisstruct.params.density_res);
    end
    
    title(['Animal: ' animal]);
    axis square
    
    logger(['Processed density map for animal: ' animal ' (' num2str(sum(idx_animal)) ' frames)'], 'INFO');
end

% Export the density figure if needed
if do_export
    logger(sprintf('Exporting %s density maps per animal to: %s', selected_condition_name, export_folder), 'INFO');
    export_name = sprintf('density_map_%s_per_animal', selected_condition);
    exportgraphics(fig_density, fullfile(export_folder, [export_name '.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
end

%% Create combined scatter plot with all animals in selected condition
logger(sprintf('Creating combined scatter plot for all animals in %s condition', selected_condition_name), 'INFO');
logger('NOTE: Showing condition data within global t-SNE coordinate space', 'INFO');

fig_scatter = figure('Name', sprintf('Scatter Plot: %s All Animals (Global Reference)', selected_condition_name), 'Color', 'w', 'Visible', visualize);
hold on;

% Plot global background data in very light gray
global_zvals = analysisstruct.zValues;
plot(global_zvals(:,1), global_zvals(:,2), '.', 'Color', [0.95, 0.95, 0.95], 'MarkerSize', 1);

% Plot each animal with its specific color
for i = 1:num_animals
    animal = unique_animals{i};
    idx_animal = strcmp(condition_animal_names, animal);
    
    if sum(idx_animal) > 0  % Check if animal has data points
        scatter(condition_zvals(idx_animal,1), condition_zvals(idx_animal,2), 15, ...
            animal_color_map(animal), 'Marker', '.', 'DisplayName', animal);
        
        logger(['Added scatter points for animal: ' animal ' (' num2str(sum(idx_animal)) ' frames)'], 'INFO');
    end
end

hold off;
legend('Location', 'best');
title(sprintf('%s Condition t-SNE Map - All Animals (Global Reference)', selected_condition_name));
xlabel('t-SNE Dimension 1 (Global Space)');
ylabel('t-SNE Dimension 2 (Global Space)');
axis equal tight;

% Export the scatter figure if needed
if do_export
    logger(sprintf('Exporting %s scatter plot per animal to: %s', selected_condition_name, export_folder), 'INFO');
    export_name = sprintf('scatter_plot_%s_per_animal', selected_condition);
    exportgraphics(fig_scatter, fullfile(export_folder, [export_name '.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
end

%% Create individual scatter plots for each animal (optional)
logger('Creating individual scatter plots for each animal', 'INFO');

% Calculate subplot layout for individual plots
fig_individual = figure('Name', sprintf('Individual Scatter Plots: %s per Animal (Global Reference)', selected_condition_name), 'Color', 'w', ...
    'Position', [200, 200, 300*n_cols, 300*n_rows], 'Visible', visualize);

for i = 1:num_animals
    animal = unique_animals{i};
    idx_animal = strcmp(condition_animal_names, animal);
    
    subplot(n_rows, n_cols, i);
    
    % Plot global background in very light gray
    plot(global_zvals(:,1), global_zvals(:,2), '.', 'Color', [0.95, 0.95, 0.95], 'MarkerSize', 0.5);
    hold on;
    
    if sum(idx_animal) > 0  % Check if animal has data points
        scatter(condition_zvals(idx_animal,1), condition_zvals(idx_animal,2), 15, ...
            animal_color_map(animal), 'Marker', '.');
    end
    
    title(['Animal: ' animal]);
    xlabel('t-SNE Dimension 1 (Global)');
    ylabel('t-SNE Dimension 2 (Global)');
    axis equal tight;
end

% Export the individual plots figure if needed
if do_export
    logger(sprintf('Exporting individual %s scatter plots per animal to: %s', selected_condition_name, export_folder), 'INFO');
    export_name = sprintf('scatter_plots_%s_individual_animals', selected_condition);
    exportgraphics(fig_individual, fullfile(export_folder, [export_name '.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
end

%% Individuality Analysis - Analyze cluster dominance by individual animals using GLOBAL data
logger('Starting individuality analysis - evaluating pose uniqueness per animal', 'INFO');
logger('NOTE: Using GLOBAL data (all conditions, all animals) for robust individuality assessment', 'INFO');

% Define individuality threshold (80% of frames in a cluster belong to one animal)
individuality_threshold = 0.8;

% Use GLOBAL cluster assignments and animal+condition identifiers for individuality calculation
% This ensures unbiased assessment across all experimental conditions
% CRITICAL: Use full animal+condition identifiers to distinguish same animals in different conditions
cluster_assignments = analysisstruct.annot_reordered{end,end}; % Final cluster assignments (ALL FRAMES)
global_animal_condition_ids = animal_list_used_after_analysis; % Full animal+condition identifiers (e.g., "Animal1_S", "Animal1_F")
unique_clusters = unique(cluster_assignments);
unique_clusters = unique_clusters(unique_clusters > 0); % Remove background/noise clusters

logger(sprintf('Found %d unique clusters in GLOBAL dataset for individuality analysis', length(unique_clusters)), 'INFO');
logger(sprintf('Using global dataset: %d total frames, %d unique animal+condition combinations across all conditions', ...
    length(cluster_assignments), length(unique(global_animal_condition_ids))), 'INFO');

% Get unique animal+condition combinations globally
global_unique_animal_condition_ids = unique(global_animal_condition_ids);
% Also get unique animal names globally (for some legacy compatibility)
global_animal_names = cellfun(@(x) x(1:find(x=='_',1)-1), global_animal_condition_ids, 'UniformOutput', false);
global_unique_animals = unique(global_animal_names);
condition_unique_animals = unique_animals; % Animals in selected condition

logger(sprintf('Selected condition (%s) contains %d/%d animals from global dataset', ...
    selected_condition_name, length(condition_unique_animals), length(global_unique_animals)), 'INFO');

% Initialize individuality analysis structure
individuality_analysis = struct();
individuality_analysis.threshold = individuality_threshold;
individuality_analysis.unique_animals = condition_unique_animals; % Animals in selected condition
individuality_analysis.global_unique_animals = global_unique_animals; % All animals globally
individuality_analysis.cluster_ids = unique_clusters;

% For each cluster, calculate animal composition using GLOBAL dataset (all conditions)
cluster_animal_composition = cell(length(unique_clusters), 1);
cluster_dominant_animal = cell(length(unique_clusters), 1);
cluster_dominance_percentage = zeros(length(unique_clusters), 1);
cluster_is_individual = false(length(unique_clusters), 1);

logger('Analyzing cluster composition using GLOBAL dataset (all conditions, all animals)...', 'INFO');

for c_idx = 1:length(unique_clusters)
    cluster_id = unique_clusters(c_idx);
    cluster_frames = cluster_assignments == cluster_id; % Use GLOBAL cluster assignments
    
    % Count frames per animal+condition combination in this cluster (GLOBAL dataset)
    animal_condition_frame_counts = containers.Map();
    total_frames_in_cluster = sum(cluster_frames);
    
    for i = 1:length(global_animal_condition_ids)
        if cluster_frames(i)
            animal_condition_id = global_animal_condition_ids{i}; % Full ID: "Animal1_S", "Animal1_F", etc.
            if isKey(animal_condition_frame_counts, animal_condition_id)
                animal_condition_frame_counts(animal_condition_id) = animal_condition_frame_counts(animal_condition_id) + 1;
            else
                animal_condition_frame_counts(animal_condition_id) = 1;
            end
        end
    end
    
    % Find dominant animal+condition combination and calculate percentage
    if total_frames_in_cluster > 0 && ~isempty(keys(animal_condition_frame_counts))
        animal_condition_ids = keys(animal_condition_frame_counts);
        counts = values(animal_condition_frame_counts);
        counts = cell2mat(counts);
        
        [max_count, max_idx] = max(counts);
        dominant_animal_condition = animal_condition_ids{max_idx};
        dominance_percentage = max_count / total_frames_in_cluster;
        
        cluster_dominant_animal{c_idx} = dominant_animal_condition; % Store full animal+condition ID
        cluster_dominance_percentage(c_idx) = dominance_percentage;
        cluster_is_individual(c_idx) = dominance_percentage >= individuality_threshold;
        
        % Store full composition with animal+condition identifiers
        composition = struct();
        for j = 1:length(animal_condition_ids)
            % Create valid field name by prefixing with 'ID_' and replacing any special chars
            valid_field_name = ['ID_' regexprep(animal_condition_ids{j}, '[^a-zA-Z0-9_]', '_')];
            composition.(valid_field_name) = counts(j) / total_frames_in_cluster;
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

% Calculate summary statistics based on GLOBAL individuality analysis
total_clusters_analyzed = sum(cluster_dominance_percentage > 0);
individual_clusters = sum(cluster_is_individual);
individual_percentage = (individual_clusters / total_clusters_analyzed) * 100;

logger(sprintf('GLOBAL Individuality Summary: %d/%d clusters (%.1f%%) show individual dominance (>%.0f%%)', ...
    individual_clusters, total_clusters_analyzed, individual_percentage, individuality_threshold*100), 'INFO');
logger(sprintf('This analysis includes %d total animal+condition combinations across all experimental conditions', ...
    length(global_unique_animal_condition_ids)), 'INFO');

% Per-animal+condition individuality statistics using GLOBAL data
animal_condition_individual_clusters = containers.Map();
animal_condition_total_clusters_present = containers.Map();
animal_condition_total_dominant_clusters = containers.Map();

% First, count all clusters where each animal+condition combination appears globally (not just dominant)
for c_idx = 1:length(unique_clusters)
    if ~isempty(cluster_animal_composition{c_idx})
        composition = cluster_animal_composition{c_idx};
        field_names = fieldnames(composition);
        
        % Check each animal+condition combination in this cluster
        for f_idx = 1:length(field_names)
            field_name = field_names{f_idx};
            % Extract animal+condition ID by removing 'ID_' prefix
            animal_condition_id = field_name(4:end); % Remove 'ID_' prefix
            % Convert back any replaced characters
            animal_condition_id = regexprep(animal_condition_id, '__', '_');
            
            % Count total clusters where this animal+condition combination appears globally
            if isKey(animal_condition_total_clusters_present, animal_condition_id)
                animal_condition_total_clusters_present(animal_condition_id) = animal_condition_total_clusters_present(animal_condition_id) + 1;
            else
                animal_condition_total_clusters_present(animal_condition_id) = 1;
            end
        end
    end
end

% Count individual and dominant clusters per animal+condition combination (global analysis)
for c_idx = 1:length(unique_clusters)
    if ~strcmp(cluster_dominant_animal{c_idx}, 'None')
        animal_condition_id = cluster_dominant_animal{c_idx}; % Full animal+condition ID
        
        % Count total dominant clusters per animal+condition combination globally
        if isKey(animal_condition_total_dominant_clusters, animal_condition_id)
            animal_condition_total_dominant_clusters(animal_condition_id) = animal_condition_total_dominant_clusters(animal_condition_id) + 1;
        else
            animal_condition_total_dominant_clusters(animal_condition_id) = 1;
        end
        
        % Count individual clusters per animal+condition combination globally
        if cluster_is_individual(c_idx)
            if isKey(animal_condition_individual_clusters, animal_condition_id)
                animal_condition_individual_clusters(animal_condition_id) = animal_condition_individual_clusters(animal_condition_id) + 1;
            else
                animal_condition_individual_clusters(animal_condition_id) = 1;
            end
        end
    end
end

% Store per-animal+condition statistics
individuality_analysis.animal_individual_clusters = animal_condition_individual_clusters;
individuality_analysis.animal_total_clusters_present = animal_condition_total_clusters_present;
individuality_analysis.animal_total_dominant_clusters = animal_condition_total_dominant_clusters;

%% Frame-level individuality analysis - Calculate proportion of frames in individual behaviors
logger('Calculating frame-level individuality statistics...', 'INFO');

% Initialize frame-level statistics containers
animal_condition_individual_frames = containers.Map();
animal_condition_total_frames = containers.Map();

% First, count total frames per animal+condition combination (all frames, individual or not)
for i = 1:length(global_animal_condition_ids)
    animal_condition_id = global_animal_condition_ids{i};
    
    if isKey(animal_condition_total_frames, animal_condition_id)
        animal_condition_total_frames(animal_condition_id) = animal_condition_total_frames(animal_condition_id) + 1;
    else
        animal_condition_total_frames(animal_condition_id) = 1;
    end
end

% Count frames in individual clusters per animal+condition combination
for c_idx = 1:length(unique_clusters)
    if cluster_is_individual(c_idx) && ~strcmp(cluster_dominant_animal{c_idx}, 'None')
        cluster_id = unique_clusters(c_idx);
        dominant_animal_condition = cluster_dominant_animal{c_idx};
        
        % Count frames belonging to the dominant animal+condition in this individual cluster
        cluster_frames = cluster_assignments == cluster_id;
        frames_for_dominant = 0;
        
        for i = 1:length(global_animal_condition_ids)
            if cluster_frames(i) && strcmp(global_animal_condition_ids{i}, dominant_animal_condition)
                frames_for_dominant = frames_for_dominant + 1;
            end
        end
        
        % Add to individual frames count for this animal+condition
        if isKey(animal_condition_individual_frames, dominant_animal_condition)
            animal_condition_individual_frames(dominant_animal_condition) = ...
                animal_condition_individual_frames(dominant_animal_condition) + frames_for_dominant;
        else
            animal_condition_individual_frames(dominant_animal_condition) = frames_for_dominant;
        end
    end
end

% Store frame-level statistics
individuality_analysis.animal_individual_frames = animal_condition_individual_frames;
individuality_analysis.animal_total_frames = animal_condition_total_frames;

% Log frame-level summary
total_global_frames = length(global_animal_condition_ids);
total_individual_frames = 0;
animal_condition_keys = keys(animal_condition_individual_frames);
for i = 1:length(animal_condition_keys)
    total_individual_frames = total_individual_frames + animal_condition_individual_frames(animal_condition_keys{i});
end

logger(sprintf('Frame-level individuality: %d/%d frames (%.1f%%) are in individual behaviors globally', ...
    total_individual_frames, total_global_frames, (total_individual_frames/total_global_frames)*100), 'INFO');

%% Create individuality visualization functions

function plot_individuality_tsne_map_condition(analysisstruct, individuality_data, condition_idx, visualize, condition_name)
    % Create figure showing individual clusters on t-SNE map using GLOBAL individuality analysis
    fig = figure('Name', 'Individuality Analysis: t-SNE Maps', 'Visible', visualize);
    set(fig, 'Position', [300, 100, 1400, 800]);
    set(fig, 'Color', 'w');
    
    % Get GLOBAL t-SNE data for reference (all conditions, all animals)
    global_zvals = analysisstruct.zValues;
    global_clusters = analysisstruct.annot_reordered{end,end};
    
    % Get condition-specific data for highlighting
    condition_zvals = global_zvals(condition_idx, :);
    condition_clusters = global_clusters(condition_idx);
    
    % Subplot 1: Individual vs non-individual clusters
    subplot(1, 2, 1);
    % Plot ALL global data in light gray as background
    plot(global_zvals(:,1), global_zvals(:,2), '.', ...
        'Color', [0.95, 0.95, 0.95], 'MarkerSize', 1);
    hold on;
    
    % Plot selected condition data in gray
    plot(condition_zvals(:,1), condition_zvals(:,2), '.', ...
        'Color', [0.7, 0.7, 0.7], 'MarkerSize', 2);
    
    % Highlight individual clusters from selected condition (based on GLOBAL individuality)
    for c_idx = 1:length(individuality_data.cluster_ids)
        if individuality_data.cluster_is_individual(c_idx)
            cluster_id = individuality_data.cluster_ids(c_idx);
            % Find frames in the selected condition that belong to this globally individual cluster
            cluster_frames = condition_clusters == cluster_id;
            if sum(cluster_frames) > 0
                plot(condition_zvals(cluster_frames,1), ...
                     condition_zvals(cluster_frames,2), '.', ...
                     'Color', [1, 0, 0], 'MarkerSize', 4);
            end
        end
    end
    
    title(sprintf('Globally Individual Clusters (>%d%% dominance) in %s', individuality_data.threshold*100, condition_name));
    xlabel('t-SNE 1 (Global Space)');
    ylabel('t-SNE 2 (Global Space)');
    axis equal;
    legend({'All data (global)', sprintf('%s condition', condition_name), 'Globally individual clusters'}, 'Location', 'best');
    
    % Subplot 2: Dominance percentage heatmap (global dominance values)
    subplot(1, 2, 2);
    
    % Plot global background
    plot(global_zvals(:,1), global_zvals(:,2), '.', ...
        'Color', [0.95, 0.95, 0.95], 'MarkerSize', 1);
    hold on;
    
    % Create a colormap based on GLOBAL dominance percentage for condition data
    dominance_values = zeros(size(condition_zvals, 1), 1);
    
    for c_idx = 1:length(individuality_data.cluster_ids)
        cluster_id = individuality_data.cluster_ids(c_idx);
        cluster_frames = condition_clusters == cluster_id;
        dominance_values(cluster_frames) = individuality_data.cluster_dominance_percentage(c_idx);
    end
    
    % Create scatter plot colored by global dominance (only for condition data)
    scatter(condition_zvals(:,1), condition_zvals(:,2), 12, dominance_values, 'filled');
    
    % Set colormap and colorbar
    colormap(hot);
    colorbar;
    caxis([0, 1]);
    
    title(sprintf('Global Dominance Percentage - %s', condition_name));
    xlabel('t-SNE 1 (Global Space)');
    ylabel('t-SNE 2 (Global Space)');
    axis equal;
    
    % Add main title
    sgtitle(sprintf('Pose Individuality Analysis - %s (Global Analysis: %d animals, %d clusters)', ...
        condition_name, length(individuality_data.global_unique_animals), length(individuality_data.cluster_ids)), ...
        'FontSize', 16, 'FontWeight', 'bold');
end

function plot_individuality_per_animal_condition(analysisstruct, individuality_data, condition_idx, ~, visualize, condition_name, selected_condition)
    % Create figure showing individual clusters for each animal using GLOBAL t-SNE coordinates
    unique_animals = individuality_data.unique_animals;
    n_animals = length(unique_animals);
    
    % Get GLOBAL t-SNE data for reference
    global_zvals = analysisstruct.zValues;
    global_clusters = analysisstruct.annot_reordered{end,end};
    
    % Get condition-specific data
    condition_zvals = global_zvals(condition_idx, :);
    condition_clusters = global_clusters(condition_idx);
    
    % Create subplot grid
    n_cols = ceil(sqrt(n_animals));
    n_rows = ceil(n_animals / n_cols);
    
    fig = figure('Name', 'Individual Clusters per Animal', 'Visible', visualize);
    set(fig, 'Position', [400, 100, 1600, 1200]);
    set(fig, 'Color', 'w');
    
    for a_idx = 1:n_animals
        animal_id = unique_animals{a_idx};
        
        subplot(n_rows, n_cols, a_idx);
        
        % Plot global background in very light gray
        plot(global_zvals(:,1), global_zvals(:,2), '.', ...
            'Color', [0.95, 0.95, 0.95], 'MarkerSize', 0.5);
        hold on;
        
        % Plot all condition points in gray
        plot(condition_zvals(:,1), condition_zvals(:,2), '.', ...
            'Color', [0.7, 0.7, 0.7], 'MarkerSize', 1);
        
        % Highlight clusters dominated by this animal+condition combination
        % Create the animal+condition identifier for this animal in the selected condition
        animal_condition_id = [animal_id '_' selected_condition];
        
        for c_idx = 1:length(individuality_data.cluster_ids)
            if strcmp(individuality_data.cluster_dominant_animal{c_idx}, animal_condition_id) && ...
               individuality_data.cluster_is_individual(c_idx)
                cluster_id = individuality_data.cluster_ids(c_idx);
                cluster_frames = condition_clusters == cluster_id;
                if sum(cluster_frames) > 0
                    plot(condition_zvals(cluster_frames,1), ...
                         condition_zvals(cluster_frames,2), '.', ...
                         'Color', [1, 0, 0], 'MarkerSize', 3);
                end
            end
        end
        
        % Calculate individual clusters for this animal+condition combination
        individual_count = 0;
        if isKey(individuality_data.animal_individual_clusters, animal_condition_id)
            individual_count = individuality_data.animal_individual_clusters(animal_condition_id);
        end
        
        title(sprintf('%s (%d individual)', animal_id, individual_count));
        xlabel('t-SNE 1 (Global)');
        ylabel('t-SNE 2 (Global)');
        axis equal;
        axis tight;
    end
    
    sgtitle(sprintf('Individual Clusters per Animal - %s (Global Individuality Analysis)', condition_name), 'FontSize', 16, 'FontWeight', 'bold');
end

function plot_individuality_statistics_condition(individuality_data, visualize, export_folder, do_export, condition_name, selected_condition, global_unique_animal_condition_ids)
    % Create comprehensive statistical plots for condition individuality
    
    fig = figure('Name', 'Individuality Statistics', 'Visible', visualize);
    set(fig, 'Position', [500, 100, 1200, 800]);
    set(fig, 'Color', 'w');
    
    % Subplot 1: Pie chart of individual vs non-individual clusters
    subplot(2, 2, 1);
    individual_count = sum(individuality_data.cluster_is_individual);
    total_count = length(individuality_data.cluster_ids);
    non_individual_count = total_count - individual_count;
    
    if individual_count > 0 || non_individual_count > 0
        pie_data = [individual_count, non_individual_count];
        pie_labels = {sprintf('Individual (%.1f%%)', (individual_count/total_count)*100), ...
                      sprintf('Shared (%.1f%%)', (non_individual_count/total_count)*100)};
        pie_handle = pie(pie_data, pie_labels);
        
        % Color the pie slices
        set(pie_handle(1), 'FaceColor', [1, 0.2, 0.2]); % Red for individual
        if length(pie_handle) >= 3
            set(pie_handle(3), 'FaceColor', [0.7, 0.7, 0.7]); % Gray for shared
        end
    end
    
    title(sprintf('Cluster Individuality (Threshold: %d%%)', individuality_data.threshold*100));
    
    % Subplot 2: Histogram of dominance percentages
    subplot(2, 2, 2);
    valid_dominance = individuality_data.cluster_dominance_percentage(individuality_data.cluster_dominance_percentage > 0);
    if ~isempty(valid_dominance)
        histogram(valid_dominance, 15, 'FaceColor', [0.3, 0.6, 1], 'EdgeColor', 'black');
        xlabel('Dominance Percentage');
        ylabel('Number of Clusters');
        title('Distribution of Cluster Dominance');
        xlim([0, 1]);
        
        % Add threshold line
        hold on;
        line([individuality_data.threshold, individuality_data.threshold], ylim, ...
            'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
        legend('Clusters', sprintf('%d%% Threshold', individuality_data.threshold*100), 'Location', 'best');
    end
    
    % Subplot 3: Per-animal+condition individual cluster counts for selected condition
    subplot(2, 2, [3, 4]);
    
    % Get animal+condition combinations for the selected condition only
    % Filter global animal+condition combinations to only those in the selected condition
    condition_animal_condition_ids = {};
    for i = 1:length(global_unique_animal_condition_ids)
        animal_condition_id = global_unique_animal_condition_ids{i};
        if endsWith(animal_condition_id, ['_' selected_condition]) % Check if ends with selected condition
            condition_animal_condition_ids{end+1} = animal_condition_id;
        end
    end
    
    % Prepare data for bar chart using animal+condition combinations
    animal_condition_combinations = condition_animal_condition_ids;
    individual_counts = zeros(length(animal_condition_combinations), 1);
    total_present_counts = zeros(length(animal_condition_combinations), 1);
    
    for i = 1:length(animal_condition_combinations)
        animal_condition_id = animal_condition_combinations{i};
        if isKey(individuality_data.animal_individual_clusters, animal_condition_id)
            individual_counts(i) = individuality_data.animal_individual_clusters(animal_condition_id);
        end
        if isKey(individuality_data.animal_total_clusters_present, animal_condition_id)
            total_present_counts(i) = individuality_data.animal_total_clusters_present(animal_condition_id);
        end
    end
    
    % Calculate percentages for normalized stacked bar chart
    individual_percentages = zeros(length(animal_condition_combinations), 1);
    shared_percentages = zeros(length(animal_condition_combinations), 1);
    
    for i = 1:length(animal_condition_combinations)
        if total_present_counts(i) > 0
            individual_percentages(i) = (individual_counts(i) / total_present_counts(i)) * 100;
            shared_percentages(i) = ((total_present_counts(i) - individual_counts(i)) / total_present_counts(i)) * 100;
        end
    end
    
    % Create normalized stacked bar chart (percentages)
    bar_data = [individual_percentages, shared_percentages];
    bar_handle = bar(bar_data, 'stacked');
    set(bar_handle(1), 'FaceColor', [1, 0.2, 0.2]); % Red for individual
    set(bar_handle(2), 'FaceColor', [0.7, 0.7, 0.7]); % Gray for shared
    
    xlabel('Animal+Condition ID');
    ylabel('Percentage of Clusters');
    title('Individual vs Shared Clusters per Animal+Condition (Normalized to 100%)');
    legend('Individual Clusters (%)', 'Shared Clusters (%)', 'Location', 'best');
    
    % Set y-axis to 0-100%
    ylim([0, 100]);
    
    % Add percentage labels on bars for individual clusters
    for i = 1:length(animal_condition_combinations)
        if individual_percentages(i) > 5 % Only show label if segment is large enough
            text(i, individual_percentages(i)/2, sprintf('%.1f%%', individual_percentages(i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'Color', 'white', 'FontWeight', 'bold', 'FontSize', 8);
        end
        % Add total cluster count as text above each bar
        text(i, 102, sprintf('n=%d', total_present_counts(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 8, 'Color', 'black');
    end
    
    % Set x-axis labels - extract just animal name for cleaner display
    display_labels = cellfun(@(x) x(1:find(x=='_',1,'last')-1), animal_condition_combinations, 'UniformOutput', false);
    set(gca, 'XTickLabel', display_labels);
    set(gca, 'TickDir', 'out');
    xtickangle(45);
    
    sgtitle(sprintf('Pose Individuality Statistical Analysis - %s (Global Analysis)', condition_name), 'FontSize', 16, 'FontWeight', 'bold');
    
    % Export if requested
    if do_export
        stats_filename = fullfile(export_folder, sprintf('%s_individuality_statistics.png', lower(condition_name)));
        saveas(fig, stats_filename);
        logger(sprintf('Saved %s individuality statistics: %s', condition_name, stats_filename), 'INFO');
    end
end

function plot_frame_individuality_analysis(individuality_data, visualize, export_folder, do_export, condition_name, selected_condition, global_unique_animal_condition_ids)
    % Create comprehensive frame-level individuality analysis plots
    
    fig = figure('Name', 'Frame-Level Individuality Analysis', 'Visible', visualize);
    set(fig, 'Position', [600, 100, 1400, 600]);
    set(fig, 'Color', 'w');
    
    % Get animal+condition combinations for the selected condition only
    condition_animal_condition_ids = {};
    for i = 1:length(global_unique_animal_condition_ids)
        animal_condition_id = global_unique_animal_condition_ids{i};
        if endsWith(animal_condition_id, ['_' selected_condition])
            condition_animal_condition_ids{end+1} = animal_condition_id;
        end
    end
    
    % Prepare data for frame-level analysis
    animal_condition_combinations = condition_animal_condition_ids;
    individual_frame_counts = zeros(length(animal_condition_combinations), 1);
    total_frame_counts = zeros(length(animal_condition_combinations), 1);
    
    for i = 1:length(animal_condition_combinations)
        animal_condition_id = animal_condition_combinations{i};
        
        % Get individual frame count
        if isKey(individuality_data.animal_individual_frames, animal_condition_id)
            individual_frame_counts(i) = individuality_data.animal_individual_frames(animal_condition_id);
        end
        
        % Get total frame count
        if isKey(individuality_data.animal_total_frames, animal_condition_id)
            total_frame_counts(i) = individuality_data.animal_total_frames(animal_condition_id);
        end
    end
    
    % Calculate percentages
    individual_frame_percentages = zeros(length(animal_condition_combinations), 1);
    shared_frame_percentages = zeros(length(animal_condition_combinations), 1);
    
    for i = 1:length(animal_condition_combinations)
        if total_frame_counts(i) > 0
            individual_frame_percentages(i) = (individual_frame_counts(i) / total_frame_counts(i)) * 100;
            shared_frame_percentages(i) = ((total_frame_counts(i) - individual_frame_counts(i)) / total_frame_counts(i)) * 100;
        end
    end
    
    % Subplot 1: Stacked bar chart showing frame percentages
    subplot(1, 2, 1);
    bar_data = [individual_frame_percentages, shared_frame_percentages];
    bar_handle = bar(bar_data, 'stacked');
    set(bar_handle(1), 'FaceColor', [1, 0.2, 0.2]); % Red for individual frames
    set(bar_handle(2), 'FaceColor', [0.7, 0.7, 0.7]); % Gray for shared frames
    
    xlabel('Animal ID');
    ylabel('Percentage of Frames');
    title('Individual vs Shared Behavioral Time per Animal+Condition');
    legend('Individual Behavior Time (%)', 'Shared Behavior Time (%)', 'Location', 'best');
    
    % Set y-axis to 0-100%
    ylim([0, 100]);
    
    % Add percentage labels on bars for individual frames
    for i = 1:length(animal_condition_combinations)
        if individual_frame_percentages(i) > 5
            text(i, individual_frame_percentages(i)/2, sprintf('%.1f%%', individual_frame_percentages(i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'Color', 'white', 'FontWeight', 'bold', 'FontSize', 8);
        end
        % Add total frame count as text above each bar
        text(i, 102, sprintf('n=%d', total_frame_counts(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 8, 'Color', 'black');
    end
    
    % Set x-axis labels
    display_labels = cellfun(@(x) x(1:find(x=='_',1,'last')-1), animal_condition_combinations, 'UniformOutput', false);
    set(gca, 'XTickLabel', display_labels);
    set(gca, 'TickDir', 'out');
    xtickangle(45);
    
    % Subplot 2: Comparison of cluster vs frame individuality
    subplot(1, 2, 2);
    
    % Get cluster percentages for comparison
    cluster_individual_counts = zeros(length(animal_condition_combinations), 1);
    cluster_total_counts = zeros(length(animal_condition_combinations), 1);
    
    for i = 1:length(animal_condition_combinations)
        animal_condition_id = animal_condition_combinations{i};
        
        if isKey(individuality_data.animal_individual_clusters, animal_condition_id)
            cluster_individual_counts(i) = individuality_data.animal_individual_clusters(animal_condition_id);
        end
        if isKey(individuality_data.animal_total_clusters_present, animal_condition_id)
            cluster_total_counts(i) = individuality_data.animal_total_clusters_present(animal_condition_id);
        end
    end
    
    cluster_individual_percentages = zeros(length(animal_condition_combinations), 1);
    for i = 1:length(animal_condition_combinations)
        if cluster_total_counts(i) > 0
            cluster_individual_percentages(i) = (cluster_individual_counts(i) / cluster_total_counts(i)) * 100;
        end
    end
    
    % Plot comparison
    x_positions = 1:length(animal_condition_combinations);
    bar_width = 0.35;
    
    bar(x_positions - bar_width/2, cluster_individual_percentages, bar_width, ...
        'FaceColor', [0.2, 0.6, 1], 'DisplayName', 'Cluster Individuality (%)');
    hold on;
    bar(x_positions + bar_width/2, individual_frame_percentages, bar_width, ...
        'FaceColor', [1, 0.2, 0.2], 'DisplayName', 'Frame Individuality (%)');
    
    xlabel('Animal ID');
    ylabel('Individuality Percentage');
    title('Cluster vs Frame Individuality Comparison');
    legend('Location', 'best');
    ylim([0, 100]);
    
    % Set x-axis labels
    set(gca, 'XTick', x_positions);
    set(gca, 'XTickLabel', display_labels);
    set(gca, 'TickDir', 'out');
    xtickangle(45);
    
    % Add main title
    sgtitle(sprintf('Frame-Level Behavioral Individuality Analysis - %s', condition_name), ...
        'FontSize', 16, 'FontWeight', 'bold');
    
    % Export if requested
    if do_export
        frame_stats_filename = fullfile(export_folder, sprintf('%s_frame_individuality_analysis.png', lower(condition_name)));
        saveas(fig, frame_stats_filename);
        logger(sprintf('Saved %s frame individuality analysis: %s', condition_name, frame_stats_filename), 'INFO');
    end
end

%% Generate individuality visualizations
logger(sprintf('Creating individuality visualizations for %s data', selected_condition_name), 'INFO');

% Create t-SNE individuality maps
plot_individuality_tsne_map_condition(analysisstruct, individuality_analysis, idx_condition, visualize, selected_condition_name);

% Create per-animal individuality maps
plot_individuality_per_animal_condition(analysisstruct, individuality_analysis, idx_condition, condition_animal_names, visualize, selected_condition_name, selected_condition);

% Create statistical analysis plots (cluster-level)
plot_individuality_statistics_condition(individuality_analysis, visualize, export_folder, do_export, selected_condition_name, selected_condition, global_unique_animal_condition_ids);

% Create frame-level individuality analysis plots
plot_frame_individuality_analysis(individuality_analysis, visualize, export_folder, do_export, selected_condition_name, selected_condition, global_unique_animal_condition_ids);

%% Create individuality summary table
logger(sprintf('Creating individuality summary table for %s condition (using global individuality metrics)', selected_condition_name), 'INFO');

% Create detailed summary table focusing on animal+condition combinations in selected condition
% but using their GLOBAL individuality metrics

% Get animal+condition combinations for the selected condition only
condition_animal_condition_ids = {};
for i = 1:length(global_unique_animal_condition_ids)
    animal_condition_id = global_unique_animal_condition_ids{i};
    if endsWith(animal_condition_id, ['_' selected_condition]) % Check if ends with selected condition
        condition_animal_condition_ids{end+1} = animal_condition_id;
    end
end

animal_condition_ids = condition_animal_condition_ids;

% Debug information
logger(sprintf('Found %d global animal+condition combinations total', length(global_unique_animal_condition_ids)), 'INFO');
logger(sprintf('Creating summary table for %d animal+condition combinations in %s condition', ...
    length(animal_condition_ids), selected_condition_name), 'INFO');

% Enhanced debugging
if isempty(animal_condition_ids)
    logger('WARNING: No animal+condition combinations found for selected condition!', 'WARNING');
    logger(sprintf('Available global combinations: %s', strjoin(global_unique_animal_condition_ids, ', ')), 'INFO');
    logger(sprintf('Looking for combinations ending with: _%s', selected_condition), 'INFO');
    
    % Try alternative matching approach - check each global combination
    logger('Checking each global combination for matches:', 'INFO');
    for i = 1:length(global_unique_animal_condition_ids)
        combo = global_unique_animal_condition_ids{i};
        logger(sprintf('  %s - ends with _%s? %s', combo, selected_condition, ...
            string(endsWith(combo, ['_' selected_condition]))), 'INFO');
    end
    return;
end

% Ensure we have valid data before creating the table
logger(sprintf('Validated %d animal+condition combinations for table creation', length(animal_condition_ids)), 'INFO');

individual_cluster_counts = zeros(length(animal_condition_ids), 1);
total_present_cluster_counts = zeros(length(animal_condition_ids), 1);
total_dominant_cluster_counts = zeros(length(animal_condition_ids), 1);
individuality_percentages = zeros(length(animal_condition_ids), 1);
individual_frame_counts = zeros(length(animal_condition_ids), 1);
total_frame_counts = zeros(length(animal_condition_ids), 1);
frame_individuality_percentages = zeros(length(animal_condition_ids), 1);

for i = 1:length(animal_condition_ids)
    animal_condition_id = animal_condition_ids{i};
    
    % Get GLOBAL individual cluster count for this animal+condition combination
    if isKey(individuality_analysis.animal_individual_clusters, animal_condition_id)
        individual_cluster_counts(i) = individuality_analysis.animal_individual_clusters(animal_condition_id);
    end
    
    % Get GLOBAL total clusters present count for this animal+condition combination
    if isKey(individuality_analysis.animal_total_clusters_present, animal_condition_id)
        total_present_cluster_counts(i) = individuality_analysis.animal_total_clusters_present(animal_condition_id);
    end
    
    % Get GLOBAL total dominant cluster count for this animal+condition combination
    if isKey(individuality_analysis.animal_total_dominant_clusters, animal_condition_id)
        total_dominant_cluster_counts(i) = individuality_analysis.animal_total_dominant_clusters(animal_condition_id);
    end
    
    % Get frame-level statistics for this animal+condition combination
    if isKey(individuality_analysis.animal_individual_frames, animal_condition_id)
        individual_frame_counts(i) = individuality_analysis.animal_individual_frames(animal_condition_id);
    end
    
    if isKey(individuality_analysis.animal_total_frames, animal_condition_id)
        total_frame_counts(i) = individuality_analysis.animal_total_frames(animal_condition_id);
    end
    
    % Calculate GLOBAL individuality percentages
    if total_present_cluster_counts(i) > 0
        individuality_percentages(i) = (individual_cluster_counts(i) / total_present_cluster_counts(i)) * 100;
    end
    
    if total_frame_counts(i) > 0
        frame_individuality_percentages(i) = (individual_frame_counts(i) / total_frame_counts(i)) * 100;
    end
end

% Create table with proper error checking and consistent array lengths
logger(sprintf('Creating table with %d rows of data', length(animal_condition_ids)), 'INFO');

% Ensure all arrays have the same length
n_combinations = length(animal_condition_ids);
individual_cluster_counts = individual_cluster_counts(1:n_combinations);
total_present_cluster_counts = total_present_cluster_counts(1:n_combinations);
total_dominant_cluster_counts = total_dominant_cluster_counts(1:n_combinations);
individuality_percentages = individuality_percentages(1:n_combinations);
individual_frame_counts = individual_frame_counts(1:n_combinations);
total_frame_counts = total_frame_counts(1:n_combinations);
frame_individuality_percentages = frame_individuality_percentages(1:n_combinations);

% Verify all arrays have the same length
logger(sprintf('Array lengths - IDs: %d, Individual Clusters: %d, Present: %d, Dominant: %d, Cluster %%: %d, Individual Frames: %d, Total Frames: %d, Frame %%: %d', ...
    length(animal_condition_ids), length(individual_cluster_counts), ...
    length(total_present_cluster_counts), length(total_dominant_cluster_counts), ...
    length(individuality_percentages), length(individual_frame_counts), ...
    length(total_frame_counts), length(frame_individuality_percentages)), 'INFO');

% Convert cell array to categorical for better table handling
animal_condition_ids_cat = categorical(animal_condition_ids);

% Create table properly with all variables at once (including frame-level statistics)
individuality_table = table(animal_condition_ids_cat', individual_cluster_counts, total_present_cluster_counts, ...
    total_dominant_cluster_counts, individuality_percentages, individual_frame_counts, total_frame_counts, frame_individuality_percentages, ...
    'VariableNames', {'Animal_Condition_ID', 'Individual_Clusters_Global', ...
    'Total_Clusters_Present_Global', 'Total_Dominant_Clusters_Global', 'Individuality_Percentage_Global', ...
    'Individual_Frames_Global', 'Total_Frames_Global', 'Frame_Individuality_Percentage_Global'});

% Add overall statistics
total_clusters_analyzed = length(individuality_analysis.cluster_ids);
total_individual_clusters = sum(individuality_analysis.cluster_is_individual);
overall_individuality_percentage = (total_individual_clusters / total_clusters_analyzed) * 100;

% Save individuality table if export is enabled
if do_export
    individuality_filename = fullfile(export_folder, sprintf('%s_pose_individuality_analysis_global.csv', lower(selected_condition_name)));
    writetable(individuality_table, individuality_filename);
    logger(sprintf('Saved %s individuality analysis table (global metrics): %s', selected_condition_name, individuality_filename), 'INFO');
end

%% Display summary statistics
logger(sprintf('=== %s CONDITION SUMMARY ===', upper(selected_condition_name)), 'INFO');
logger(sprintf('Total animals in %s condition: %d', selected_condition_name, num_animals), 'INFO');
logger(sprintf('Total frames in %s condition: %d', selected_condition_name, sum(idx_condition)), 'INFO');
logger(sprintf('NOTE: Individuality metrics calculated using GLOBAL dataset (%d animals, %d total frames)', ...
    length(global_unique_animals), length(cluster_assignments)), 'INFO');

for i = 1:num_animals
    animal = unique_animals{i};
    idx_animal = strcmp(condition_animal_names, animal);
    frames_count = sum(idx_animal);
    percentage = (frames_count / sum(idx_condition)) * 100;
    logger(sprintf('Animal %s: %d frames (%.1f%%)', animal, frames_count, percentage), 'INFO');
end

% Display individuality results
logger('=== POSE INDIVIDUALITY ANALYSIS RESULTS (GLOBAL) ===', 'INFO');
logger(sprintf('Threshold: %.0f%% dominance', individuality_analysis.threshold * 100), 'INFO');
logger(sprintf('Total clusters analyzed (global): %d', total_clusters_analyzed), 'INFO');
logger(sprintf('Individual clusters (global): %d (%.1f%%)', total_individual_clusters, overall_individuality_percentage), 'INFO');
logger(sprintf('Shared clusters (global): %d (%.1f%%)', total_clusters_analyzed - total_individual_clusters, ...
    100 - overall_individuality_percentage), 'INFO');
logger(' ', 'INFO');
logger(sprintf('Per-animal+condition individuality for %s condition (based on global analysis):', selected_condition_name), 'INFO');
logger('NOTE: Each row represents a unique animal+condition combination (e.g., Animal1_S, Animal2_S)', 'INFO');
logger('This ensures same animals in different conditions are analyzed separately', 'INFO');
logger(' ', 'INFO');
logger('TABLE COLUMNS EXPLANATION:', 'INFO');
logger('- Individual_Clusters_Global: Number of clusters uniquely dominated by this animal+condition', 'INFO');
logger('- Individuality_Percentage_Global: Percentage of clusters that are individual (cluster-level metric)', 'INFO');
logger('- Individual_Frames_Global: Number of frames spent in individual behaviors', 'INFO');
logger('- Frame_Individuality_Percentage_Global: Percentage of time spent in individual behaviors (frame-level metric)', 'INFO');
logger(' ', 'INFO');
disp(individuality_table);

logger(sprintf('%s analysis complete - individuality assessed using global dataset for robust metrics', selected_condition_name), 'INFO');
logger(' ', 'INFO');
logger('=== CRITICAL FIX IMPLEMENTED ===', 'INFO');
logger('✓ Animal+condition combinations now treated as separate entities', 'INFO');
logger('✓ Same animal in different conditions (e.g., Animal1_S vs Animal1_F) analyzed independently', 'INFO');
logger('✓ Bar charts and statistics show animal+condition specific individuality', 'INFO');
logger('✓ Global analysis prevents bias from sparse conditions', 'INFO');
logger('✓ Frame-level analysis shows proportion of behavioral time in individual behaviors', 'INFO');
logger('✓ Cluster vs Frame individuality comparison reveals behavioral patterns', 'INFO');
logger('Analysis complete!', 'INFO');
