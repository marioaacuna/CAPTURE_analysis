%% ========================================================================
%                   CROSS-CONDITION POSE INDIVIDUALITY ANALYSIS
% ========================================================================
%
% DESCRIPTION:
% This script performs comprehensive cross-condition individuality analysis 
% across ALL experimental conditions simultaneously. It calculates cluster-level 
% and frame-level individuality metrics for each animal+condition combination, 
% performs rigorous statistical comparisons across conditions, and generates 
% publication-quality visualizations to understand how behavioral individuality 
% varies across experimental paradigms.
%
% MAIN OBJECTIVES:
% 1. Calculate individuality metrics across ALL conditions (B, S, F, H, N) simultaneously
% 2. Store cluster-level and frame-level individuality data per animal+condition combination
% 3. Perform comprehensive statistical analysis comparing individuality across conditions
% 4. Generate publication-quality visualizations and statistical comparison plots
% 5. Export comprehensive datasets for external statistical analysis and reporting
%
% KEY CONCEPTS:
% - Individual Clusters: Behavioral clusters where ≥80% of frames belong to one animal+condition
% - Shared Clusters: Behavioral clusters with mixed animal+condition representation
% - Cross-Condition Analysis: Systematic comparison of individuality patterns across experimental conditions
% - Global Reference Frame: Uses complete dataset (all animals, all conditions) for unbiased individuality assessment
% - Animal+Condition Identity: Treats each animal in each condition as unique entity (e.g., "Animal1_B" vs "Animal1_F")
%
% KEY METRICS CALCULATED:
% Cluster-Level Metrics (per animal+condition combination):
% - Individual_Clusters_Global: Number of clusters with ≥80% dominance globally
% - Total_Clusters_Present_Global: Total clusters where animal+condition appears globally
% - Total_Dominant_Clusters_Global: Clusters where animal+condition is most frequent globally
% - Cluster_Individuality_Percentage_Global: (Individual_Clusters / Total_Clusters_Present) × 100
%
% Frame-Level Metrics (per animal+condition combination):
% - Individual_Frames_Global: Number of frames spent in globally individual behaviors
% - Total_Frames_Global: Total frames for this animal+condition combination globally
% - Frame_Individuality_Percentage_Global: (Individual_Frames / Total_Frames) × 100
%
% STATISTICAL ANALYSES:
% - One-way ANOVA comparing cluster individuality percentages across conditions
% - One-way ANOVA comparing frame individuality percentages across conditions
% - Post-hoc multiple comparison tests for pairwise condition differences
% - Effect size calculations and confidence intervals
% - Correlation analysis between cluster and frame individuality measures
% - Individual animal trajectory analysis across conditions (where applicable)
%
% BIOLOGICAL INTERPRETATION:
% High individuality indicates an animal+condition combination exhibits many unique 
% behavioral patterns not commonly shared with other animal+condition combinations.
% Cross-condition analysis reveals:
% - Whether experimental manipulations increase or decrease behavioral individuality
% - Condition-specific effects on behavioral diversity and stereotypy
% - Individual animal responses to different experimental paradigms
% - Population-level vs. individual-level behavioral adaptations
% - Baseline individual differences vs. treatment-induced changes
%
% Cluster vs Frame Individuality Interpretation:
% - Cluster Individuality: Measures diversity of unique behavioral types per animal+condition
% - Frame Individuality: Measures time allocation to unique vs. shared behaviors
% - Cross-condition patterns reveal how treatments affect behavioral repertoire diversity
% - Helps distinguish between changes in behavioral variety vs. behavioral frequency
%
% ANALYSIS WORKFLOW:
% 1. Global Data Preparation
%    - Load complete dataset with all animals and conditions
%    - Extract animal+condition identifiers and cluster assignments
%    - Prepare global reference frame for unbiased individuality calculation
%
% 2. Global Individuality Analysis
%    - Calculate cluster composition using entire dataset (all conditions)
%    - Identify individual vs. shared clusters based on global dominance patterns
%    - Count cluster-level and frame-level statistics per animal+condition combination
%    - CRITICAL: Uses global reference to prevent condition-specific bias
%
% 3. Condition-Wise Data Organization
%    - Group animal+condition combinations by experimental condition
%    - Calculate summary statistics (means, standard errors) per condition
%    - Prepare data structures for statistical analysis and visualization
%
% 4. Statistical Analysis
%    - Perform one-way ANOVA comparing individuality across conditions
%    - Conduct post-hoc tests for pairwise condition comparisons
%    - Calculate effect sizes and statistical significance measures
%    - Generate comprehensive statistical summary reports
%
% 5. Visualization Generation
%    - Bar graphs with error bars showing mean individuality per condition
%    - Box plots showing distribution of individuality within conditions
%    - Correlation plots between cluster and frame individuality measures
%    - Individual animal trajectory plots across conditions (when applicable)
%    - Statistical significance annotations and effect size indicators
%
% 6. Data Export and Reporting
%    - Export comprehensive summary tables with all metrics
%    - Save statistical analysis results and effect sizes
%    - Generate publication-ready figures in multiple formats
%    - Create detailed analysis logs and methodology documentation
%
% VISUALIZATIONS GENERATED:
% Statistical Comparison Plots:
% - Mean individuality bar charts with error bars and significance annotations
% - Box plots showing individuality distributions within each condition
% - Correlation scatter plots between cluster and frame individuality measures
% - Individual animal trajectory plots across multiple conditions
%
% Publication-Quality Outputs:
% - High-resolution figures suitable for manuscripts (PNG, PDF formats)
% - Statistical summary plots with effect sizes and confidence intervals
% - Cross-condition comparison matrices and heatmaps
% - Comprehensive figure legends and statistical annotations
%
% DATA OUTPUTS:
% Comprehensive Summary Tables:
% - cross_condition_individuality_summary.csv: Complete per-animal+condition metrics
% - statistical_analysis_results.mat: ANOVA results, effect sizes, post-hoc tests
% - condition_means_and_statistics.csv: Per-condition summary statistics
%
% Analysis Documentation:
% - Detailed console logging with analysis progress and methodology
% - Statistical test results with interpretation guidelines
% - Complete individuality_analysis structure saved to workspace
%
% IMPORTANT METHODOLOGICAL NOTES:
% Global Reference Frame: All individuality calculations use the complete dataset 
% (all animals, all conditions) as reference. This ensures unbiased comparison 
% across conditions and prevents inflation of individuality scores in conditions 
% with fewer animals or limited behavioral repertoires.
%
% Animal+Condition Identity: Each animal in each condition is treated as a unique 
% entity (e.g., "Animal1_B", "Animal1_F", "Animal1_S"). This approach captures 
% how the same animal's behavioral individuality changes across different 
% experimental conditions, which is critical for understanding treatment effects.
%
% Statistical Robustness: Multiple statistical approaches (ANOVA, post-hoc tests, 
% effect sizes) provide comprehensive assessment of condition differences while 
% controlling for multiple comparisons and ensuring reproducible results.
%
% CONFIGURATION OPTIONS:
% - individuality_threshold: Dominance threshold for individual classification (default: 0.8)
% - visualize: Control figure display ('on' for debug mode, 'off' for production)
% - do_export: Enable/disable comprehensive figure and data export
% - export_folder: Destination for all generated outputs and analysis results
%
% USAGE SCENARIOS:
% Research Applications:
% - Comparing behavioral individuality across different experimental treatments
% - Assessing how pharmacological interventions affect behavioral diversity
% - Studying individual differences in response to environmental manipulations
% - Quantifying treatment effects on behavioral stereotypy vs. flexibility
%
% Quality Control:
% - Validating experimental design effects on behavioral measurements
% - Identifying potential batch effects or technical artifacts in data
% - Assessing completeness and balance of experimental conditions
%
% DEPENDENCIES:
% - CAPTURE analysis framework with complete individuality analysis pipeline
% - Statistics and Machine Learning Toolbox (for ANOVA and post-hoc tests)
% - Complete behavioral clustering and t-SNE embedding results
% - Properly formatted animal+condition identifier arrays
%
% USAGE EXAMPLE:
% % Ensure all required data structures are available
% GC = general_configs;
% % Run complete cross-condition analysis (user will be prompted for analysis mode)
% run('script__XX_analysis_individuality.m');
%
% AUTHORS: CAPTURE Analysis Team
% CREATED: December 2025
% LAST MODIFIED: June 2025 - Added comprehensive cross-condition individuality analysis
%                            with statistical comparisons and publication-quality outputs
%
% RELATED SCRIPTS:
% - plot_pose_individuality_analysis.m: Single-condition detailed individuality analysis
% - plot_tsnemap_difference_controls.m: Comparative analysis between specific conditions
% - compute_tsne_features.m: Feature extraction for behavioral clustering
% - cluster_tsne_maps.m: Behavioral sequence clustering and annotation
%
% ========================================================================

%% User Input Dialog - Analysis Mode Selection
clc;
logger('Starting cross-condition pose individuality analysis script', 'INFO');

% User input dialog for production/debug mode selection
analysis_mode_options = {'Production Mode', 'Debug Mode'};
[analysis_mode_idx, tf] = listdlg('PromptString', {'Select analysis mode:', '', ...
    'Production Mode: Full analysis with all conditions and export', ...
    'Debug Mode: Quick analysis with limited conditions for testing'}, ...
    'SelectionMode', 'single', ...
    'ListString', analysis_mode_options, ...
    'InitialValue', 1, ...
    'Name', 'Analysis Mode Selection', ...
    'ListSize', [400, 80]);

% Check if user cancelled
if tf == 0
    logger('Analysis cancelled by user', 'WARNING');
    return;
end

% Set mode and configure accordingly
is_production_mode = (analysis_mode_idx == 1);
analysis_mode = analysis_mode_options{analysis_mode_idx};

%% Initialization
close all;

GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Configure analysis based on selected mode
if is_production_mode
    logger('PRODUCTION MODE selected: Full analysis with all conditions', 'INFO');
    available_conditions = {'B', 'S', 'F', 'H', 'N'};
    condition_names = {'Baseline', 'Saline', 'Formalin', 'Sham', 'Neuropathic'};
    individuality_threshold = 0.80;  % 80% dominance threshold
    visualize = 'on';  % Show figures
    do_export = true;  % Export results
    export_high_res = true;  % High resolution exports
    perform_statistical_tests = true;  % Full statistical analysis
else
    logger('DEBUG MODE selected: Quick analysis for testing', 'INFO');
    available_conditions = {'B', 'S'};  % Only baseline and saline for quick testing
    condition_names = {'Baseline', 'Saline'};
    individuality_threshold = 0.75;  % Lower threshold for testing
    visualize = 'on';  % Show figures
    do_export = false;  % Skip exports in debug mode
    export_high_res = false;  % Skip high res exports
    perform_statistical_tests = false;  % Skip heavy statistical computations
end

logger(sprintf('Analysis mode: %s', analysis_mode), 'INFO');
logger(sprintf('Conditions to analyze: %s', strjoin(available_conditions, ', ')), 'INFO');

% Export folder
export_folder = fullfile(GC.temp_root, 'figs_presentation_painAI', 'cross_condition_individuality');
if ~exist(export_folder, 'dir')
    mkdir(export_folder);
end

logger('Configuration: analyzing all conditions with 80% individuality threshold', 'INFO');

%% Load Data
logger('Loading data', 'INFO');
% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');

% Load predictions
load(GC.filename_predictions, 'predictions', 'animal_condition_identifier');

% Load ratception structure
load(GC.filename_ratception, 'ratception_struct');

%% Prepare Global Data
logger('Preparing global data for individuality analysis', 'INFO');

% Extract animal list and conditions
upsamplig_factor = GC.repfactor;
long_animal_frames_identifier = repelem(animal_condition_identifier, upsamplig_factor);
animal_list_used_after_analysis = long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});

% Use GLOBAL cluster assignments and animal+condition identifiers
cluster_assignments = analysisstruct.annot_reordered{end,end}; % Final cluster assignments (ALL FRAMES)
global_animal_condition_ids = animal_list_used_after_analysis; % Full animal+condition identifiers
unique_clusters = unique(cluster_assignments);
unique_clusters = unique_clusters(unique_clusters > 0); % Remove background/noise clusters

% Get unique animal+condition combinations globally
global_unique_animal_condition_ids = unique(global_animal_condition_ids);

logger(sprintf('Found %d unique clusters and %d unique animal+condition combinations', ...
    length(unique_clusters), length(global_unique_animal_condition_ids)), 'INFO');

%% Global Individuality Analysis
logger('Performing GLOBAL individuality analysis (all conditions, all animals)', 'INFO');

% For each cluster, calculate animal+condition composition using GLOBAL dataset
cluster_animal_composition = cell(length(unique_clusters), 1);
cluster_dominant_animal = cell(length(unique_clusters), 1);
cluster_dominance_percentage = zeros(length(unique_clusters), 1);
cluster_is_individual = false(length(unique_clusters), 1);

logger('Analyzing cluster composition using GLOBAL dataset...', 'INFO');

for c_idx = 1:length(unique_clusters)
    cluster_id = unique_clusters(c_idx);
    cluster_frames = cluster_assignments == cluster_id;
    
    % Count frames per animal+condition combination in this cluster
    animal_condition_frame_counts = containers.Map();
    total_frames_in_cluster = sum(cluster_frames);
    
    for i = 1:length(global_animal_condition_ids)
        if cluster_frames(i)
            animal_condition_id = global_animal_condition_ids{i};
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
        
        cluster_dominant_animal{c_idx} = dominant_animal_condition;
        cluster_dominance_percentage(c_idx) = dominance_percentage;
        cluster_is_individual(c_idx) = dominance_percentage >= individuality_threshold;
        
        % Store full composition
        composition = struct();
        for j = 1:length(animal_condition_ids)
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

% Calculate summary statistics
total_clusters_analyzed = sum(cluster_dominance_percentage > 0);
individual_clusters = sum(cluster_is_individual);
individual_percentage = (individual_clusters / total_clusters_analyzed) * 100;

logger(sprintf('GLOBAL Summary: %d/%d clusters (%.1f%%) show individual dominance (>%.0f%%)', ...
    individual_clusters, total_clusters_analyzed, individual_percentage, individuality_threshold*100), 'INFO');

%% Count Clusters and Frames per Animal+Condition Combination
logger('Calculating per-animal+condition individuality statistics...', 'INFO');

% Initialize containers for cluster-level statistics
animal_condition_individual_clusters = containers.Map();
animal_condition_total_clusters_present = containers.Map();
animal_condition_total_dominant_clusters = containers.Map();

% Count all clusters where each animal+condition combination appears
for c_idx = 1:length(unique_clusters)
    if ~isempty(cluster_animal_composition{c_idx})
        composition = cluster_animal_composition{c_idx};
        field_names = fieldnames(composition);
        
        for f_idx = 1:length(field_names)
            field_name = field_names{f_idx};
            animal_condition_id = field_name(4:end); % Remove 'ID_' prefix
            animal_condition_id = regexprep(animal_condition_id, '__', '_');
            
            if isKey(animal_condition_total_clusters_present, animal_condition_id)
                animal_condition_total_clusters_present(animal_condition_id) = animal_condition_total_clusters_present(animal_condition_id) + 1;
            else
                animal_condition_total_clusters_present(animal_condition_id) = 1;
            end
        end
    end
end

% Count individual and dominant clusters per animal+condition combination
for c_idx = 1:length(unique_clusters)
    if ~strcmp(cluster_dominant_animal{c_idx}, 'None')
        animal_condition_id = cluster_dominant_animal{c_idx};
        
        % Count total dominant clusters
        if isKey(animal_condition_total_dominant_clusters, animal_condition_id)
            animal_condition_total_dominant_clusters(animal_condition_id) = animal_condition_total_dominant_clusters(animal_condition_id) + 1;
        else
            animal_condition_total_dominant_clusters(animal_condition_id) = 1;
        end
        
        % Count individual clusters
        if cluster_is_individual(c_idx)
            if isKey(animal_condition_individual_clusters, animal_condition_id)
                animal_condition_individual_clusters(animal_condition_id) = animal_condition_individual_clusters(animal_condition_id) + 1;
            else
                animal_condition_individual_clusters(animal_condition_id) = 1;
            end
        end
    end
end

%% Frame-level Analysis
logger('Calculating frame-level individuality statistics...', 'INFO');

% Initialize frame-level statistics containers
animal_condition_individual_frames = containers.Map();
animal_condition_total_frames = containers.Map();

% Count total frames per animal+condition combination
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
        
        % Add to individual frames count
        if isKey(animal_condition_individual_frames, dominant_animal_condition)
            animal_condition_individual_frames(dominant_animal_condition) = ...
                animal_condition_individual_frames(dominant_animal_condition) + frames_for_dominant;
        else
            animal_condition_individual_frames(dominant_animal_condition) = frames_for_dominant;
        end
    end
end

%% Organize Data by Condition
logger('Organizing data by experimental condition...', 'INFO');

% Initialize data structure for each condition
condition_data = struct();
for cond_idx = 1:length(available_conditions)
    condition = available_conditions{cond_idx};
    condition_data.(condition) = struct();
    condition_data.(condition).animal_ids = {};
    condition_data.(condition).cluster_individuality = [];
    condition_data.(condition).frame_individuality = [];
    condition_data.(condition).individual_clusters = [];
    condition_data.(condition).total_clusters_present = [];
    condition_data.(condition).individual_frames = [];
    condition_data.(condition).total_frames = [];
end

% Process each animal+condition combination
all_animal_condition_ids = keys(animal_condition_total_frames);

for i = 1:length(all_animal_condition_ids)
    animal_condition_id = all_animal_condition_ids{i};
    
    % Extract condition from animal+condition ID
    condition = animal_condition_id(end);
    
    % Skip if this condition is not in our analysis list
    if ~ismember(condition, available_conditions)
        continue;
    end
    
    % Extract animal ID (everything before the last '_')
    underscore_pos = find(animal_condition_id == '_', 1, 'last');
    if ~isempty(underscore_pos)
        animal_id = animal_condition_id(1:underscore_pos-1);
    else
        animal_id = animal_condition_id;
    end
    
    % Get statistics for this animal+condition combination
    individual_clusters = 0;
    if isKey(animal_condition_individual_clusters, animal_condition_id)
        individual_clusters = animal_condition_individual_clusters(animal_condition_id);
    end
    
    total_clusters_present = 0;
    if isKey(animal_condition_total_clusters_present, animal_condition_id)
        total_clusters_present = animal_condition_total_clusters_present(animal_condition_id);
    end
    
    individual_frames = 0;
    if isKey(animal_condition_individual_frames, animal_condition_id)
        individual_frames = animal_condition_individual_frames(animal_condition_id);
    end
    
    total_frames = 0;
    if isKey(animal_condition_total_frames, animal_condition_id)
        total_frames = animal_condition_total_frames(animal_condition_id);
    end
    
    % Calculate percentages
    cluster_individuality = 0;
    if total_clusters_present > 0
        cluster_individuality = (individual_clusters / total_clusters_present) * 100;
    end
    
    frame_individuality = 0;
    if total_frames > 0
        frame_individuality = (individual_frames / total_frames) * 100;
    end
    
    % Store in condition data
    condition_data.(condition).animal_ids{end+1} = animal_id;
    condition_data.(condition).cluster_individuality(end+1) = cluster_individuality;
    condition_data.(condition).frame_individuality(end+1) = frame_individuality;
    condition_data.(condition).individual_clusters(end+1) = individual_clusters;
    condition_data.(condition).total_clusters_present(end+1) = total_clusters_present;
    condition_data.(condition).individual_frames(end+1) = individual_frames;
    condition_data.(condition).total_frames(end+1) = total_frames;
end

% Log summary for each condition
for cond_idx = 1:length(available_conditions)
    condition = available_conditions{cond_idx};
    n_animals = length(condition_data.(condition).animal_ids);
    if n_animals > 0
        mean_cluster = mean(condition_data.(condition).cluster_individuality);
        mean_frame = mean(condition_data.(condition).frame_individuality);
        logger(sprintf('Condition %s: %d animals, Mean cluster individuality: %.1f%%, Mean frame individuality: %.1f%%', ...
            condition, n_animals, mean_cluster, mean_frame), 'INFO');
    else
        logger(sprintf('Condition %s: No animals found', condition), 'WARNING');
    end
end

%% Statistical Analysis
if perform_statistical_tests
    logger('Performing statistical analysis across conditions...', 'INFO');
else
    logger('Skipping statistical analysis (Debug Mode)', 'INFO');
end

% Prepare data for statistical tests
conditions_with_data = {};
cluster_individuality_by_condition = {};
frame_individuality_by_condition = {};
all_cluster_values = [];
all_frame_values = [];
condition_labels_cluster = {};
condition_labels_frame = {};

for cond_idx = 1:length(available_conditions)
    condition = available_conditions{cond_idx};
    if length(condition_data.(condition).cluster_individuality) > 0
        conditions_with_data{end+1} = condition;
        cluster_individuality_by_condition{end+1} = condition_data.(condition).cluster_individuality;
        frame_individuality_by_condition{end+1} = condition_data.(condition).frame_individuality;
        
        % For ANOVA (only if statistical tests are enabled)
        if perform_statistical_tests
            all_cluster_values = [all_cluster_values, condition_data.(condition).cluster_individuality];
            all_frame_values = [all_frame_values, condition_data.(condition).frame_individuality];
            
            n_animals = length(condition_data.(condition).cluster_individuality);
            condition_labels_cluster = [condition_labels_cluster, repmat({condition}, 1, n_animals)];
            condition_labels_frame = [condition_labels_frame, repmat({condition}, 1, n_animals)];
        end
    end
end

% Perform statistical tests if we have multiple conditions with data and statistical tests are enabled
stats_results = struct();
if length(conditions_with_data) > 1 && perform_statistical_tests
    logger('Performing ANOVA and post-hoc tests...', 'INFO');
    
    % One-way ANOVA for cluster individuality
    try
        [p_cluster, tbl_cluster, stats_cluster] = anova1(all_cluster_values, condition_labels_cluster, 'off');
        stats_results.cluster_anova_p = p_cluster;
        stats_results.cluster_anova_table = tbl_cluster;
        stats_results.cluster_anova_stats = stats_cluster;
        
        % Post-hoc tests for cluster individuality
        if p_cluster < 0.05
            [c_cluster, m_cluster] = multcompare(stats_cluster, 'Display', 'off');
            stats_results.cluster_posthoc = c_cluster;
            stats_results.cluster_means = m_cluster;
        end
        
        logger(sprintf('Cluster individuality ANOVA: F=%.3f, p=%.4f', tbl_cluster{2,5}, p_cluster), 'INFO');
    catch ME
        logger(sprintf('Error in cluster individuality ANOVA: %s', ME.message), 'WARNING');
    end
    
    % One-way ANOVA for frame individuality
    try
        [p_frame, tbl_frame, stats_frame] = anova1(all_frame_values, condition_labels_frame, 'off');
        stats_results.frame_anova_p = p_frame;
        stats_results.frame_anova_table = tbl_frame;
        stats_results.frame_anova_stats = stats_frame;
        
        % Post-hoc tests for frame individuality
        if p_frame < 0.05
            [c_frame, m_frame] = multcompare(stats_frame, 'Display', 'off');
            stats_results.frame_posthoc = c_frame;
            stats_results.frame_means = m_frame;
        end
        
        logger(sprintf('Frame individuality ANOVA: F=%.3f, p=%.4f', tbl_frame{2,5}, p_frame), 'INFO');
    catch ME
        logger(sprintf('Error in frame individuality ANOVA: %s', ME.message), 'WARNING');
    end
else
    logger('Insufficient conditions with data for statistical comparisons', 'WARNING');
end

%% Create Visualizations
logger('Creating visualizations...', 'INFO');

% Function to create main bar graph comparison
function create_condition_comparison_plots(condition_data, available_conditions, condition_names, stats_results, export_folder, do_export, perform_statistical_tests, export_high_res)
    
    % Calculate means and SEMs for each condition
    means_cluster = [];
    sems_cluster = [];
    means_frame = [];
    sems_frame = [];
    condition_labels = {};
    
    for cond_idx = 1:length(available_conditions)
        condition = available_conditions{cond_idx};
        if length(condition_data.(condition).cluster_individuality) > 0
            condition_labels{end+1} = condition_names{cond_idx};
            
            % Cluster individuality
            cluster_vals = condition_data.(condition).cluster_individuality;
            means_cluster(end+1) = mean(cluster_vals);
            sems_cluster(end+1) = std(cluster_vals) / sqrt(length(cluster_vals));
            
            % Frame individuality
            frame_vals = condition_data.(condition).frame_individuality;
            means_frame(end+1) = mean(frame_vals);
            sems_frame(end+1) = std(frame_vals) / sqrt(length(frame_vals));
        end
    end
    
    % Create main comparison figure
    fig = figure('Name', 'Cross-Condition Individuality Analysis', 'Position', [100, 100, 1400, 800], 'Color', 'w');
    
    % Subplot 1: Cluster Individuality Comparison
    subplot(2, 3, 1);
    bar_handle = bar(means_cluster, 'FaceColor', [0.3, 0.6, 1], 'EdgeColor', 'black', 'LineWidth', 1);
    hold on;
    errorbar(1:length(means_cluster), means_cluster, sems_cluster, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'XTickLabel', condition_labels, 'TickDir', 'out');
    ylabel('Cluster Individuality (%)');
    title('Cluster Individuality Across Conditions');
    xtickangle(45);
    ylim([0, max(means_cluster + sems_cluster) * 1.2]);
    
    % Add statistical significance annotations if available and statistical tests were performed
    if perform_statistical_tests && isfield(stats_results, 'cluster_anova_p') && stats_results.cluster_anova_p < 0.05
        text(0.5, max(means_cluster + sems_cluster) * 1.1, sprintf('ANOVA: p=%.4f', stats_results.cluster_anova_p), ...
            'FontSize', 10, 'FontWeight', 'bold');
    elseif ~perform_statistical_tests
        text(0.5, max(means_cluster + sems_cluster) * 1.1, 'Statistical tests disabled (Debug Mode)', ...
            'FontSize', 8, 'FontStyle', 'italic', 'Color', [0.5 0.5 0.5]);
    end
    
    % Subplot 2: Frame Individuality Comparison
    subplot(2, 3, 2);
    bar_handle = bar(means_frame, 'FaceColor', [1, 0.4, 0.4], 'EdgeColor', 'black', 'LineWidth', 1);
    hold on;
    errorbar(1:length(means_frame), means_frame, sems_frame, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'XTickLabel', condition_labels, 'TickDir', 'out');
    ylabel('Frame Individuality (%)');
    title('Frame Individuality Across Conditions');
    xtickangle(45);
    ylim([0, max(means_frame + sems_frame) * 1.2]);
    
    % Add statistical significance annotations if available and statistical tests were performed
    if perform_statistical_tests && isfield(stats_results, 'frame_anova_p') && stats_results.frame_anova_p < 0.05
        text(0.5, max(means_frame + sems_frame) * 1.1, sprintf('ANOVA: p=%.4f', stats_results.frame_anova_p), ...
            'FontSize', 10, 'FontWeight', 'bold');
    elseif ~perform_statistical_tests
        text(0.5, max(means_frame + sems_frame) * 1.1, 'Statistical tests disabled (Debug Mode)', ...
            'FontSize', 8, 'FontStyle', 'italic', 'Color', [0.5 0.5 0.5]);
    end
    
    % Subplot 3: Correlation between cluster and frame individuality
    subplot(2, 3, 3);
    all_cluster_vals = [];
    all_frame_vals = [];
    colors = lines(length(available_conditions));
    
    % Create mapping between available conditions and condition labels
    condition_label_idx = 1;
    for cond_idx = 1:length(available_conditions)
        condition = available_conditions{cond_idx};
        if length(condition_data.(condition).cluster_individuality) > 0
            cluster_vals = condition_data.(condition).cluster_individuality;
            frame_vals = condition_data.(condition).frame_individuality;
            
            scatter(cluster_vals, frame_vals, 50, colors(cond_idx, :), 'filled', 'DisplayName', condition_labels{condition_label_idx});
            hold on;
            
            all_cluster_vals = [all_cluster_vals, cluster_vals];
            all_frame_vals = [all_frame_vals, frame_vals];
            condition_label_idx = condition_label_idx + 1;
        end
    end
    
    % Add mean points with error bars for each condition
    condition_label_idx = 1;
    for cond_idx = 1:length(available_conditions)
        condition = available_conditions{cond_idx};
        if length(condition_data.(condition).cluster_individuality) > 0
            cluster_vals = condition_data.(condition).cluster_individuality;
            frame_vals = condition_data.(condition).frame_individuality;
            
            % Calculate means and SEMs
            mean_cluster = mean(cluster_vals);
            mean_frame = mean(frame_vals);
            sem_cluster = std(cluster_vals) / sqrt(length(cluster_vals));
            sem_frame = std(frame_vals) / sqrt(length(frame_vals));
            
            % Plot mean point with error bars
            errorbar(mean_cluster, mean_frame, sem_frame, sem_frame, sem_cluster, sem_cluster, ...
                'o', 'Color', colors(cond_idx, :), 'MarkerSize', 8, 'LineWidth', 2, ...
                'MarkerFaceColor', colors(cond_idx, :), 'MarkerEdgeColor', 'k', ...
                'HandleVisibility', 'off'); % Don't show in legend
            
            condition_label_idx = condition_label_idx + 1;
        end
    end
    
    % Add correlation line if we have data
    if length(all_cluster_vals) > 2
        [rho, p_corr] = corr(all_cluster_vals', all_frame_vals');
        coeffs = polyfit(all_cluster_vals, all_frame_vals, 1);
        x_fit = linspace(min(all_cluster_vals), max(all_cluster_vals), 100);
        y_fit = polyval(coeffs, x_fit);
        plot(x_fit, y_fit, 'k--', 'LineWidth', 2);
        
        text(0.1, 0.9, sprintf('r=%.3f, p=%.4f', rho, p_corr), 'Units', 'normalized', ...
            'FontSize', 10, 'FontWeight', 'bold');
    end
    
    xlabel('Cluster Individuality (%)');
    ylabel('Frame Individuality (%)');
    title('Cluster vs Frame Individuality');
    set(gca, 'TickDir', 'out');
    legend('Location', 'best');
    
    % Subplot 4: Individual animal trajectories (if multiple conditions per animal)
    subplot(2, 3, [4, 5, 6]);
    
    % Find animals that appear in multiple conditions
    all_animals = {};
    for cond_idx = 1:length(available_conditions)
        condition = available_conditions{cond_idx};
        all_animals = [all_animals, condition_data.(condition).animal_ids];
    end
    
    unique_animals = unique(all_animals);
    animals_multi_condition = {};
    
    for animal_idx = 1:length(unique_animals)
        animal_id = unique_animals{animal_idx};
        conditions_count = 0;
        
        for cond_idx = 1:length(available_conditions)
            condition = available_conditions{cond_idx};
            if ismember(animal_id, condition_data.(condition).animal_ids)
                conditions_count = conditions_count + 1;
            end
        end
        
        if conditions_count > 1
            animals_multi_condition{end+1} = animal_id;
        end
    end
    
    % Plot trajectories for animals with multiple conditions
    if ~isempty(animals_multi_condition)
        colors = lines(length(animals_multi_condition));
        
        for animal_idx = 1:min(10, length(animals_multi_condition)) % Limit to 10 animals for clarity
            animal_id = animals_multi_condition{animal_idx};
            frame_vals = [];
            condition_indices = [];
            
            for cond_idx = 1:length(available_conditions)
                condition = available_conditions{cond_idx};
                animal_position = find(strcmp(condition_data.(condition).animal_ids, animal_id));
                
                if ~isempty(animal_position)
                    frame_vals(end+1) = condition_data.(condition).frame_individuality(animal_position);
                    condition_indices(end+1) = cond_idx;
                end
            end
            
            if length(frame_vals) > 1
                plot(condition_indices, frame_vals, 'o-', 'Color', colors(animal_idx, :), ...
                    'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', animal_id);
                hold on;
            end
        end
        
        set(gca, 'XTick', 1:length(available_conditions), 'TickDir', 'out');
        % Only set labels for conditions that have data
        all_condition_labels = cell(1, length(available_conditions));
        label_idx = 1;
        for cond_idx = 1:length(available_conditions)
            condition = available_conditions{cond_idx};
            if length(condition_data.(condition).cluster_individuality) > 0
                all_condition_labels{cond_idx} = condition_labels{label_idx};
                label_idx = label_idx + 1;
            else
                all_condition_labels{cond_idx} = available_conditions{cond_idx};
            end
        end
        set(gca, 'XTickLabel', all_condition_labels);
        ylabel('Frame Individuality (%)');
        title('Individual Animal Frame Individuality Trajectories Across Conditions');
        legend('Location', 'best');
        xtickangle(45);
    else
        text(0.5, 0.5, 'No animals found across multiple conditions', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Units', 'normalized');
        title('Individual Animal Trajectories Across Conditions');
    end
    
    % Add main title
    sgtitle('Cross-Condition Pose Individuality Analysis', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Export if requested (production mode)
    if do_export
        fig_filename = fullfile(export_folder, 'cross_condition_individuality_analysis.png');
        if export_high_res
            % High resolution export for production
            saveas(fig, fig_filename);
            % Also save as PDF
            fig_filename_pdf = fullfile(export_folder, 'cross_condition_individuality_analysis.pdf');
            exportgraphics(fig, fig_filename_pdf, 'ContentType', 'vector', 'BackgroundColor', 'none');
            logger(sprintf('Saved high-res cross-condition analysis figure: %s', fig_filename), 'INFO');
        else
            saveas(fig, fig_filename);
            logger(sprintf('Saved cross-condition analysis figure: %s', fig_filename), 'INFO');
        end
    else
        logger('Figure export skipped (Debug Mode)', 'INFO');
    end
end

% Create the main visualization
create_condition_comparison_plots(condition_data, available_conditions, condition_names, stats_results, export_folder, do_export, perform_statistical_tests, export_high_res);

%% Additional Statistical Plots
if perform_statistical_tests && length(conditions_with_data) > 1
    logger('Creating additional statistical plots...', 'INFO');
else
    logger('Skipping additional statistical plots (Debug Mode or insufficient data)', 'INFO');
end

% Create detailed statistical comparison figure (only if statistical tests enabled)
if perform_statistical_tests && length(conditions_with_data) > 1
    fig_stats = figure('Name', 'Statistical Analysis Results', 'Position', [200, 200, 1200, 600], 'Color', 'w');
    
    % Subplot 1: Box plots for cluster individuality
    subplot(1, 2, 1);
    
    % Prepare data for boxplot
    cluster_data_for_boxplot = [];
    cluster_group_labels = [];
    
    for cond_idx = 1:length(conditions_with_data)
        condition = conditions_with_data{cond_idx};
        vals = condition_data.(condition).cluster_individuality;
        cluster_data_for_boxplot = [cluster_data_for_boxplot, vals];
        cluster_group_labels = [cluster_group_labels, repmat({condition}, 1, length(vals))];
    end
    
    boxplot(cluster_data_for_boxplot, cluster_group_labels);
    ylabel('Cluster Individuality (%)');
    title('Cluster Individuality Distribution by Condition');
    set(gca, 'TickDir', 'out');
    
    % Add sample sizes
    for cond_idx = 1:length(conditions_with_data)
        condition = conditions_with_data{cond_idx};
        n = length(condition_data.(condition).cluster_individuality);
        text(cond_idx, -5, sprintf('n=%d', n), 'HorizontalAlignment', 'center');
    end
    
    % Subplot 2: Box plots for frame individuality
    subplot(1, 2, 2);
    
    % Prepare data for boxplot
    frame_data_for_boxplot = [];
    frame_group_labels = [];
    
    for cond_idx = 1:length(conditions_with_data)
        condition = conditions_with_data{cond_idx};
        vals = condition_data.(condition).frame_individuality;
        frame_data_for_boxplot = [frame_data_for_boxplot, vals];
        frame_group_labels = [frame_group_labels, repmat({condition}, 1, length(vals))];
    end
    
    boxplot(frame_data_for_boxplot, frame_group_labels);
    ylabel('Frame Individuality (%)');
    title('Frame Individuality Distribution by Condition');
    set(gca, 'TickDir', 'out');
    
    % Add sample sizes
    for cond_idx = 1:length(conditions_with_data)
        condition = conditions_with_data{cond_idx};
        n = length(condition_data.(condition).frame_individuality);
        text(cond_idx, -5, sprintf('n=%d', n), 'HorizontalAlignment', 'center');
    end
    
    % Export if requested (production mode)
    if do_export
        stats_filename = fullfile(export_folder, 'cross_condition_statistical_analysis.png');
        if export_high_res
            % High resolution export for production
            saveas(fig_stats, stats_filename);
            logger(sprintf('Saved high-res statistical analysis figure: %s', stats_filename), 'INFO');
        else
            saveas(fig_stats, stats_filename);
            logger(sprintf('Saved statistical analysis figure: %s', stats_filename), 'INFO');
        end
    else
        logger('Statistical figure export skipped (Debug Mode)', 'INFO');
    end
end

%% Export Summary Data
if do_export
    logger('Exporting summary data...', 'INFO');
else
    logger('Data export skipped (Debug Mode)', 'INFO');
end

% Create comprehensive summary table
summary_table = table();
all_animal_condition_ids_export = {};
all_conditions_export = {};
all_animal_ids_export = {};
all_cluster_individuality_export = [];
all_frame_individuality_export = [];
all_individual_clusters_export = [];
all_total_clusters_present_export = [];
all_individual_frames_export = [];
all_total_frames_export = [];

for cond_idx = 1:length(available_conditions)
    condition = available_conditions{cond_idx};
    condition_name = condition_names{cond_idx};
    
    n_animals = length(condition_data.(condition).animal_ids);
    for animal_idx = 1:n_animals
        animal_id = condition_data.(condition).animal_ids{animal_idx};
        animal_condition_id = [animal_id '_' condition];
        
        all_animal_condition_ids_export{end+1} = animal_condition_id;
        all_conditions_export{end+1} = condition_name;
        all_animal_ids_export{end+1} = animal_id;
        all_cluster_individuality_export(end+1) = condition_data.(condition).cluster_individuality(animal_idx);
        all_frame_individuality_export(end+1) = condition_data.(condition).frame_individuality(animal_idx);
        all_individual_clusters_export(end+1) = condition_data.(condition).individual_clusters(animal_idx);
        all_total_clusters_present_export(end+1) = condition_data.(condition).total_clusters_present(animal_idx);
        all_individual_frames_export(end+1) = condition_data.(condition).individual_frames(animal_idx);
        all_total_frames_export(end+1) = condition_data.(condition).total_frames(animal_idx);
    end
end

% Create table
summary_table.Animal_Condition_ID = all_animal_condition_ids_export';
summary_table.Condition = all_conditions_export';
summary_table.Animal_ID = all_animal_ids_export';
summary_table.Cluster_Individuality_Percentage = all_cluster_individuality_export';
summary_table.Frame_Individuality_Percentage = all_frame_individuality_export';
summary_table.Individual_Clusters = all_individual_clusters_export';
summary_table.Total_Clusters_Present = all_total_clusters_present_export';
summary_table.Individual_Frames = all_individual_frames_export';
summary_table.Total_Frames = all_total_frames_export';

% Export summary table (production mode only)
if do_export
    summary_filename = fullfile(export_folder, 'cross_condition_individuality_summary.csv');
    writetable(summary_table, summary_filename);
    logger(sprintf('Saved summary table: %s', summary_filename), 'INFO');
    
    % Export statistical results if available
    if ~isempty(fieldnames(stats_results))
        stats_filename = fullfile(export_folder, 'cross_condition_statistical_results.mat');
        save(stats_filename, 'stats_results', 'condition_data');
        logger(sprintf('Saved statistical results: %s', stats_filename), 'INFO');
    end
else
    logger('Summary data export skipped (Debug Mode)', 'INFO');
end

%% Display Summary Results
logger('=== CROSS-CONDITION INDIVIDUALITY ANALYSIS RESULTS ===', 'INFO');
logger(sprintf('Individuality threshold: %.0f%% dominance', individuality_threshold * 100), 'INFO');
logger(sprintf('Total clusters analyzed: %d', total_clusters_analyzed), 'INFO');
logger(sprintf('Global individual clusters: %d (%.1f%%)', individual_clusters, individual_percentage), 'INFO');
logger(' ', 'INFO');

% Display condition-wise summary
logger('Condition-wise summary:', 'INFO');
for cond_idx = 1:length(available_conditions)
    condition = available_conditions{cond_idx};
    condition_name = condition_names{cond_idx};
    
    n_animals = length(condition_data.(condition).animal_ids);
    if n_animals > 0
        mean_cluster = mean(condition_data.(condition).cluster_individuality);
        std_cluster = std(condition_data.(condition).cluster_individuality);
        mean_frame = mean(condition_data.(condition).frame_individuality);
        std_frame = std(condition_data.(condition).frame_individuality);
        
        logger(sprintf('%s (%s): n=%d animals', condition_name, condition, n_animals), 'INFO');
        logger(sprintf('  Cluster individuality: %.1f ± %.1f%%', mean_cluster, std_cluster), 'INFO');
        logger(sprintf('  Frame individuality: %.1f ± %.1f%%', mean_frame, std_frame), 'INFO');
    else
        logger(sprintf('%s (%s): No animals', condition_name, condition), 'INFO');
    end
end

% Display statistical results
if perform_statistical_tests
    if isfield(stats_results, 'cluster_anova_p')
        logger(' ', 'INFO');
        logger('Statistical analysis results:', 'INFO');
        logger(sprintf('Cluster individuality ANOVA: p=%.4f', stats_results.cluster_anova_p), 'INFO');
        
        if stats_results.cluster_anova_p < 0.05
            logger('Significant differences found in cluster individuality across conditions', 'INFO');
        else
            logger('No significant differences in cluster individuality across conditions', 'INFO');
        end
    end

    if isfield(stats_results, 'frame_anova_p')
        logger(sprintf('Frame individuality ANOVA: p=%.4f', stats_results.frame_anova_p), 'INFO');
        
        if stats_results.frame_anova_p < 0.05
            logger('Significant differences found in frame individuality across conditions', 'INFO');
        else
            logger('No significant differences in frame individuality across conditions', 'INFO');
        end
    end
else
    logger(' ', 'INFO');
    logger('Statistical analysis skipped (Debug Mode)', 'INFO');
end

logger('Cross-condition pose individuality analysis complete', 'INFO');
logger(sprintf('Results exported to: %s', export_folder), 'INFO');
