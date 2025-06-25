%% ========================================================================
%                   MATHEMATICAL PAIN POSES IDENTIFICATION AND PLOTTING
% ========================================================================
%
% DESCRIPTION:
% This script identifies "pain poses" by determining clusters that show 
% statistically significant differences in occurrence frequency between 
% pain conditions and their respective controls using mathematical approaches.
% It then generates visualizations of these pain-specific poses.
%
% PAIN POSE DEFINITION:
% A pain pose is a behavioral cluster where ≥90% of frames come from the pain 
% condition, and this dominance is statistically significant compared to random 
% chance (determined by permutation testing).
%
% MATHEMATICAL FORMULATION:
% Pain Pose = Cluster i where:
% 1. pain_dominance = pain_frames / (pain_frames + control_frames) ≥ 0.90
% 2. P(observed_dominance | random_assignment) < α (via permutation test)
%
% STATISTICAL APPROACH:
% - Binomial permutation test: simulates random frame assignment
% - Each frame has probability = global_pain_proportion of being "pain"
% - Tests if observed cluster dominance significantly exceeds chance expectation
% - P-value = proportion of simulations achieving ≥ observed dominance
% - Significance threshold: α = 0.01
%
% CONDITION COMPARISONS:
% - Formalin (F) vs Saline (S) - Chemical pain model
% - Neuropathic (N) vs Sham (H) - Neuropathic pain model
%
% WORKFLOW:
% 1. Load behavioral clustering data
% 2. Extract cluster frequencies for each condition
% 3. Perform statistical tests (Chi-square, effect size)
% 4. Apply multiple comparison corrections
% 5. Identify significantly different clusters as "pain poses"
% 6. Generate pose visualizations
% 7. Export results and figures
%
% OUTPUTS:
% - List of pain-specific clusters for each condition
% - Statistical test results
% - Pose visualization figures
% - Summary statistics and export files
%
% AUTHORS: CAPTURE Analysis Team
% CREATED: June 2025
% LAST MODIFIED: June 2025
%
% ========================================================================

%% Initialization
clc; close all;
logger('Starting mathematical pain poses identification', 'INFO');

% Load configurations
GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Statistical thresholds
alpha_level = 0.01;  % Significance level for permutation test
pain_dominance_threshold = 0.90;  % 90% dominance threshold for pain poses
n_permutations = 1000;  % Number of random shuffles for permutation test

% Export settings
export_folder = fullfile(GC.temp_root, 'figs_presentation_painAI', 'pain_poses_mathematical');
if ~exist(export_folder, 'dir')
    mkdir(export_folder);
end

logger('Configuration: dominance≥90%, permutations=1000, α=0.01', 'INFO');

%% Load Data
logger('Loading behavioral clustering data', 'INFO');

% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');

% Load predictions and animal identifiers
load(GC.filename_predictions, 'predictions', 'animal_condition_identifier');

% Extract cluster assignments and animal conditions
upsampling_factor = GC.repfactor;
long_animal_frames_identifier = repelem(animal_condition_identifier, upsampling_factor);
animal_list_used = long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});

% Get cluster assignments (final clustering result)
cluster_assignments = analysisstruct.annot_reordered{end,end};
unique_clusters = unique(cluster_assignments);
unique_clusters = unique_clusters(unique_clusters > 0); % Remove background/noise

logger(sprintf('Loaded %d frames across %d clusters', length(cluster_assignments), length(unique_clusters)), 'INFO');

%% Extract Condition-Specific Data
logger('Extracting condition-specific cluster frequencies', 'INFO');

% Define condition pairs for comparison
condition_pairs = {
    {'F', 'S', 'Formalin vs Saline'};
    {'N', 'H', 'Neuropathic vs Sham'}
};

% Initialize results structure
pain_poses_results = struct();

for pair_idx = 1:size(condition_pairs, 1)
    pain_condition = condition_pairs{pair_idx}{1};
    control_condition = condition_pairs{pair_idx}{2};
    comparison_name = condition_pairs{pair_idx}{3};
    
    logger(sprintf('Analyzing %s', comparison_name), 'INFO');
    
    % Extract frames for each condition
    pain_frames = contains(animal_list_used, ['_' pain_condition]);
    control_frames = contains(animal_list_used, ['_' control_condition]);
    
    % Get cluster assignments for each condition
    pain_clusters = cluster_assignments(pain_frames);
    control_clusters = cluster_assignments(control_frames);
    
    logger(sprintf('Pain condition (%s): %d frames, Control condition (%s): %d frames', ...
        pain_condition, sum(pain_frames), control_condition, sum(control_frames)), 'INFO');
    
    % Calculate cluster frequencies
    pain_cluster_counts = histcounts(pain_clusters, [unique_clusters, max(unique_clusters)+1]);
    control_cluster_counts = histcounts(control_clusters, [unique_clusters, max(unique_clusters)+1]);
    
    % Calculate proportions
    pain_total = sum(pain_cluster_counts);
    control_total = sum(control_cluster_counts);
    pain_proportions = pain_cluster_counts / pain_total;
    control_proportions = control_cluster_counts / control_total;    %% Statistical Testing - Permutation Test for Pain Dominance
    logger('Performing permutation tests for pain pose identification', 'INFO');
    
    % Pain pose criteria
    pain_dominance_threshold = 0.90;  % 90% of frames must be from pain condition
    n_permutations = 1000;            % Number of random shuffles
    
    % Initialize results arrays
    n_clusters = length(unique_clusters);
    observed_pain_dominance = zeros(n_clusters, 1);
    p_values = ones(n_clusters, 1);  % Initialize as non-significant
    effect_sizes = zeros(n_clusters, 1);  % Cohen's d for proportion differences
    is_pain_pose = false(n_clusters, 1);
    
    logger(sprintf('Pain dominance threshold: %.1f%%, Permutations: %d', ...
        pain_dominance_threshold*100, n_permutations), 'INFO');
    
    for c_idx = 1:n_clusters
        cluster_id = unique_clusters(c_idx);
        
        % Get counts for this cluster
        pain_frames = pain_cluster_counts(c_idx);
        control_frames = control_cluster_counts(c_idx);
        total_cluster_frames = pain_frames + control_frames;
        
        % Skip if cluster has too few frames
        if total_cluster_frames < 10
            observed_pain_dominance(c_idx) = 0;
            p_values(c_idx) = 1;
            effect_sizes(c_idx) = 0;  % No meaningful effect size for low-count clusters
            continue;
        end
        
        % Calculate observed pain dominance ratio
        observed_ratio = pain_frames / total_cluster_frames;
        observed_pain_dominance(c_idx) = observed_ratio;
        
        % Calculate effect size (Cohen's d for proportion differences)
        % Compare cluster proportion vs global proportion
        % Cohen's d interpretation: 0.2=small, 0.5=medium, 0.8=large effect
        global_pain_proportion = pain_total / (pain_total + control_total);
        
        % Calculate pooled standard error for proportions
        p1 = observed_ratio;  % Cluster proportion (pain/(pain+control) in this cluster)
        p2 = global_pain_proportion;  % Global proportion (overall pain/(pain+control))
        
        % Standard Cohen's d for proportions:
        % d = (p1 - p2) / sqrt((p1*(1-p1) + p2*(1-p2))/2)
        pooled_variance = (p1*(1-p1) + p2*(1-p2)) / 2;
        
        if pooled_variance > 0
            cohens_d = (p1 - p2) / sqrt(pooled_variance);
        else
            cohens_d = 0;  % No variance, no effect
        end
        
        effect_sizes(c_idx) = cohens_d;
        
        % Check if this cluster meets the pain pose criteria
        if observed_ratio >= pain_dominance_threshold
            
            % Perform permutation test
            % H0: Pain dominance is due to chance (random assignment of frames)
            % H1: Pain dominance is significantly higher than chance
              logger(sprintf('  Testing cluster %d: %.2f%% pain dominance (%d/%d frames)', ...
                cluster_id, observed_ratio*100, pain_frames, total_cluster_frames), 'DEBUG');
            
            % Calculate global pain proportion (baseline expectation)
            global_pain_proportion = pain_total / (pain_total + control_total);
            
            % Perform permutations using binomial sampling
            % H0: Each frame has probability = global_pain_proportion of being "pain"
            % H1: This cluster has significantly higher pain proportion than global average
            shuffled_pain_ratios = zeros(n_permutations, 1);
            
            for perm = 1:n_permutations
                % Simulate random assignment: each frame has global_pain_proportion 
                % chance of being labeled as "pain"
                random_pain_frames = binornd(total_cluster_frames, global_pain_proportion);
                shuffled_pain_ratios(perm) = random_pain_frames / total_cluster_frames;
            end

            % Calculate p-value: proportion of shuffles >= observed ratio
            % Combine exact ranking with conservative adjustment
            sorted_ratios = sort(shuffled_pain_ratios);
            exact_rank = sum(sorted_ratios < observed_ratio);
            ties = sum(sorted_ratios == observed_ratio);

            % Handle ties properly and add conservative adjustment
            p_val = (n_permutations - exact_rank - ties + 1) / (n_permutations + 1);
            p_values(c_idx) = p_val;
            
            % Mark as pain pose if significant
            is_pain_pose(c_idx) = (p_val < alpha_level);
              if p_val < alpha_level
                logger(sprintf('    → PAIN POSE: p=%.4f, observed=%.2f%%, expected=%.2f%%, mean_shuffle=%.2f%%', ...
                    p_val, observed_ratio*100, global_pain_proportion*100, mean(shuffled_pain_ratios)*100), 'INFO');
            end
            
        else
            % Doesn't meet dominance threshold
            p_values(c_idx) = 1;            is_pain_pose(c_idx) = false;
        end
    end
    
    % Identify pain poses
    significant_clusters = is_pain_pose;
    pain_clusters_identified = unique_clusters(significant_clusters);
      logger(sprintf('Found %d pain poses for %s (dominance≥%.0f%%, p<%.2f)', ...
        sum(significant_clusters), comparison_name, pain_dominance_threshold*100, alpha_level), 'INFO');
    
    % Calculate global pain proportion for this comparison
    global_pain_proportion = pain_total / (pain_total + control_total);
    logger(sprintf('Global %s proportion: %.2f%% (%d/%d total frames)', ...
        pain_condition, global_pain_proportion*100, pain_total, pain_total + control_total), 'INFO');
    
    % Store results
    results_field = [pain_condition '_vs_' control_condition];
    pain_poses_results.(results_field) = struct();
    pain_poses_results.(results_field).pain_condition = pain_condition;
    pain_poses_results.(results_field).control_condition = control_condition;
    pain_poses_results.(results_field).comparison_name = comparison_name;
    pain_poses_results.(results_field).pain_clusters = pain_clusters_identified;
    pain_poses_results.(results_field).all_clusters = unique_clusters;
    pain_poses_results.(results_field).observed_pain_dominance = observed_pain_dominance;
    pain_poses_results.(results_field).p_values = p_values;
    pain_poses_results.(results_field).effect_sizes = effect_sizes;
    pain_poses_results.(results_field).is_pain_pose = is_pain_pose;
    pain_poses_results.(results_field).significant_clusters = significant_clusters;
    pain_poses_results.(results_field).pain_proportions = pain_proportions;    pain_poses_results.(results_field).control_proportions = control_proportions;
    pain_poses_results.(results_field).pain_dominance_threshold = pain_dominance_threshold;
    pain_poses_results.(results_field).n_permutations = n_permutations;
    pain_poses_results.(results_field).global_pain_proportion = global_pain_proportion;
      % Log top significant clusters
    if sum(significant_clusters) > 0
        [sorted_dominance, sort_idx] = sort(observed_pain_dominance(significant_clusters), 'descend');
        sig_cluster_ids = unique_clusters(significant_clusters);
        sig_cluster_ids = sig_cluster_ids(sort_idx);
        
        logger(sprintf('Top pain poses for %s:', comparison_name), 'INFO');
        for i = 1:min(5, length(sig_cluster_ids))
            cluster_id = sig_cluster_ids(i);
            cluster_idx = find(unique_clusters == cluster_id);
            logger(sprintf('  Cluster %d: p=%.4f, d=%.2f, dominance=%.1f%%, pain_freq=%.3f, control_freq=%.3f', ...
                cluster_id, p_values(cluster_idx), effect_sizes(cluster_idx), observed_pain_dominance(cluster_idx)*100, ...
                pain_proportions(cluster_idx), control_proportions(cluster_idx)), 'INFO');
        end
    end
end

%% Identify Common Pain Poses Across Modalities
logger('Identifying common pain poses across pain modalities', 'INFO');

% Get pain poses from each condition
field_names = fieldnames(pain_poses_results);
F_vs_S_poses = [];
N_vs_H_poses = [];

for field_idx = 1:length(field_names)
    field_name = field_names{field_idx};
    
    % Skip the common_across_modalities field as it has different structure
    if strcmp(field_name, 'common_across_modalities')
        continue;
    end
    
    results = pain_poses_results.(field_name);
    
    if strcmp(results.pain_condition, 'F')
        F_vs_S_poses = results.pain_clusters;
        logger(sprintf('Formalin pain poses: %s', mat2str(F_vs_S_poses)), 'INFO');
    elseif strcmp(results.pain_condition, 'N')
        N_vs_H_poses = results.pain_clusters;
        logger(sprintf('Neuropathic pain poses: %s', mat2str(N_vs_H_poses)), 'INFO');
    end
end

% Find common pain poses
common_pain_poses = intersect(F_vs_S_poses, N_vs_H_poses);

if ~isempty(common_pain_poses)
    logger(sprintf('Found %d COMMON PAIN POSES across both modalities: %s', ...
        length(common_pain_poses), mat2str(common_pain_poses)), 'INFO');
    
    % Store common pain poses results
    pain_poses_results.common_across_modalities = struct();
    pain_poses_results.common_across_modalities.common_pain_poses = common_pain_poses;
    pain_poses_results.common_across_modalities.formalin_poses = F_vs_S_poses;
    pain_poses_results.common_across_modalities.neuropathic_poses = N_vs_H_poses;
    pain_poses_results.common_across_modalities.n_common = length(common_pain_poses);
    pain_poses_results.common_across_modalities.n_formalin_only = length(setdiff(F_vs_S_poses, common_pain_poses));
    pain_poses_results.common_across_modalities.n_neuropathic_only = length(setdiff(N_vs_H_poses, common_pain_poses));
    
    % Log detailed analysis
    logger('=== CROSS-MODALITY PAIN POSE ANALYSIS ===', 'INFO');
    logger(sprintf('Total Formalin pain poses: %d', length(F_vs_S_poses)), 'INFO');
    logger(sprintf('Total Neuropathic pain poses: %d', length(N_vs_H_poses)), 'INFO');
    logger(sprintf('Common pain poses: %d', length(common_pain_poses)), 'INFO');
    logger(sprintf('Formalin-specific poses: %d', length(setdiff(F_vs_S_poses, common_pain_poses))), 'INFO');
    logger(sprintf('Neuropathic-specific poses: %d', length(setdiff(N_vs_H_poses, common_pain_poses))), 'INFO');
    
    if ~isempty(setdiff(F_vs_S_poses, common_pain_poses))
        logger(sprintf('Formalin-only poses: %s', mat2str(setdiff(F_vs_S_poses, common_pain_poses))), 'INFO');
    end
    if ~isempty(setdiff(N_vs_H_poses, common_pain_poses))
        logger(sprintf('Neuropathic-only poses: %s', mat2str(setdiff(N_vs_H_poses, common_pain_poses))), 'INFO');
    end
    
else
    logger('No common pain poses found across modalities', 'WARNING');
    pain_poses_results.common_across_modalities = struct();
    pain_poses_results.common_across_modalities.common_pain_poses = [];
    pain_poses_results.common_across_modalities.formalin_poses = F_vs_S_poses;
    pain_poses_results.common_across_modalities.neuropathic_poses = N_vs_H_poses;
    pain_poses_results.common_across_modalities.n_common = 0;
    pain_poses_results.common_across_modalities.n_formalin_only = length(F_vs_S_poses);
    pain_poses_results.common_across_modalities.n_neuropathic_only = length(N_vs_H_poses);
end

%% Generate Pain Pose Visualizations
logger('Generating pain pose visualizations', 'INFO');

% Plot pain poses for each condition comparison
field_names = fieldnames(pain_poses_results);

for field_idx = 1:length(field_names)
    field_name = field_names{field_idx};
    
    % Skip the common_across_modalities field as it has different structure
    if strcmp(field_name, 'common_across_modalities')
        continue;
    end
    
    results = pain_poses_results.(field_name);
    
    if isempty(results.pain_clusters)
        logger(sprintf('No significant pain poses found for %s', results.comparison_name), 'WARNING');
        continue;
    end
    
    % Create figure for this comparison
    pain_condition = results.pain_condition;
    cls = results.pain_clusters;
    
    fig_name = sprintf('Pain_Poses_%s', field_name);
    fig_poses = figure('Name', fig_name, 'Position', [10, 300, 1500, 1900]);
    
    nclus = numel(cls);
    n_rows = ceil(sqrt(nclus));
    n_cols = ceil(sqrt(nclus));
    
    logger(sprintf('Plotting %d pain poses for %s in %dx%d grid', ...
        nclus, results.comparison_name, n_rows, n_cols), 'INFO');
    
    for ic = 1:numel(cls)
        subplot(n_rows, n_cols, ic);
        this_cls = cls(ic);
        
        % Find frames belonging to this cluster
        cluster_frames = find(cluster_assignments == this_cls);
        
        % Get statistics for this cluster
        cluster_idx = find(results.all_clusters == this_cls);
        p_val = results.p_values(cluster_idx);
        effect_size = results.effect_sizes(cluster_idx);
        pain_freq = results.pain_proportions(cluster_idx);
        control_freq = results.control_proportions(cluster_idx);
        
        % Plot the mean cluster pose
        try
            plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1}, ...
                cluster_frames, ['Cluster ' num2str(this_cls)]);
        catch ME
            logger(sprintf('Error plotting cluster %d: %s', this_cls, ME.message), 'WARNING');
            % Create empty plot with error message
            text(0.5, 0.5, sprintf('Error plotting\nCluster %d', this_cls), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        
        % Add detailed title with statistics
        title_str = sprintf('Cluster %d\np=%.3f, d=%.2f\nPain: %.3f, Ctrl: %.3f', ...
            this_cls, p_val, effect_size, pain_freq, control_freq);
        title(title_str, 'FontSize', 10);
    end
    
    % Add main title
    main_title = sprintf('Pain Poses: %s (n=%d poses)', results.comparison_name, nclus);
    sgtitle(main_title, 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save figure
    cluster_poses_figure_filename = fullfile(export_folder, [fig_name '.pdf']);
    try
        exportgraphics(fig_poses, cluster_poses_figure_filename, 'Resolution', 300);
        logger(sprintf('Saved pain poses figure: %s', cluster_poses_figure_filename), 'INFO');
    catch ME
        logger(sprintf('Error saving figure: %s', ME.message), 'WARNING');
    end
    
    % % Also save as PNG
    % png_filename = fullfile(export_folder, [fig_name '.png']);
    % try
    %     exportgraphics(fig_poses, png_filename, 'Resolution', 300);
    % catch ME
    %     logger(sprintf('Error saving PNG: %s', ME.message), 'WARNING');    end
end

%% Generate Common Pain Poses Visualization
if isfield(pain_poses_results, 'common_across_modalities') && ...
   ~isempty(pain_poses_results.common_across_modalities.common_pain_poses)
    
    logger('Generating common pain poses visualization', 'INFO');
    
    common_poses = pain_poses_results.common_across_modalities.common_pain_poses;
    
    % Create figure for common pain poses
    fig_common = figure('Name', 'Common_Pain_Poses', 'Position', [10, 300, 1500, 1200]);
    
    nclus = length(common_poses);
    n_rows = 6;%ceil(sqrt(nclus));
    n_cols = 3;%ceil(sqrt(nclus));
    
    logger(sprintf('Plotting %d common pain poses in %dx%d grid', nclus, n_rows, n_cols), 'INFO');
    
    for ic = 1:nclus
        subplot(n_rows, n_cols, ic);
        this_cls = common_poses(ic);
        
        % Find frames belonging to this cluster
        cluster_frames = find(cluster_assignments == this_cls);
        
        % Plot the mean cluster pose
        try
            plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1}, ...
                cluster_frames, ['Common Pain Pose ' num2str(this_cls)]);
        catch ME
            logger(sprintf('Error plotting common pose %d: %s', this_cls, ME.message), 'WARNING');
            % Create empty plot with error message
            text(0.5, 0.5, sprintf('Error plotting\nCommon Pose %d', this_cls), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        axis square
        
        % Add title
        title(sprintf('Cluster %d\n(Common Pain Pose)', this_cls), 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    % Add main title
    main_title = sprintf('COMMON PAIN POSES (n=%d)\nPresent in both Formalin and Neuropathic Pain', nclus);
    sgtitle(main_title, 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'red');
    
    % Save common pain poses figure
    common_figure_filename = fullfile(export_folder, 'Common_Pain_Poses.pdf');
    try
        exportgraphics(fig_common, common_figure_filename, 'Resolution', 300);
        logger(sprintf('Saved common pain poses figure: %s', common_figure_filename), 'INFO');
    catch ME
        logger(sprintf('Error saving common poses figure: %s', ME.message), 'WARNING');
    end
    
    % % Also save as PNG
    % common_png_filename = fullfile(export_folder, 'Common_Pain_Poses.png');
    % try
    %     exportgraphics(fig_common, common_png_filename, 'Resolution', 300);
% Also save as PNG
    % common_png_filename = fullfile(export_folder, 'Common_Pain_Poses.png');
    % try
    %     exportgraphics(fig_common, common_png_filename, 'Resolution', 300);
    % catch ME
    %     logger(sprintf('Error saving common poses PNG: %s', ME.message), 'WARNING');
    % endfik    % catch ME
    %     logger(sprintf('Error saving common poses PNG: %s', ME.message), 'WARNING');
    % end
    
else
    logger('No common pain poses to visualize', 'INFO');
end

%% Generate Summary Statistics and Export
logger('Generating summary statistics and exports', 'INFO');

% Create comprehensive summary table
summary_table = table();
row_idx = 1;

for field_idx = 1:length(field_names)
    field_name = field_names{field_idx};
    
    % Skip the common_across_modalities field as it has different structure
    if strcmp(field_name, 'common_across_modalities')
        continue;
    end
    
    results = pain_poses_results.(field_name);
    
    for c_idx = 1:length(results.all_clusters)
        cluster_id = results.all_clusters(c_idx);
          summary_table.Comparison{row_idx} = results.comparison_name;
        summary_table.Pain_Condition{row_idx} = results.pain_condition;
        summary_table.Control_Condition{row_idx} = results.control_condition;
        summary_table.Cluster_ID(row_idx) = cluster_id;
        summary_table.Pain_Dominance_Ratio(row_idx) = results.observed_pain_dominance(c_idx);
        summary_table.P_Value(row_idx) = results.p_values(c_idx);
        summary_table.Effect_Size_Cohens_d(row_idx) = results.effect_sizes(c_idx);
        summary_table.Pain_Frequency(row_idx) = results.pain_proportions(c_idx);        summary_table.Control_Frequency(row_idx) = results.control_proportions(c_idx);
        summary_table.Is_Pain_Pose(row_idx) = results.is_pain_pose(c_idx);
        summary_table.Dominance_Threshold(row_idx) = results.pain_dominance_threshold;
        summary_table.N_Permutations(row_idx) = results.n_permutations;
        
        row_idx = row_idx + 1;
    end
end

% Save summary table
summary_filename = fullfile(export_folder, 'pain_poses_statistical_summary.csv');
writetable(summary_table, summary_filename);
logger(sprintf('Saved summary table: %s', summary_filename), 'INFO');

% Save results structure
results_filename = fullfile(export_folder, 'pain_poses_results.mat');
save(results_filename, 'pain_poses_results', 'summary_table', 'alpha_level', ...
    'pain_dominance_threshold', 'n_permutations');
logger(sprintf('Saved results structure: %s', results_filename), 'INFO');

%% Generate Statistical Summary Report
logger('Generating statistical summary report', 'INFO');

% Create text report
report_filename = fullfile(export_folder, 'pain_poses_statistical_report.txt');
fid = fopen(report_filename, 'w');

fprintf(fid, 'MATHEMATICAL PAIN POSES IDENTIFICATION - STATISTICAL REPORT\n');
fprintf(fid, '==========================================================\n\n');
fprintf(fid, 'Analysis Date: %s\n', datestr(now));
fprintf(fid, 'Method: Permutation Test for Pain Dominance\n');
fprintf(fid, 'Statistical Parameters:\n');
fprintf(fid, '  - Significance level (α): %.3f\n', alpha_level);
fprintf(fid, '  - Pain dominance threshold: %.0f%%\n', pain_dominance_threshold*100);
fprintf(fid, '  - Number of permutations: %d\n', n_permutations);
fprintf(fid, '  - Effect size metric: Cohen''s d (0.2=small, 0.5=medium, 0.8=large)\n');
fprintf(fid, '\n');

for field_idx = 1:length(field_names)
    field_name = field_names{field_idx};
    
    % Skip the common_across_modalities field as it has different structure
    if strcmp(field_name, 'common_across_modalities')
        continue;
    end
    
    results = pain_poses_results.(field_name);
      fprintf(fid, 'COMPARISON: %s\n', results.comparison_name);
    fprintf(fid, '----------------------------------------\n');
    fprintf(fid, 'Pain condition: %s\n', results.pain_condition);
    fprintf(fid, 'Control condition: %s\n', results.control_condition);
    fprintf(fid, 'Pain dominance threshold: %.0f%%\n', results.pain_dominance_threshold*100);
    fprintf(fid, 'Number of permutations: %d\n', results.n_permutations);
    fprintf(fid, 'Total clusters tested: %d\n', length(results.all_clusters));
    fprintf(fid, 'Pain poses identified: %d\n', sum(results.is_pain_pose));
    fprintf(fid, '\n');
    
    if sum(results.is_pain_pose) > 0
        fprintf(fid, 'IDENTIFIED PAIN POSES:\n');
        sig_indices = find(results.is_pain_pose);
        for i = 1:length(sig_indices)
            idx = sig_indices(i);
            cluster_id = results.all_clusters(idx);
            fprintf(fid, '  Cluster %d: p=%.4f, d=%.2f, dominance=%.1f%%, pain_freq=%.3f, control_freq=%.3f\n', ...
                cluster_id, results.p_values(idx), results.effect_sizes(idx), results.observed_pain_dominance(idx)*100, ...
                results.pain_proportions(idx), results.control_proportions(idx));
        end
    else
        fprintf(fid, 'No pain poses identified.\n');    end
    fprintf(fid, '\n');
end

% Add common pain poses analysis to report
if isfield(pain_poses_results, 'common_across_modalities')
    fprintf(fid, '=======================================================\n');
    fprintf(fid, 'CROSS-MODALITY PAIN POSE ANALYSIS\n');
    fprintf(fid, '=======================================================\n\n');
    
    common_results = pain_poses_results.common_across_modalities;
    
    fprintf(fid, 'Summary:\n');
    fprintf(fid, '  - Formalin pain poses: %d\n', length(common_results.formalin_poses));
    fprintf(fid, '  - Neuropathic pain poses: %d\n', length(common_results.neuropathic_poses));
    fprintf(fid, '  - Common pain poses: %d\n', common_results.n_common);
    fprintf(fid, '  - Formalin-specific poses: %d\n', common_results.n_formalin_only);
    fprintf(fid, '  - Neuropathic-specific poses: %d\n', common_results.n_neuropathic_only);
    fprintf(fid, '\n');
    
    if ~isempty(common_results.common_pain_poses)
        fprintf(fid, 'COMMON PAIN POSES (present in both modalities):\n');
        for i = 1:length(common_results.common_pain_poses)
            fprintf(fid, '  Cluster %d\n', common_results.common_pain_poses(i));
        end
        fprintf(fid, '\n');
    end
    
    if common_results.n_formalin_only > 0
        formalin_only = setdiff(common_results.formalin_poses, common_results.common_pain_poses);
        fprintf(fid, 'FORMALIN-SPECIFIC PAIN POSES:\n');
        for i = 1:length(formalin_only)
            fprintf(fid, '  Cluster %d\n', formalin_only(i));
        end
        fprintf(fid, '\n');
    end
    
    if common_results.n_neuropathic_only > 0
        neuropathic_only = setdiff(common_results.neuropathic_poses, common_results.common_pain_poses);
        fprintf(fid, 'NEUROPATHIC-SPECIFIC PAIN POSES:\n');
        for i = 1:length(neuropathic_only)
            fprintf(fid, '  Cluster %d\n', neuropathic_only(i));
        end
        fprintf(fid, '\n');
    end
end

fclose(fid);
logger(sprintf('Saved statistical report: %s', report_filename), 'INFO');

%% Final Summary
logger('=== PAIN POSES ANALYSIS COMPLETE ===', 'INFO');
logger(sprintf('Results saved to: %s', export_folder), 'INFO');

total_pain_poses = 0;
for field_idx = 1:length(field_names)
    field_name = field_names{field_idx};
    if strcmp(field_name, 'common_across_modalities')
        continue; % Skip the common analysis summary
    end
    results = pain_poses_results.(field_name);
    n_poses = sum(results.is_pain_pose);
    total_pain_poses = total_pain_poses + n_poses;
    logger(sprintf('%s: %d pain poses identified', results.comparison_name, n_poses), 'INFO');
end

logger(sprintf('TOTAL PAIN POSES IDENTIFIED: %d', total_pain_poses), 'INFO');

% Add common pain poses summary
if isfield(pain_poses_results, 'common_across_modalities')
    common_results = pain_poses_results.common_across_modalities;
    logger('=== CROSS-MODALITY ANALYSIS ===', 'INFO');
    logger(sprintf('Common pain poses (both modalities): %d', common_results.n_common), 'INFO');
    logger(sprintf('Formalin-specific poses: %d', common_results.n_formalin_only), 'INFO');
    logger(sprintf('Neuropathic-specific poses: %d', common_results.n_neuropathic_only), 'INFO');
    
    if ~isempty(common_results.common_pain_poses)
        logger(sprintf('ROBUST PAIN POSES: %s', mat2str(common_results.common_pain_poses)), 'INFO');
    end
end
logger('Files generated:', 'INFO');
logger(sprintf('  - Pain pose figures: %s', export_folder), 'INFO');
logger(sprintf('  - Statistical summary: %s', summary_filename), 'INFO');
logger(sprintf('  - Results structure: %s', results_filename), 'INFO');
logger(sprintf('  - Statistical report: %s', report_filename), 'INFO');

%% Logger Function
function logger(message, level)
    % Simple logging function
    if nargin < 2
        level = 'INFO';
    end
    
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('[%s] %s: %s\n', timestamp, level, message);
end
