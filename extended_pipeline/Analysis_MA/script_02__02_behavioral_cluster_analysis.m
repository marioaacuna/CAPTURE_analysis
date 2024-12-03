%% Behavioral Cluster Analysis Script
% This script is intended for the analysis of behavioral clusters derived from
% motion capture data. It takes the analysisstruct as input and computes all
% necessary statistics for comprehensive behavioral analysis.
%
% The script performs the following analyses:
% 1. Cluster composition analysis
% 2. Temporal analysis of cluster occurrences
% 3. Pose analysis within clusters
% 4. Kinematic feature analysis
% 5. Cluster stability analysis
% 6. Dimensionality reduction and visualization
% 7. Time series analysis of cluster sequences
% 8. Machine learning classification of conditions based on clusters
% 9. Network analysis of behavioral transitions
% 10. Entropy and complexity measures of behavioral sequences
%
% Inputs:
%   - analysisstruct: Structure containing t-SNE results and cluster assignments
%   - predictions: Structure containing raw motion capture predictions
%   - ratception_struct: Structure containing preprocessed motion capture data
%
% Outputs:
%   - Various figures, tables, and statistical results saved to the specified output directory

%% Initialization
clear;
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;


% File paths
% analysis_filename = fullfile(rootpath, 'raw_concat_analysis.mat');
% predictions_filename = fullfile(rootpath, 'agg_predictions.mat');
% ratception_filename = fullfile(rootpath, 'ratception_prediction.mat');

%% Load Data
% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');

% Load predictions
load(GC.filename_predictions, 'predictions', 'animal_condition_identifier');

% Load ratception structure
load(GC.filename_ratception, 'ratception_struct');

%% 1. Cluster Composition Analysis
disp('Performing Cluster Composition Analysis...');
%% Prepare Data
% Upsample animal_condition_identifier
input_params.repfactor = GC.repfactor;
% input_params = analysisstruct.input_params;
upsampled_identifiers = repelem(animal_condition_identifier, input_params.repfactor);

% Get the frames with good tracking
good_frames = analysisstruct.frames_with_good_tracking{1, 1};

% Extract cluster IDs and corresponding identifiers for good frames
cluster_ids = analysisstruct.annot_reordered{end};
frame_identifiers = upsampled_identifiers(good_frames);

% Extract unique animal IDs and conditions
[animal_ids, ~, animal_indices] = unique(cellfun(@(x) x(1:end-2), frame_identifiers, 'UniformOutput', false));
conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);

%% Analyze Cluster Proportions
unique_clusters = unique(cluster_ids);
num_clusters = length(unique_clusters);
num_animals = length(animal_ids);
unique_conditions = unique(conditions);
num_conditions = length(unique_conditions);

% Initialize cell array to store cluster proportions for each condition
cluster_proportions = cell(num_conditions, 1);
for i = 1:num_conditions
    cluster_proportions{i} = nan(num_animals, num_clusters);
end

% Calculate cluster proportions for each animal and condition
for i = 1:num_animals
    for c = 1:num_conditions
        animal_frames = animal_indices == i & strcmp(conditions, unique_conditions{c});
        total_frames = sum(animal_frames);
        
        if total_frames > 0  % Only calculate if there are frames for this condition
            for j = 1:num_clusters
                cluster_frames = cluster_ids(animal_frames) == unique_clusters(j);
                cluster_proportions{c}(i, j) = sum(cluster_frames) / total_frames;
            end
        end
    end
end

%% Perform Statistical Analysis
p_values_all = zeros(1, num_clusters);
f_stats = zeros(1, num_clusters);
condition_effects = cell(1, num_clusters);

% Define colors for each condition
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
          [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};

for j = 1:num_clusters
    % Prepare data for ANOVA
    cluster_data = [];
    group_labels = [];
    for c = 1:num_conditions
        valid_data = cluster_proportions{c}(:, j);
        valid_data = valid_data(~isnan(valid_data));
        cluster_data = [cluster_data; valid_data];
        group_labels = [group_labels; repmat(unique_conditions(c), length(valid_data), 1)];
    end
    
    % Perform one-way ANOVA
    [p_values_all(j), tbl] = anova1(cluster_data, group_labels, 'off');
    f_stats(j) = tbl{2,5};
    
    % Perform post-hoc multiple comparisons
    [~, ~, stats] = anova1(cluster_data, group_labels, 'off');
    condition_effects{j} = multcompare(stats, 'Display', 'off');
end

% Correct for multiple comparisons
[~, ~, ~, adj_p] = fdr_bh(p_values_all);

%% Visualize Proportions Results
fig_prop = figure('color', 'w', 'Position', [100 100 1500 700]);

barWidth = 0.15; % Adjusted for multiple conditions
gapWidth = 1;
currentX = 1;

for i = 1:length(unique_clusters)
    % Plot bars for each condition
    means = zeros(1, num_conditions);
    sems = zeros(1, num_conditions);
    
    for c = 1:num_conditions
        data = cluster_proportions{c}(:, i);
        means(c) = mean(data, 'omitnan');
        sems(c) = std(data, 'omitnan') / sqrt(sum(~isnan(data)));
        
        x_pos = currentX + (c-1)*barWidth;
        bar(x_pos, means(c), barWidth, 'FaceColor', colors{c});
        hold on;
        errorbar(x_pos, means(c), sems(c), 'k', 'LineStyle', 'none');
    end
    
    % Add significance markers if ANOVA shows significance
    if p_values_all(i) < 0.05
        y_max = max(means + sems) * 1.1;
        text(currentX + (num_conditions*barWidth)/2, y_max, '*', ...
             'HorizontalAlignment', 'center', 'FontSize', 22);
    end
    
    currentX = currentX + num_conditions*barWidth + gapWidth;
end

% Customize the plot
xlabel('Cluster ID', 'FontSize', 14);
ylabel('Mean Proportion', 'FontSize', 14);
title('Cluster Proportions Across Conditions', 'FontSize', 16);
legend(unique_conditions, 'Location', 'northeast', 'FontSize', 12);
set(gca, 'XTick', 1:barWidth*num_conditions+gapWidth:currentX-gapWidth, ...
    'XTickLabel', unique_clusters, 'FontSize', 12);
xlim([0, currentX-gapWidth]);
ylim([0, max(cell2mat(cellfun(@(x) max(max(x)), cluster_proportions, 'UniformOutput', false))) * 1.3]);
set(gca, 'TickDir', 'out');
grid off;
box off;

% Save the figure
saveas(fig_prop, fullfile(GC.figure_folder, 'cluster_proportions_comparison_all_conditions.fig'));

%% 2. Visualization of significant clusters
to_take = clusters(p_values_all < 0.05 & mean_diff > 0);
fig_predominant = figure('pos', [10,300,1500,1900]);
n_rows = ceil(sqrt(numel(to_take)));
n_cols = ceil(sqrt(numel(to_take)));

for ic = 1:numel(to_take)
    subplot(n_rows, n_cols, ic)
    this_cls = to_take(ic);
    fprintf('ic = %i - \n', this_cls)
    plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1},...
        find(analysisstruct.annot_reordered{end}==this_cls),['cl nr :  ', num2str(this_cls)]);
    title(this_cls)
end

%% Calculate predominant frames
% Initialize storage for cluster density per condition
cluster_density = zeros(length(unique_clusters), num_conditions);

% Initialize matrix to store results
clusterComposition = zeros(num_clusters, num_conditions);

for c = 1:length(unique_clusters)
    cluster_index = unique_clusters(c);
    
    % Find all frames belonging to the current cluster
    frames_in_cluster = find(cluster_ids == cluster_index);
    
    % Identify the animal-condition combinations for these frames
    animal_conditions_in_cluster = frame_identifiers(frames_in_cluster);
    
    % Count frames for each condition
    for cond = 1:num_conditions
        condition = unique_conditions{cond};
        cluster_density(c, cond) = sum(cellfun(@(x) endsWith(x, ['_' condition]), animal_conditions_in_cluster));
        clusterComposition(c, cond) = sum(cellfun(@(x) endsWith(x, ['_' condition]), animal_conditions_in_cluster));
    end
end

% Create a bar graph showing relative differences
figure('Position', [100 100 1200 600], 'Color', 'w');
b = bar(cluster_density, 'stacked');
for i = 1:num_conditions
    b(i).FaceColor = colors{i};
end
title('Cluster Frame Distribution Across Conditions');
xlabel('Cluster');
ylabel('Number of Frames');
legend(unique_conditions, 'Location', 'northeastoutside');
set(gca, 'XTick', 1:length(unique_clusters));
box off

% Calculate proportions
totalFrames = sum(clusterComposition, 2);
clusterProportions = clusterComposition ./ totalFrames;

% Identify predominantly associated clusters for each condition
threshold = 0.75;  % Adjusted threshold for multiple conditions
predominant_clusters = cell(num_conditions, 1);

for cond = 1:num_conditions
    predominant_clusters{cond} = find(clusterProportions(:, cond) >= threshold);
    
    % Create figure for predominant clusters of each condition
    if ~isempty(predominant_clusters{cond})
        fig_predominant = figure('pos', [10,10,2056,1350], 'color','w');
        to_take = predominant_clusters{cond};
        n_rows = ceil(sqrt(numel(to_take)));
        n_cols = ceil(sqrt(numel(to_take)));
        
        for ic = 1:numel(to_take)
            subplot(n_rows, n_cols, ic)
            this_cls = to_take(ic);
            fprintf('Condition %s - Cluster %i\n', unique_conditions{cond}, this_cls)
            try
                plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1},...
                    find(analysisstruct.annot_reordered{end}==this_cls),...
                    sprintf('Cluster %d (%s)', this_cls, unique_conditions{cond}));
                title(sprintf('Cluster %d\nProportion: %.2f', this_cls, clusterProportions(this_cls, cond)))
            catch
                continue
            end
        end
        
        sgtitle(sprintf('Predominant Clusters for Condition %s', unique_conditions{cond}))
        predominant_figure_name = fullfile(GC.figure_folder, sprintf('predominant_%s.pdf', unique_conditions{cond}));
        % export_fig(predominant_figure_name, '-pdf', fig_predominant)
    end
end

% Update results structure
results.cluster_density = cluster_density;
results.clusterProportions = clusterProportions;
results.predominant_clusters = predominant_clusters;

% Save updated results
save(fullfile(rootpath, 'cluster_analysis_results.mat'), 'results');


%% 2.1 feature analyses
%% Cluster Distribution Analysis Across Conditions
disp('Analyzing cluster distributions across conditions...');

% Load zvals
roothpath_CAPTURE = fullfile(rootpath, 'kinematics');
zvals_filename = fullfile(roothpath_CAPTURE, 'zvals.mat');
load(zvals_filename, 'zvals');

% Use zvals for analysis
features_reduced = zvals;
clusters = cluster_ids(:);

% Initialize structure for condition-specific data
condition_data = struct();

for c = 1:num_conditions
    cond_mask = strcmp(conditions, unique_conditions{c});
    
    condition_data(c).label = unique_conditions{c};
    condition_data(c).features = features_reduced(cond_mask, :);
    condition_data(c).clusters = clusters(cond_mask);
end
%% Visualization
figure('Position', [100 100 1200 800]);

% 1. t-SNE space distribution
subplot(2,2,1)
hold on;
colors = lines(num_conditions);
for c = 1:num_conditions
    scatter(condition_data(c).features(:,1), condition_data(c).features(:,2), ...
        5, colors(c,:), 'filled', 'AlphaData', 0.3);
end
xlabel('t-SNE 1'); ylabel('t-SNE 2');
title('Distribution in t-SNE Space');
legend(unique_conditions, 'Location', 'best');

% 2. Cluster frequency by condition
subplot(2,2,2)
cluster_freq = zeros(num_conditions, num_clusters);
for c = 1:num_conditions
    temp = histcounts(condition_data(c).clusters, 1:num_clusters+1);
    cluster_freq(c,:) = temp / sum(temp);
end
bar(cluster_freq', 'stacked');
xlabel('Cluster ID');
ylabel('Normalized Frequency');
title('Cluster Distribution by Condition');
legend(unique_conditions, 'Location', 'best');

% 3. Cluster similarity matrix
subplot(2,2,[3,4])
similarity_matrix = zeros(num_conditions);
for i = 1:num_conditions
    for j = 1:num_conditions
        p = cluster_freq(i,:);
        q = cluster_freq(j,:);
        similarity_matrix(i,j) = 1 - sqrt(0.5 * ...
            (kldiv(p + eps, q + eps) + kldiv(q + eps, p + eps)));
    end
end

% Display heatmap with values
heatmap_handle = heatmap(similarity_matrix, 'Colormap', jet, 'ColorbarVisible', 'on');
heatmap_handle.XDisplayLabels = unique_conditions;
heatmap_handle.YDisplayLabels = unique_conditions;
heatmap_handle.CellLabelFormat = '%.2f';
title('Condition Similarity Matrix');

% Helper function for KL divergence
function d = kldiv(p, q)
    d = sum(p .* log2(p./q));
end

%% 2.2 Transition Matrix Analysis
disp('Performing Transition Matrix Analysis...');

unique_animals = unique(animal_ids);

% Initialize a cell array to store transition matrices per condition
transition_matrices = cell(num_conditions, 1);
avg_offdiag_probs = cell(num_conditions, 1);


for cond_idx = 1:num_conditions
    condition = unique_conditions{cond_idx};
    condition_matrices = [];
    condition_avg_probs = [];

    
    for animal_idx = 1:length(unique_animals)
        animal = unique_animals{animal_idx};
        
        % Get indices for this animal and condition
        animal_cond_mask = strcmp(conditions, condition) & strcmp(animal_ids(animal_indices), animal);
        
        % Get cluster sequence for this animal and condition
        clusters_seq = cluster_ids(animal_cond_mask);
        
        if length(clusters_seq) > 1
            % Compute transition matrix for this sequence
            num_clusters = max(unique_clusters);
            T = zeros(num_clusters, num_clusters);
            
            for i = 1:length(clusters_seq)-1
                from = clusters_seq(i);
                to = clusters_seq(i+1);
                T(from, to) = T(from, to) + 1;
            end
            
            % Normalize transition matrix
            T = T ./ sum(T, 2);
            T(isnan(T)) = 0;

            % Extract off-diagonal elements excluding zeros
            off_diag_indices = ~eye(num_clusters);
            off_diag_elements = T(off_diag_indices & T > 0);
            if ~isempty(off_diag_elements)
                avg_prob = mean(off_diag_elements);
                % Store average off-diagonal probability
                condition_avg_probs = [condition_avg_probs; avg_prob];
            end
            
            % Store transition matrix
            condition_matrices = cat(3, condition_matrices, T);
        end
    end
    
    % Store all matrices for this condition
    transition_matrices{cond_idx} = condition_matrices;
    avg_offdiag_probs{cond_idx} = condition_avg_probs;

end

% Average transition matrices for each condition
avg_transition_matrices = cell(num_conditions,1);
for cond_idx = 1:num_conditions
    avg_transition_matrices{cond_idx} = mean(transition_matrices{cond_idx}, 3, 'omitnan');
end

% Visualize average transition matrices for each condition
figure('Position', [100 100 1200 600]);
for cond_idx = 1:num_conditions
    subplot(1, num_conditions, cond_idx);
    imagesc(avg_transition_matrices{cond_idx});
    colorbar;
    title(['Avg Transition Matrix - ' unique_conditions{cond_idx}]);
    xlabel('To Cluster');
    ylabel('From Cluster');
    set(gca, 'XTick', 1:num_clusters, 'YTick', 1:num_clusters);
    axis square;
    caxis([0 0.05])
end


% Perform statistical comparison across conditions
% Combine data and group labels
all_avg_probs = [];
group_labels = [];

for cond_idx = 1:num_conditions
    all_avg_probs = [all_avg_probs; avg_offdiag_probs{cond_idx}];
    group_labels = [group_labels; repmat(unique_conditions(cond_idx), length(avg_offdiag_probs{cond_idx}), 1)];
end

% Perform Kruskal-Wallis test
[p, tbl, stats] = kruskalwallis(all_avg_probs, group_labels, 'off');


% Display results
fprintf('Kruskal-Wallis test p-value: %.4f\n', p);

% Create figure for transition probability comparison
figure('Position', [100 100 800 600], 'Color', 'w');

% Calculate means and SEMs for each condition
means = zeros(1, num_conditions);
sems = zeros(1, num_conditions);
for cond_idx = 1:num_conditions
    means(cond_idx) = mean(avg_offdiag_probs{cond_idx});
    sems(cond_idx) = std(avg_offdiag_probs{cond_idx}) / sqrt(length(avg_offdiag_probs{cond_idx}));
end

% Create bar plot with error bars
b = bar(means, 'FaceColor', 'flat');
hold on;
errorbar(1:num_conditions, means, sems, 'k', 'LineStyle', 'none', 'CapSize', 10);

% Customize plot
for cond_idx = 1:num_conditions
    b.CData(cond_idx,:) = colors(cond_idx,:);
end

xlabel('Condition');
ylabel('Average Transition Probability');
title('Mean Transition Probabilities Across Conditions');
set(gca, 'XTick', 1:num_conditions, 'XTickLabel', unique_conditions);
box off;
set(gca, 'TickDir', 'out');

% Add significance marker if test is significant
if p < 0.05
    plot(1:num_conditions, max(means + sems) * 1.1 * ones(1,num_conditions), 'k-');
    text(mean(1:num_conditions), max(means + sems) * 1.15, sprintf('p = %.3f', p), ...
        'HorizontalAlignment', 'center');
end


%% 2. Temporal Analysis
disp('Performing Temporal Analysis...');
% TODO: Implement temporal analysis
% - Analyze temporal distribution of cluster occurrences
% - Compute transition probabilities between clusters
% - Perform Markov chain analysis on behavioral sequences

%% 3. Pose Analysis
disp('Performing Pose Analysis...');
% TODO: Implement pose analysis
% - Compute mean pose and variance for each cluster and condition
% - Perform statistical tests comparing poses between conditions within clusters
% - Visualize pose differences using vector fields or heatmaps

%% 4. Kinematic Feature Analysis
disp('Performing Kinematic Feature Analysis...');
% TODO: Implement kinematic feature analysis
% - Extract kinematic features (joint angles, velocities) for each cluster
% - Compare features between conditions using t-tests or ANOVAs
% - Perform discriminant analysis to identify distinguishing features

%% 5. Cluster Stability Analysis
disp('Performing Cluster Stability Analysis...');
% TODO: Implement cluster stability analysis
% - Assess stability of clusters across different animals within each condition
% - Use measures like Adjusted Rand Index or Normalized Mutual Information

%% 6. Dimensionality Reduction and Visualization
disp('Performing Dimensionality Reduction and Visualization...');
% TODO: Implement dimensionality reduction and visualization
% - Apply t-SNE or UMAP to visualize high-dimensional pose data
% - Color points by cluster and condition
% - Analyze distribution and overlap of conditions in reduced space

%% 7. Time Series Analysis
disp('Performing Time Series Analysis...');
% TODO: Implement time series analysis
% - Perform autocorrelation and cross-correlation analysis on cluster sequences
% - Apply change point detection algorithms to identify behavioral shifts

%% 8. Machine Learning Classification
disp('Performing Machine Learning Classification...');
% TODO: Implement machine learning classification
% - Train classifier to distinguish between conditions based on cluster occurrences
% - Analyze feature importance to identify key discriminative behaviors

%% 9. Network Analysis
disp('Performing Network Analysis...');
% TODO: Implement network analysis
% - Construct behavioral networks (nodes: clusters, edges: transitions)
% - Compare network properties between conditions

%% 10. Entropy and Complexity Measures
disp('Calculating Entropy and Complexity Measures...');
% TODO: Implement entropy and complexity analysis
% - Calculate entropy of cluster distributions for each condition
% - Apply complexity measures to cluster sequences

%% Save Results
disp('Saving Results...');
% TODO: Save all results, figures, and tables to output directory

disp('Analysis Complete!');