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

% Initialize matrices to store cluster proportions
cluster_proportions_S = zeros(num_animals, num_clusters);
cluster_proportions_F = zeros(num_animals, num_clusters);

% Calculate cluster proportions for each animal and condition
for i = 1:num_animals
    animal_frames_S = animal_indices == i & strcmp(conditions, 'S');
    animal_frames_F = animal_indices == i & strcmp(conditions, 'F');
    
    total_frames_S = sum(animal_frames_S);
    total_frames_F = sum(animal_frames_F);
    
    for j = 1:num_clusters
        cluster_frames_S = cluster_ids(animal_frames_S) == unique_clusters(j);
        cluster_frames_F = cluster_ids(animal_frames_F) == unique_clusters(j);
        
        cluster_proportions_S(i, j) = sum(cluster_frames_S) / total_frames_S;
        cluster_proportions_F(i, j) = sum(cluster_frames_F) / total_frames_F;
    end
end


%% Perform Statistical Analysis
p_values_all = zeros(1, num_clusters);
mean_diff = zeros(1, num_clusters);

for j = 1:num_clusters
    [h, p_values_all(j), ci, stats] = ttest(cluster_proportions_S(:, j), cluster_proportions_F(:, j));
    mean_diff(j) = mean(cluster_proportions_F(:, j) - cluster_proportions_S(:, j));
end

% Correct for multiple comparisons
% [h_corrected, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p_values);

%% Visualize Results
f1 = figure;
bar(mean_diff);
hold on;
% errorbar(1:num_clusters, mean_diff, adj_p, 'k', 'LineStyle', 'none');
xlabel('Cluster ID');
ylabel('Mean Difference (F - S)');
title('Difference in Cluster Proportions (F - S)');
sig_clusters = find(p_values_all < 0.05);
plot(sig_clusters, mean_diff(sig_clusters), 'r*');
hold off;


%% Visualize Proportions Results
%%
fig_prop = figure('color', 'w', 'Position', [100 100 1500 700]); % Create a rectangular figure

% Define some constants for plotting aesthetics
barWidth = 0.75; % Width of the bars
gapWidth = 1; % Gap between each cluster
currentX = 1; % Starting x position for the first bar

p_values = struct();
t_statistics = struct();

for i = 1:length(unique_clusters)
    cluster_name = ['id_',num2str(unique_clusters(i))];
    data_S = cluster_proportions_S(:, i);
    data_F = cluster_proportions_F(:, i);

    % Perform paired t-test (since we have before-after data for each animal)
    [h, p, ci, stats] = ttest(data_S, data_F);

    % Store p-value and t-statistic
    p_values.(cluster_name) = p;
    t_statistics.(cluster_name) = stats.tstat;

    % Calculate mean and SEM for S
    mean_S = mean(data_S);
    SEM_S = std(data_S) / sqrt(length(data_S));

    % Calculate mean and SEM for F
    mean_F = mean(data_F);
    SEM_F = std(data_F) / sqrt(length(data_F));

    % Plot bar for S
    bar(currentX, mean_S, barWidth, 'b');
    hold on;

    % Plot error bar for S
    errorbar(currentX, mean_S, SEM_S, 'k', 'LineStyle', 'none');

    % Plot bar for F right next to S
    bar(currentX + barWidth, mean_F, barWidth, 'r');

    % Plot error bar for F
    errorbar(currentX + barWidth, mean_F, SEM_F, 'k', 'LineStyle', 'none');

    % Add significance stars
    if p < 0.05 && p >= 0.01
        text(currentX + barWidth/2, max(mean_F, mean_S) + max(mean_F, mean_S) * 0.20, '*', 'HorizontalAlignment', 'center', 'FontSize',22)
    elseif p < 0.01 && p >= 0.001
        text(currentX + barWidth/2, max(mean_F, mean_S) + max(mean_F, mean_S) * 0.20, '**', 'HorizontalAlignment', 'center','FontSize',22)
    elseif p < 0.001
        text(currentX + barWidth/2, max(mean_F, mean_S) + max(mean_F, mean_S) * 0.20, '***', 'HorizontalAlignment', 'center','FontSize',22)
    end

    % Update x position for the next cluster
    currentX = currentX + barWidth * 2 + gapWidth;
end

% Customize the plot
xlabel('Cluster ID', 'FontSize', 14);
ylabel('Mean Proportion', 'FontSize', 14);
title('Cluster Proportions in Control (S) and Pain (F) Conditions', 'FontSize', 16);
% legend({'Control (S)', 'Pain (F)'}, 'Location', 'Best', 'FontSize', 12);
set(gca, 'XTick', 1:barWidth*2+gapWidth:currentX-gapWidth, 'XTickLabel', unique_clusters, 'FontSize', 12);
xlim([0, currentX-gapWidth]);
ylim([0, max(max(mean(cluster_proportions_S)), max(mean(cluster_proportions_F))) * 1.3]);  % Adjust y-axis to accommodate stars
set(gca, 'TickDir', 'out');
grid off;
box off;


% Adjust figure properties for better visibility
set(gcf, 'Color', 'w');  % Set figure background to white

hold off;

% Save the figure
% export_fig(fullfile(rootpath,'figures', 'cluster_proportions_comparison.pdf'), '-pdf', fig_prop)
saveas(fig_prop, fullfile(GC.figure_folder, 'cluster_proportions_comparison.fig'));



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

%% CAlculate predominant frames

% Initialize storage for the density of each cluster per condition
cluster_density_condition_1 = zeros(size(unique_clusters));
cluster_density_condition_0 = zeros(size(unique_clusters));

% Initialize matrices to store results
clusterComposition = zeros(num_clusters, 2);  % [S_count, F_count]

% animal_frames_ids = frame_identifiers(good_frames);
    
for c = 1:length(unique_clusters)
    cluster_index = unique_clusters(c);
    
    % Find all frames belonging to the current cluster
    frames_in_cluster = find(cluster_ids == cluster_index);
    
    % Identify the animal-condition combinations for these frames
    animal_conditions_in_cluster = frame_identifiers(frames_in_cluster);
    
    % Count frames for each condition
    cluster_density_condition_1(c) = sum(cellfun(@(x) endsWith(x, '_F'), animal_conditions_in_cluster));
    cluster_density_condition_0(c) = sum(cellfun(@(x) endsWith(x, '_S'), animal_conditions_in_cluster));
    
    
    clusterComposition(c, 1) =  sum(cellfun(@(x) endsWith(x, '_F'), animal_conditions_in_cluster));
    clusterComposition(c, 2) =  sum(cellfun(@(x) endsWith(x, '_S'), animal_conditions_in_cluster));
end

% Calculate the difference in number of frames between condition 1 and 0 for each cluster
difference_frames = cluster_density_condition_1 - cluster_density_condition_0;

% Create a bar graph
figure;
bar(difference_frames);
title('Difference in Frame Counts: Condition 1 vs. Condition 0');
xlabel('Cluster');
ylabel('Difference in Frame Counts');
set(gca, 'XTick', 0:10:length(unique_clusters), 'XTickLabel', arrayfun(@num2str, 0:10:length(unique_clusters), 'UniformOutput', false));
box off

totalFrames = sum(clusterComposition, 2);
clusterProportions = clusterComposition ./ totalFrames;

% Identify predominantly associated clusters
threshold = 0.99;  % Define threshold for "predominant" association
predominantF = find(clusterProportions(:, 1) >= threshold);
predominantS = find(clusterProportions(:, 2) >= threshold);

to_take = predominantF;
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


%% Save Results
results = struct();
results.cluster_proportions_S = cluster_proportions_S;
results.cluster_proportions_F = cluster_proportions_F;
results.p_values = p_values;
results.mean_differences = mean_diff;
results.predominantF = predominantF;
results.predominantS = predominantS;


save(fullfile(rootpath, 'cluster_analysis_results.mat'), 'results');


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