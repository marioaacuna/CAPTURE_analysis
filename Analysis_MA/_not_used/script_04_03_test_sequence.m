%% Initialization
clear;
close all;
clc;

% Set root path based on operating system
if ispc
    rootpath = 'D:\CAPTURE';
else
    rootpath = '/Users/mario/Library/Mobile Documents/com~apple~CloudDocs/_todo/SpontaneousPainProject/240704_data';
end

% File paths
analysis_filename = fullfile(rootpath, 'raw_concat_analysis.mat');
predictions_filename = fullfile(rootpath, 'agg_predictions.mat');
ratception_filename = fullfile(rootpath, 'ratception_prediction.mat');

%% Load Data
% Load analysis structure
load(analysis_filename, 'analysisstruct');
% Load predictions
load(predictions_filename, 'predictions', 'animal_condition_identifier');
% Load ratception structure
load(ratception_filename, 'ratception_struct');

%% Prepare Data
% Upsample animal_condition_identifier
input_params.repfactor = 3;
upsampled_identifiers = repelem(animal_condition_identifier, input_params.repfactor);

% Get the frames with good tracking
good_frames = analysisstruct.frames_with_good_tracking{1, 1};

% Extract cluster IDs and corresponding identifiers for good frames
cluster_ids = analysisstruct.annot_reordered{end};
frame_identifiers = upsampled_identifiers(good_frames);

% Extract unique animal IDs and conditions
[animal_ids, ~, animal_indices] = unique(cellfun(@(x) x(1:end-2), frame_identifiers, 'UniformOutput', false));
conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);

%% Sequence and State Analysis
disp('Starting Sequence and State Analysis...');

% Parameters for sequence analysis
params = struct();
params.do_show_pdistmatrix = 1;   % Show distance matrix visualization
params.decimation_factor = 10;    % No additional downsampling
params.doclustering = 1;          % Perform clustering
params.corr_threshold = 0.2;      % Correlation threshold for sequence identification
params.clustercutoff = 0.65;      % Cutoff for hierarchical clustering
params.timescales = [1/4, 2];     % Timescales in minutes (15 seconds and 2 minutes)

% Separate cluster IDs for 'S' and 'F' conditions
cluster_ids_S = cluster_ids(strcmp(conditions, 'S'));
cluster_ids_F = cluster_ids(strcmp(conditions, 'F'));

% Create a 2x1 cell array for annot_reordered
analysisstruct_seq.annot_reordered = {cluster_ids_S; cluster_ids_F};

% Update condition_inds to reflect the new structure
analysisstruct_seq.condition_inds = [ones(size(cluster_ids_S)), 2 * ones(size(cluster_ids_F))];

% Update conditionnames if necessary
analysisstruct_seq.conditionnames = {'S', 'F'};

% Ensure other fields are correctly set
analysisstruct_seq.tsnegranularity = input_params.repfactor;
analysisstruct_seq.zValues = analysisstruct.zValues;
analysisstruct_seq.density_objects = max(cluster_ids);
analysisstruct_seq.mocapstruct_reduced_agg = {analysisstruct.mocapstruct_reduced_agg{1}};

% Set annotation_choose to use both conditions
annotation_choose = [1, 2];

% Run the sequence analysis
hierarchystruct = find_sequences_states_demo(analysisstruct_seq, annotation_choose, params);
save(fullfile(rootpath, 'sequence_hyerarchystruct.mat'), "hierarchystruct")
%% Visualize and analyze results
unique_conditions =unique(conditions, 'stable');
% Plot transition probabilities for each timescale
for i = 1:length(params.timescales)
    figure;
    imagesc(hierarchystruct.integrated_correlation_histograms{i});
    colorbar;
    title(sprintf('Transition Probabilities (Timescale: %.2f min)', params.timescales(i)));
    xlabel('To Condition');
    ylabel('From Condition');
    set(gca, 'XTick', 1:length(unique_conditions), 'XTickLabel', unique_conditions);
    set(gca, 'YTick', 1:length(unique_conditions), 'YTickLabel', unique_conditions);
end

% Analyze behavioral composition for each condition
for i = 1:length(params.timescales)
    figure;
    bar(hierarchystruct.behavior_composition_vector{i}');
    title(sprintf('Behavioral Composition (Timescale: %.2f min)', params.timescales(i)));
    xlabel('Sequence Cluster ID');
    ylabel('Proportion');
    legend(unique_conditions, 'Location', 'eastoutside');
end

% Analyze sequences
for i = 1:length(params.timescales)
    seq_clusters = hierarchystruct.clustered_behavior{i};
    conditions = hierarchystruct.base_annotation_inds;  % Assuming this contains condition labels

    % Get unique clusters and conditions
    unique_clusters = unique(seq_clusters);
    unique_conditions = unique(conditions);

    % Initialize transition matrices for each condition
    transitions = zeros(length(unique_clusters), length(unique_clusters), length(unique_conditions));

    % Count transitions for each condition
    for cond = 1:length(unique_conditions)
        cond_mask = conditions == unique_conditions(cond);
        cond_seq = seq_clusters(cond_mask);

        for j = 1:length(cond_seq)-1
            from = find(unique_clusters == cond_seq(j));
            to = find(unique_clusters == cond_seq(j+1));
            if ~isempty(from) && ~isempty(to)
                transitions(from, to, cond) = transitions(from, to, cond) + 1;
            end
        end
    end

    % Visualize transitions for each condition
    figure('Position', [100, 100, 1200, 500]);
    for cond = 1:length(unique_conditions)
        subplot(1, length(unique_conditions), cond);
        imagesc(transitions(:,:,cond));
        colorbar;
        title(sprintf('Transition Counts - Condition %s (Timescale: %.2f min)', ...
            char(unique_conditions(cond)), params.timescales(i)));
        xlabel('To Sequence Cluster');
        ylabel('From Sequence Cluster');
        set(gca, 'XTick', 1:length(unique_clusters), 'XTickLabel', unique_clusters);
        set(gca, 'YTick', 1:length(unique_clusters), 'YTickLabel', unique_clusters);
    end

    % Compute and visualize transition probability differences
    if length(unique_conditions) == 2  % Assuming two conditions (e.g., S and F)
        prob_diff = (transitions(:,:,2) ./ sum(transitions(:,:,2), 2)) - ...
            (transitions(:,:,1) ./ sum(transitions(:,:,1), 2));

        figure;
        imagesc(prob_diff);
        colorbar;
        title(sprintf('Transition Probability Difference (Cond2 - Cond1, Timescale: %.2f min)', ...
            params.timescales(i)));
        xlabel('To Sequence Cluster');
        ylabel('From Sequence Cluster');
        set(gca, 'XTick', 1:length(unique_clusters), 'XTickLabel', unique_clusters);
        set(gca, 'YTick', 1:length(unique_clusters), 'YTickLabel', unique_clusters);
    end

    % Analyze most common transitions
    [sorted_trans, sorted_idx] = sort(transitions(:), 'descend');
    top_n = 5;  % Number of top transitions to display

    fprintf('Top %d transitions for Timescale %.2f min:\n', top_n, params.timescales(i));
    for k = 1:top_n
        [from, to, cond] = ind2sub(size(transitions), sorted_idx(k));
        fprintf('Condition %s: %d -> %d (Count: %d)\n', ...
            num2str(unique_conditions(cond)), unique_clusters(from), unique_clusters(to), sorted_trans(k));
    end
    fprintf('\n');
end

% Compare sequence patterns between conditions
for i = 1:length(params.timescales)
    seq_clusters = hierarchystruct.clustered_behavior{i};
    conditions = hierarchystruct.base_annotation_inds;
    unique_clusters = unique(seq_clusters);

    % Separate sequences by condition
    seq_S = seq_clusters(conditions == 1);  % Assuming 1 is for 'S'
    seq_F = seq_clusters(conditions == 2);  % Assuming 2 is for 'F'

    % Compute transition probabilities for each condition
    trans_prob_S = compute_transition_probabilities(seq_S, unique_clusters);
    trans_prob_F = compute_transition_probabilities(seq_F, unique_clusters);

    % Visualize difference in transition probabilities
    figure('Position', [100, 100, 1200, 900]);

    % Subplot 1: Transition probabilities for S condition
    subplot(2, 2, 1);
    imagesc(trans_prob_S);
    colorbar;
    title(sprintf('Transition Probabilities - S Condition (Timescale: %.2f min)', params.timescales(i)));
    xlabel('To Sequence Cluster');
    ylabel('From Sequence Cluster');
    set(gca, 'XTick', 1:length(unique_clusters), 'XTickLabel', unique_clusters);
    set(gca, 'YTick', 1:length(unique_clusters), 'YTickLabel', unique_clusters);

    % Subplot 2: Transition probabilities for F condition
    subplot(2, 2, 2);
    imagesc(trans_prob_F);
    colorbar;
    title(sprintf('Transition Probabilities - F Condition (Timescale: %.2f min)', params.timescales(i)));
    xlabel('To Sequence Cluster');
    ylabel('From Sequence Cluster');
    set(gca, 'XTick', 1:length(unique_clusters), 'XTickLabel', unique_clusters);
    set(gca, 'YTick', 1:length(unique_clusters), 'YTickLabel', unique_clusters);

    % Subplot 3: Difference in transition probabilities
    subplot(2, 2, 3);
    imagesc(trans_prob_F - trans_prob_S);
    colorbar;
    title(sprintf('Difference in Transition Probabilities (F - S) (Timescale: %.2f min)', params.timescales(i)));
    xlabel('To Sequence Cluster');
    ylabel('From Sequence Cluster');
    set(gca, 'XTick', 1:length(unique_clusters), 'XTickLabel', unique_clusters);
    set(gca, 'YTick', 1:length(unique_clusters), 'YTickLabel', unique_clusters);

    % Subplot 4: Significant differences
    subplot(2, 2, 4);
    sig_diff = (trans_prob_F - trans_prob_S) .* (abs(trans_prob_F - trans_prob_S) > 0.15);  % Threshold can be adjusted
    imagesc(sig_diff);
    colorbar;
    title(sprintf('Significant Differences in Transition Probabilities (Timescale: %.2f min)', params.timescales(i)));
    xlabel('To Sequence Cluster');
    ylabel('From Sequence Cluster');
    set(gca, 'XTick', 1:length(unique_clusters), 'XTickLabel', unique_clusters);
    set(gca, 'YTick', 1:length(unique_clusters), 'YTickLabel', unique_clusters);

    sgtitle(sprintf('Transition Probability Analysis (Timescale: %.2f min)', params.timescales(i)));
end

% Helper function to compute transition probabilities
function trans_prob = compute_transition_probabilities(sequence, unique_clusters)
    n_clusters = length(unique_clusters);
    trans_count = zeros(n_clusters);
    for j = 1:length(sequence)-1
        from = find(unique_clusters == sequence(j));
        to = find(unique_clusters == sequence(j+1));
        if ~isempty(from) && ~isempty(to)
            trans_count(from, to) = trans_count(from, to) + 1;
        end
    end
    row_sums = sum(trans_count, 2);
    trans_prob = trans_count ./ row_sums;
    trans_prob(isnan(trans_prob)) = 0;  % Handle division by zero
end

%% Plot graph motifs
% Choose a timescale (you can modify this to loop through all timescales)
timescale_index = 1; % For the short-term timescale

seq_clusters = hierarchystruct.clustered_behavior{timescale_index};
conditions = hierarchystruct.base_annotation_inds;

unique_clusters = unique(seq_clusters);
unique_conditions = unique(conditions);

% Initialize transition count matrix
trans_count = zeros(length(unique_clusters), length(unique_clusters), length(unique_conditions));

for cond = 1:length(unique_conditions)
    cond_mask = conditions == unique_conditions(cond);
    cond_seq = seq_clusters(cond_mask);

    for j = 1:length(cond_seq)-1
        from = find(unique_clusters == cond_seq(j), 1);
        to = find(unique_clusters == cond_seq(j+1), 1);
        if ~isempty(from) && ~isempty(to)
            trans_count(from, to, cond) = trans_count(from, to, cond) + 1;
        end
    end
end

% Compute transition probabilities
trans_prob = zeros(size(trans_count));
for cond = 1:length(unique_conditions)
    row_sums = sum(trans_count(:,:,cond), 2);
    for row = 1:size(trans_count, 1)
        if row_sums(row) > 0
            trans_prob(row,:,cond) = trans_count(row,:,cond) / row_sums(row);
        else
            trans_prob(row,:,cond) = 0;  % Set to 0 for rows with no transitions
        end
    end
end

% Compute difference in transition probabilities
prob_diff = trans_prob(:,:,2) - trans_prob(:,:,1); % Assuming F is 2 and S is 1

% Create a graph object
G = digraph;

% Add nodes and edges
threshold = 0.05; % Adjust this to show more or fewer transitions
for i = 1:length(unique_clusters)
    for j = 1:length(unique_clusters)
        if abs(prob_diff(i,j)) > threshold
            G = addedge(G, i, j, abs(prob_diff(i,j)));
        end
    end
end

% Plot the graph
figure('Position', [100, 100, 1000, 800]);
p = plot(G, 'Layout', 'force');

% Customize node appearance
p.NodeColor = [0.7 0.7 0.7];
% p.NodeLabel = G.;
p.MarkerSize = 1 + 5*sqrt(outdegree(G)); % Node size based on outgoing connections

% Customize edge appearance
p.EdgeColor = [0.5 0.5 0.5];
p.LineWidth = 2*G.Edges.Weight/max(G.Edges.Weight);

% Color edges based on which condition has higher probability
edge_colors = zeros(G.numedges, 3);
for e = 1:G.numedges
    if prob_diff(G.Edges.EndNodes(e,1), G.Edges.EndNodes(e,2)) > 0
        edge_colors(e,:) = [1 0 0]; % Red for higher in F condition
    else
        edge_colors(e,:) = [0 0 1]; % Blue for higher in S condition
    end
end
p.EdgeColor = edge_colors;

% Add a colorbar legend
colormap([0 0 1; 1 0 0]);
c = colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'S', 'F'});
c.Label.String = 'Dominant Condition';

% Set title and adjust layout
title(sprintf('Motif Transition Network (Timescale: %.2f min)', params.timescales(timescale_index)));
layout(p, 'force', 'WeightEffect', 'inverse');

% Add annotations
annotation('textbox', [0.15, 0.95, 0.7, 0.05], 'String', ...
    'Node size represents cluster prevalence. Edge thickness represents transition strength.', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none');

%% Save results
save(fullfile(rootpath, 'sequence_analysis_results.mat'), 'hierarchystruct', 'params');

disp('Sequence and State Analysis Complete!');

% % Helper function to compute transition probabilities
% function trans_prob = compute_transition_probabilities(sequence)
% max_cluster = max(sequence);
% trans_count = zeros(max_cluster);
% for j = 1:length(sequence)-1
%     from = sequence(j);
%     to = sequence(j+1);
%     if from > 0 && to > 0
%         trans_count(from, to) = trans_count(from, to) + 1;
%     end
% end
% trans_prob = trans_count ./ sum(trans_count, 2);
% trans_prob(isnan(trans_prob)) = 0;
% end