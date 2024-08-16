%% Preamble
% This script reads labeled facial pain frames and determines which clusters 
% they belong to. The process is performed on a per-animal basis. The script 
% first creates a cluster vector that indicates the cluster to which each 
% frame belongs. It then iterates through the animals, matches labeled frames 
% to clusters, and calculates the proportion of pain and no-pain frames per cluster.

%% INIT

clear; close all; clc;
global GC

% Load cluster results
fps = GC.frame_rate;
rootpath = GC.preprocessing_rootpath;
clusters_struct_file = fullfile(rootpath, 'clusters_struct.mat');
load(fullfile(rootpath, 'cluster_analysis_results.mat'), 'results');

% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');

% Load cluster vectors for all the animals
load(clusters_struct_file, 'clusters_struct');

%% Iterate through animals
animal_ids = fieldnames(clusters_struct);
animal_ids(ismember(animal_ids, 'conditions')) = [];

predominantF = results.predominantF;
predominantS = results.predominantS;
n_predominantF = length(predominantF);
n_predominantS = length(predominantS);

% Create empty cell array to fill up the clusters per condition per animal
clusters_pain_no_pain = cell(length(animal_ids), 2); % 2 columns for pain and no pain proportions
prop_clusters_pain_no_pain = clusters_pain_no_pain;
concat_coi_pain = [];

for iid = 1:length(animal_ids)
    animal_ID = animal_ids{iid}; % Use the animal ID from the list

    % Calibrate condition
    if ~endsWith(animal_ID, 'F')
        continue
    end

    %% Match labeled frames to clusters
    % Load pain frames
    facial_rootpath = fullfile(rootpath, 'facial_labels');
    labels_file = fullfile(facial_rootpath, [animal_ID '_labels.csv']);
    T = readtable(labels_file);
    pain_frames = find(T.pain == 1);

    % Clusters of interest
    coi = clusters_struct.(animal_ID);
    coi_pain = coi(pain_frames);
    coi_pain = coi_pain(~isnan(coi_pain));
    coi_pain = unique(coi_pain, 'stable');
    coi_pain_inter = intersect(coi_pain, predominantF); % Pain clusters
    proportion_pain = 100 * numel(coi_pain_inter) / numel(coi_pain);

    % No pain clusters
    coi_no_pain = intersect(coi_pain, predominantS);
    proportion_no_pain = 100 * numel(coi_no_pain) / numel(coi_pain);

    % Store results
    prop_clusters_pain_no_pain{iid, 1} = proportion_pain;
    prop_clusters_pain_no_pain{iid, 2} = proportion_no_pain;
    clusters_pain_no_pain{iid, 1} = coi_pain_inter;
    clusters_pain_no_pain{iid, 2} = coi_no_pain;
    concat_coi_pain = cat(2, concat_coi_pain, coi_pain);
end

%% Do statistics and plots
% Convert cell array to matrix for plotting
proportions_matrix = cell2mat(prop_clusters_pain_no_pain);

% Generate bar plots
figure;
bar(proportions_matrix);
legend('Pain Proportion', 'No Pain Proportion');
title('Proportions of Pain and No Pain Clusters');
xlabel('Animal ID');
ylabel('Proportion (%)');

% Plot the facial expression clusters
clusters_to_plot = cell2mat(clusters_pain_no_pain(:, 1));
plot_identified_clusters(clusters_to_plot, analysisstruct);
plot_identified_clusters(unique(concat_coi_pain), analysisstruct);

%% Enhanced Bar Plot
figure('Position', [100, 100, 1200, 600]);
b = bar(proportions_matrix, 'grouped');
b(1).FaceColor = [0.8, 0.2, 0.2];  % Red for Pain
b(2).FaceColor = [0, 0.4, 0.8];    % Blue for No Pain

title('Proportions of Pain and No Pain Clusters per Animal', 'FontSize', 16);
xlabel('Animal ID', 'FontSize', 14);
ylabel('Proportion (%)', 'FontSize', 14);
legend('Pain Proportion', 'No Pain Proportion', 'Location', 'bestoutside');

% Add value labels on top of each bar
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1, ytips1, labels1, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2, ytips2, labels2, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Customize x-axis
xticks(1:length(animal_ids));
xticklabels(animal_ids);
xtickangle(45);

% Add a reference line at 100%
yline(100, '--k', '100% Reference');
grid off;

%% Venn Diagram
figure('Position', [100, 100, 800, 600]);
[~, ~] = venn([n_predominantF, length(unique(concat_coi_pain)), length(intersect(predominantF, unique(concat_coi_pain)))], 'ErrMinMode', 'None');
title('Overlap of Manual and Automatic Pain Clusters', 'FontSize', 16);
legend({'CAPTURE Pain Clusters', 'Grimace Pain Clusters'});

%% Cluster frequency and duration analysis
% Initialize arrays to store results
frequency_pain_clusters = zeros(length(animal_ids), 1);
mean_duration_pain_clusters = zeros(length(animal_ids), 1);
total_duration_pain_clusters = zeros(length(animal_ids), 1);

for i = 1:length(animal_ids)
    animal_ID = animal_ids{i};
    
    % Skip if not a formalin-treated animal
    if ~endsWith(animal_ID, 'F')
        continue;
    end
    
    % Get clusters for this animal
    animal_clusters = clusters_struct.(animal_ID);
    
    % Find pain-related clusters
    pain_cluster_indices = ismember(animal_clusters', predominantF);
    
    % Calculate frequency (number of pain cluster occurrences)
    frequency_pain_clusters(i) = sum(pain_cluster_indices) / length(pain_cluster_indices);
    
    % Find durations of continuous pain cluster segments
    diff_indices = diff([0; pain_cluster_indices; 0]);
    start_indices = find(diff_indices == 1);
    end_indices = find(diff_indices == -1) - 1;
    durations = end_indices - start_indices + 1;
    durations = durations / fps; % Convert to time
    
    % Calculate mean and total duration
    mean_duration_pain_clusters(i) = median(durations);
    total_duration_pain_clusters(i) = sum(durations);
end

% Remove zeros (non-formalin animals)
frequency_pain_clusters(frequency_pain_clusters == 0) = [];
mean_duration_pain_clusters(mean_duration_pain_clusters == 0) = [];
total_duration_pain_clusters(total_duration_pain_clusters == 0) = [];

% Visualize results
figure('Color', 'w');

% Frequency plot
subplot(3, 1, 1);
bar(frequency_pain_clusters);
title('Frequency of CAPTURE Pain Clusters per Animal');
xlabel('Animal Index');
ylabel('Frequency');
box off;

% Mean Duration plot
subplot(3, 1, 2);
bar(mean_duration_pain_clusters);
title('Mean Duration of CAPTURE Pain Clusters per Animal');
xlabel('Animal Index');
ylabel('Mean Duration (s)');
box off;

% Total Duration plot
subplot(3, 1, 3);
bar(total_duration_pain_clusters);
title('Total Duration of CAPTURE Pain Clusters per Animal');
xlabel('Animal Index');
ylabel('Total Duration (s)');
box off;

% Overall statistics
disp(['Mean Frequency: ' num2str(mean(frequency_pain_clusters))]);
disp(['Mean of Mean Durations: ' num2str(mean(mean_duration_pain_clusters))]);
disp(['Mean Total Duration: ' num2str(mean(total_duration_pain_clusters))]);