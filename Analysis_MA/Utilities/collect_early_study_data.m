%% Early Study Data Collection and Analysis
% This script processes and analyzes data from a cross-sectional study conducted in Oct 2023.
% The study includes miniscope recordings of animals that received either saline (S) or formalin (F) injections.
% The data was collected by Jun between X-Y 2023.
% 
% Author: [Your Name]
% Created: 30 Aug 2024
% Last Modified: [Current Date]

%% Initialize
clear;
close all;
global GC

%% Define study parameters
xs_animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
xs_conditions = {'F', 'F', 'S', 'S', 'N', 'S', 'F', 'F', 'F'};

%% Set paths
data_in_server_path = fullfile(GC.preprocessing_rootpath, 'miniscope_batch_behavior_data');

%% Load data
try
    load(fullfile(data_in_server_path, "analysisstruct_clusters.mat"), 'analysisstruct');
    load(fullfile(data_in_server_path, "agg_predictions.mat"), 'predictions', 'animal_condition_identifier');
catch ME
    error('Error loading data files: %s', ME.message);
end

%% Create animal_condition_identifier if not exists
if ~exist('animal_condition_identifier', 'var')
    animal_condition_identifier = cellfun(@(x) strcat(x, '_', xs_conditions{strcmp(x, xs_animal_list)}), ...
        animal_frames_identifier, 'UniformOutput', false);
    save(fullfile(data_in_server_path, "agg_predictions.mat"), 'predictions', 'animal_condition_identifier', '-append');
end

%% Prepare data for analysis
long_animal_frames_identifier = repelem(animal_condition_identifier, GC.repfactor);
animal_list_used_after_analysis = long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});
animal_list = unique(animal_list_used_after_analysis, 'stable');

% Create condition indices
cond_inds = zeros(1, length(analysisstruct.condition_inds));
for iid = 1:length(animal_list)
    animal_ID = animal_list{iid};
    if ~endsWith(animal_ID, '_N')
        idx = ismember(animal_list_used_after_analysis, animal_ID);
        cond_inds(idx) = iid;
    end
end

% Create condition indices for actual conditions (F or S)% conditions ie. actua conditions, concatenating animals in same cond
condition_inds = cond_inds; % sorting per condition
for iff = 1:length(animal_list_used_after_analysis)
     animal_ID = animal_list_used_after_analysis{iff};
     if endsWith(animal_ID, '_N')
         continue
     end
    condition_inds(iff) = double(endsWith(animal_ID, '_F') +1);
   
end

conds = unique(condition_inds);

%% Plot density maps
for icond = 2:3 % cond2: sal, cond3: F
    figure('Color', 'w', 'Name', sprintf('Condition %d', icond-1));
    this_cond = conds(icond);
    Zvals = analysisstruct.zValues(ismember(condition_inds, this_cond), :);
    plotdensitymaps({Zvals}, 1, gcf, analysisstruct.params.density_width, ...
        max(analysisstruct.zValues(:)) * analysisstruct.params.expansion_factor, ...
        analysisstruct.params.density_res);
end

%% Analyze cluster composition
cluster_ids = analysisstruct.annot_reordered{end};
frame_identifiers = long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});

unique_clusters = unique(cluster_ids);
num_clusters = length(unique_clusters);

cluster_density = zeros(num_clusters, 2);  % [F_count, S_count]

for c = 1:num_clusters
    cluster_index = unique_clusters(c);
    frames_in_cluster = cluster_ids == cluster_index;
    animal_conditions_in_cluster = frame_identifiers(frames_in_cluster);
    
    cluster_density(c, 1) = sum(cellfun(@(x) endsWith(x, '_F'), animal_conditions_in_cluster));
    cluster_density(c, 2) = sum(cellfun(@(x) endsWith(x, '_S'), animal_conditions_in_cluster));
end

%% Plot difference in frame counts between conditions
difference_frames = cluster_density(:, 1) - cluster_density(:, 2);

figure('Name', 'Difference in Frame Counts: Formalin vs. Saline');
bar(difference_frames);
title('Difference in Frame Counts: Formalin vs. Saline');
xlabel('Cluster');
ylabel('Difference in Frame Counts');
set(gca, 'XTick', 0:10:num_clusters, 'XTickLabel', arrayfun(@num2str, 0:10:num_clusters, 'UniformOutput', false));
box off;

%% Identify predominantly associated clusters
total_frames = sum(cluster_density, 2);
cluster_proportions = cluster_density ./ total_frames;

threshold = 0.99;
predominant_F = find(cluster_proportions(:, 1) >= threshold);
predominant_S = find(cluster_proportions(:, 2) >= threshold);

%% Plot predominant clusters
plot_predominant_clusters(analysisstruct, predominant_F, 'Predominant Formalin Clusters');

%% Helper function
function plot_predominant_clusters(analysisstruct, predominant_clusters, title_text)
    global GC
    fig = figure('Position', [10, 10, 2056, 1350], 'Color', 'w', 'Name', title_text);
    n_rows = ceil(sqrt(numel(predominant_clusters)));
    n_cols = ceil(sqrt(numel(predominant_clusters)));

    for ic = 1:numel(predominant_clusters)
        subplot(n_rows, n_cols, ic);
        this_cls = predominant_clusters(ic);
        plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1}, ...
            find(analysisstruct.annot_reordered{end} == this_cls), ...
            sprintf('Cluster %d', this_cls));
        title(sprintf('Cluster %d', this_cls));
    end

    sgtitle(title_text);
    
    % Save figure
    figure_name = fullfile(GC.figure_folder, [strrep(lower(title_text), ' ', '_'), '.pdf']);
    %export_fig(figure_name, '-pdf', fig);
end
