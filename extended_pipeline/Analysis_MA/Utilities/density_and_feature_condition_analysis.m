
clc;
logger('Starting pose individuality analysis script', 'INFO');
clear;
close all;

GC = general_configs;
rootpath = GC.preprocessing_rootpath;

%% Define experimental groups and their conditions
% Formalin experiment group
formalin_experiment = struct();
formalin_experiment.name = 'Formalin Experiment';
formalin_experiment.conditions = {'S', 'F', 'C'};
formalin_experiment.condition_names = {'Saline', 'Formalin', 'Formalin+Carprofen'};
formalin_experiment.comparisons = {{'S', 'F'}, {'F', 'C'}}; % baseline vs treatment, treatment vs intervention

% Neuropathic experiment group
neuropathic_experiment = struct();
neuropathic_experiment.name = 'Neuropathic Experiment';
neuropathic_experiment.conditions = {'H', 'N', 'G'};
neuropathic_experiment.condition_names = {'Control', 'Neuropathic', 'Gabapentin'};
neuropathic_experiment.comparisons = {{'H', 'N'}, {'N', 'G'}}; % control vs pathology, pathology vs treatment

% Store all experiments
experiments = {formalin_experiment, neuropathic_experiment};

%% User Input Dialog for Experiment and Analysis Mode Selection
% Create experiment selection dialog
experiment_options = cellfun(@(x) x.name, experiments, 'UniformOutput', false);
[experiment_idx, ok] = listdlg('PromptString', 'Select experiment type:', ...
                              'SelectionMode', 'single', ...
                              'ListString', experiment_options, ...
                              'Name', 'Experiment Selection', ...
                              'ListSize', [300, 150]);

if ~ok
    logger('Analysis cancelled by user', 'WARNING');
    return;
end

selected_experiment = experiments{experiment_idx};
logger(sprintf('Selected experiment: %s', selected_experiment.name), 'INFO');

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
    logger(sprintf('Running in DEBUG MODE for %s', selected_experiment.name), 'INFO');
else
    debugging = false;
    visualize = 'off';  % Don't show figures in production mode
    do_export = true;   % Export figures in production mode
    logger(sprintf('Running in PRODUCTION MODE for %s', selected_experiment.name), 'INFO');
end

% Export folder
export_folder = fullfile(GC.temp_root, 'figs_paper');
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

% Extract condition identifiers and animal names
conditions = cellfun(@(x) x(end), animal_list_used_after_analysis, 'UniformOutput', false);
animal_names = cellfun(@(x) x(1:find(x=='_',1)-1), animal_list_used_after_analysis, 'UniformOutput', false);

% Check which conditions from the selected experiment are available in the data
available_experiment_conditions = intersect(selected_experiment.conditions, unique(conditions));
logger(sprintf('Available conditions for %s: %s', selected_experiment.name, strjoin(available_experiment_conditions, ', ')), 'INFO');

if length(available_experiment_conditions) < 2
    logger(sprintf('Insufficient conditions available for %s. Need at least 2 conditions.', selected_experiment.name), 'ERROR');
    return;
end

% Filter comparisons to only include those with available conditions
valid_comparisons = {};
for i = 1:length(selected_experiment.comparisons)
    comparison = selected_experiment.comparisons{i};
    if all(ismember(comparison, available_experiment_conditions))
        valid_comparisons{end+1} = comparison;
    else
        logger(sprintf('Skipping comparison %s vs %s - conditions not available in data', comparison{1}, comparison{2}), 'WARNING');
    end
end

if isempty(valid_comparisons)
    logger('No valid comparisons available for the selected experiment', 'ERROR');
    return;
end

logger(sprintf('Valid comparisons for %s: %d', selected_experiment.name, length(valid_comparisons)), 'INFO');

%% Calculate normalized difference maps for all valid comparisons
for comp_idx = 1:length(valid_comparisons)
    comparison = valid_comparisons{comp_idx};
    condition1 = comparison{1};
    condition2 = comparison{2};
    
    % Get condition names for display
    cond1_idx = find(strcmp(selected_experiment.conditions, condition1));
    cond2_idx = find(strcmp(selected_experiment.conditions, condition2));
    condition1_name = selected_experiment.condition_names{cond1_idx};
    condition2_name = selected_experiment.condition_names{cond2_idx};
    
    logger(sprintf('Generating normalized difference maps: %s (%s) vs %s (%s)', ...
        condition2_name, condition2, condition1_name, condition1), 'INFO');
    
    % Identify frames for each condition
    condition1_mask = strcmp(conditions, condition1);
    condition2_mask = strcmp(conditions, condition2);
    
    logger(sprintf('Found %d frames for %s condition', sum(condition1_mask), condition1_name), 'INFO');
    logger(sprintf('Found %d frames for %s condition', sum(condition2_mask), condition2_name), 'INFO');
    
    % Extract animal names for each condition
    condition1_animal_names = animal_names(condition1_mask);
    condition2_animal_names = animal_names(condition2_mask);
    
    % Extract z-values for each condition
    condition1_zvals = analysisstruct.zValues(condition1_mask, :);
    condition2_zvals = analysisstruct.zValues(condition2_mask, :);
    
    % Get unique animals in both conditions
    condition1_unique_animals = unique(condition1_animal_names);
    condition2_unique_animals = unique(condition2_animal_names);
    
    % Find common animals between both conditions
    common_animals = intersect(condition1_unique_animals, condition2_unique_animals);
    num_common_animals = length(common_animals);
    
    if num_common_animals == 0
        logger(sprintf('No common animals found between %s and %s conditions', condition1_name, condition2_name), 'WARNING');
        continue;
    end
    
    logger(sprintf('Found %d common animals between %s and %s conditions', num_common_animals, condition1_name, condition2_name), 'INFO');
    
    % Calculate density map difference for this comparison
    calculate_density_difference_maps(condition1, condition2, condition1_name, condition2_name, ...
        condition1_mask, condition2_mask, condition1_animal_names, condition2_animal_names, ...
        condition1_zvals, condition2_zvals, common_animals, analysisstruct, ...
        visualize, do_export, export_folder, selected_experiment.name);
end


%% Feature Similarity Analysis for all valid comparisons
logger('Starting feature similarity analysis for selected experiment', 'INFO');

% Combine joint features and extra joint features
feature_matrix = [analysisstruct.jt_features, analysisstruct.extra_jt_features];

% Get cluster assignments for fine behaviors if available
if isfield(analysisstruct, 'cluster_assignments')
    cluster_assignments = analysisstruct.cluster_assignments;
else
    % If no cluster assignments, create dummy assignment (single cluster)
    logger('No cluster assignments found, treating all data as a single behavior', 'WARNING');
    cluster_assignments = ones(size(feature_matrix, 1), 1);
end

% Perform feature similarity analysis for each valid comparison
for comp_idx = 1:length(valid_comparisons)
    comparison = valid_comparisons{comp_idx};
    condition1 = comparison{1};
    condition2 = comparison{2};
    
    % Get condition names for display
    cond1_idx = find(strcmp(selected_experiment.conditions, condition1));
    cond2_idx = find(strcmp(selected_experiment.conditions, condition2));
    condition1_name = selected_experiment.condition_names{cond1_idx};
    condition2_name = selected_experiment.condition_names{cond2_idx};
    
    logger(sprintf('Analyzing feature similarity: %s (%s) vs %s (%s)', ...
        condition1_name, condition1, condition2_name, condition2), 'INFO');
    
    % Perform the feature similarity analysis
    perform_feature_similarity_analysis(condition1, condition2, condition1_name, condition2_name, ...
        conditions, animal_names, feature_matrix, cluster_assignments, ...
        visualize, do_export, export_folder, selected_experiment.name);
end

logger('Feature similarity analysis completed', 'INFO');

%% Helper Functions

function calculate_density_difference_maps(condition1, condition2, condition1_name, condition2_name, ...
    condition1_mask, condition2_mask, condition1_animal_names, condition2_animal_names, ...
    condition1_zvals, condition2_zvals, common_animals, analysisstruct, ...
    visualize, do_export, export_folder, experiment_name)
    
    % Use global t-SNE parameters for consistent density estimation
    global_zvals = analysisstruct.zValues;
    
    num_common_animals = length(common_animals);
    
    % Calculate subplot layout for the difference maps
    n_cols = ceil(sqrt(num_common_animals));
    n_rows = ceil(num_common_animals / n_cols);
    
    % Create figure for the difference maps
    fig_diff = figure('Name', sprintf('Normalized Difference Maps: (%s - %s) / (%s + %s)', ...
        condition2_name, condition1_name, condition1_name, condition2_name), 'Color', 'w', ...
        'Position', [300, 300, 300*n_cols, 300*n_rows], 'Visible', visualize);
    
    % Initialize colormap for difference visualization
    % Red-White-Blue colormap for difference visualization (-1 to 1)
    diff_cmap = [flipud(hot(128)); winter(128)];
    diff_cmap = diff_cmap(32:224, :); % Trim extreme colors for better visualization
    
    % Loop through common animals and create difference maps
    for i = 1:num_common_animals
        animal = common_animals{i};
        logger(sprintf('Processing difference map for animal: %s', animal), 'INFO');
        
        % Get indices for this animal in each condition
        condition1_idx = strcmp(condition1_animal_names, animal);
        condition2_idx = strcmp(condition2_animal_names, animal);
        
        % Skip if not enough data points in either condition
        if sum(condition1_idx) < 10 || sum(condition2_idx) < 10
            logger(sprintf('Skipping animal %s - insufficient data points (%s:%d, %s:%d)', ...
                animal, condition1, sum(condition1_idx), condition2, sum(condition2_idx)), 'WARNING');
            continue;
        end
        
        % Get z-values for this animal in each condition
        animal_condition1_zvals = condition1_zvals(condition1_idx, :);
        animal_condition2_zvals = condition2_zvals(condition2_idx, :);
        
        % Create subplot for this animal
        subplot(n_rows, n_cols, i);
        h_ax = gca;
        set(h_ax, 'Color', 'w');
        
        % Generate density maps for both conditions
        density_width = analysisstruct.params.density_width;
        density_max = max(global_zvals(:))*analysisstruct.params.expansion_factor;
        density_res = analysisstruct.params.density_res;
        
        % Generate density maps using plotdensitymaps function
        [condition1_density_maps, xx, yy, ~, ~] = plotdensitymaps({animal_condition1_zvals}, 1, h_ax, ...
            density_width, density_max, density_res, 'Blues9');
        
        cla(h_ax);
        
        [condition2_density_maps, ~, ~, ~, ~] = plotdensitymaps({animal_condition2_zvals}, 1, h_ax, ...
            density_width, density_max, density_res, 'Reds9');
        
        cla(h_ax);
        
        % Extract the actual density maps
        condition1_map = condition1_density_maps{1};
        condition2_map = condition2_density_maps{1};
        
        % Handle potential NaN or Inf values
        condition1_map(isnan(condition1_map) | isinf(condition1_map)) = 0;
        condition2_map(isnan(condition2_map) | isinf(condition2_map)) = 0;
        
        % Apply small epsilon to avoid division by zero
        epsilon = 1e-10;
        
        % Calculate normalized difference map: (condition2-condition1)/(condition1+condition2)
        diff_map = (condition2_map - condition1_map) ./ (condition1_map + condition2_map + epsilon);
        
        % Clip extreme values for better visualization
        diff_map = max(min(diff_map, 1), -1);
        
        % Display the difference map
        imagesc(diff_map);
        colormap(h_ax, diff_cmap);
        caxis([-1, 1]); % Set color scale to -1 to 1
        
        % Add colorbar
        c = colorbar;
        c.Label.String = sprintf('(%s-%s)/(%s+%s)', condition2, condition1, condition1, condition2);
        c.Label.FontSize = 12;
        c.FontSize = 10;
        
        % Add title with animal name
        title(['Animal: ' animal], 'FontSize', 12);
        
        % Turn off axis
        axis off;
        
        logger(sprintf('Completed difference map for animal: %s', animal), 'INFO');
    end
    
    % Add overall title
    sgtitle(sprintf('Normalized Difference Maps: (%s - %s) / (%s + %s)', ...
        condition2_name, condition1_name, condition1_name, condition2_name), 'FontSize', 14);
    
    % Export the difference maps figure if needed
    if do_export
        logger('Exporting normalized difference maps to export folder', 'INFO');
        export_name = sprintf('normalized_difference_maps_%s_%s_vs_%s', ...
            strrep(experiment_name, ' ', '_'), condition2, condition1);
        exportgraphics(fig_diff, fullfile(export_folder, [export_name '.pdf']), ...
            'ContentType', 'vector', 'BackgroundColor', 'none');
    end
    
    % Create combined difference map (all animals averaged)
    create_combined_difference_map(condition1, condition2, condition1_name, condition2_name, ...
        condition1_animal_names, condition2_animal_names, condition1_zvals, condition2_zvals, ...
        common_animals, analysisstruct, visualize, do_export, export_folder, experiment_name);
end

function create_combined_difference_map(condition1, condition2, condition1_name, condition2_name, ...
    condition1_animal_names, condition2_animal_names, condition1_zvals, condition2_zvals, ...
    common_animals, analysisstruct, visualize, do_export, export_folder, experiment_name)
    
    logger('Creating combined difference map (average of all animals)', 'INFO');
    
    % Figure for combined difference map
    fig_combined_diff = figure('Name', sprintf('Combined Normalized Difference Map (%s vs %s)', ...
        condition2_name, condition1_name), 'Color', 'w', ...
        'Position', [400, 400, 500, 500], 'Visible', visualize);
    
    % Get density parameters
    density_width = analysisstruct.params.density_width;
    global_zvals = analysisstruct.zValues;
    density_max = max(global_zvals(:))*analysisstruct.params.expansion_factor;
    density_res = analysisstruct.params.density_res;
    
    % Initialize arrays to store accumulated density maps
    accumulated_condition1_map = zeros(density_res, density_res);
    accumulated_condition2_map = zeros(density_res, density_res);
    animal_count = 0;
    
    % Loop through common animals to accumulate density maps
    for i = 1:length(common_animals)
        animal = common_animals{i};
        
        % Get indices for this animal in each condition
        condition1_idx = strcmp(condition1_animal_names, animal);
        condition2_idx = strcmp(condition2_animal_names, animal);
        
        % Skip if not enough data points in either condition
        if sum(condition1_idx) < 10 || sum(condition2_idx) < 10
            continue;
        end
        
        % Get z-values for this animal in each condition
        animal_condition1_zvals = condition1_zvals(condition1_idx, :);
        animal_condition2_zvals = condition2_zvals(condition2_idx, :);
        
        % Temporary invisible axes for generating density maps
        temp_ax = axes('Visible', 'off');
        
        % Generate density maps
        [condition1_density_maps, xx, yy, ~, ~] = plotdensitymaps({animal_condition1_zvals}, 1, temp_ax, ...
            density_width, density_max, density_res);
        
        cla(temp_ax);
        
        [condition2_density_maps, ~, ~, ~, ~] = plotdensitymaps({animal_condition2_zvals}, 1, temp_ax, ...
            density_width, density_max, density_res);
        
        % Delete temporary axes
        delete(temp_ax);
        
        % Extract density maps
        condition1_map = condition1_density_maps{1};
        condition2_map = condition2_density_maps{1};
        
        % Handle potential NaN or Inf values
        condition1_map(isnan(condition1_map) | isinf(condition1_map)) = 0;
        condition2_map(isnan(condition2_map) | isinf(condition2_map)) = 0;
        
        % Accumulate density maps
        accumulated_condition1_map = accumulated_condition1_map + condition1_map;
        accumulated_condition2_map = accumulated_condition2_map + condition2_map;
        animal_count = animal_count + 1;
        
        logger(sprintf('Added animal %s to combined map (%d of %d)', animal, animal_count, length(common_animals)), 'INFO');
    end
    
    % Calculate average density maps and combined difference map
    if animal_count > 0
        avg_condition1_map = accumulated_condition1_map;
        avg_condition2_map = accumulated_condition2_map;
        
        % Calculate combined difference map
        epsilon = 1e-10;
        combined_diff_map = (avg_condition2_map - avg_condition1_map) ./ (avg_condition1_map + avg_condition2_map + epsilon);
        
        % Display the combined difference map
        h_ax = axes('Parent', fig_combined_diff);
        imagesc(flipud(combined_diff_map));
        caxis([-1, 1]);
        
        % Add colorbar
        c = colorbar;
        c.Label.String = sprintf('(%s-%s)/(%s+%s)', condition2, condition1, condition1, condition2);
        c.Label.FontSize = 14;
        c.FontSize = 12;
        
        % Add title
        title(sprintf('Combined Normalized Difference Map: (%s - %s) / (%s + %s)', ...
            condition2_name, condition1_name, condition1_name, condition2_name), 'FontSize', 14);
        subtitle(sprintf('Average of %d animals', animal_count), 'FontSize', 12);
        
        % Turn off axis
        axis off;
        
        % Export the combined difference map if needed
        if do_export
            logger('Exporting combined normalized difference map to export folder', 'INFO');
            export_name = sprintf('combined_normalized_difference_map_%s_%s_vs_%s', ...
                strrep(experiment_name, ' ', '_'), condition2, condition1);
            exportgraphics(fig_combined_diff, fullfile(export_folder, [export_name '.pdf']), ...
                'ContentType', 'vector', 'BackgroundColor', 'none');
        end
    else
        logger('No valid animals to create combined difference map', 'ERROR');
    end
end

function perform_feature_similarity_analysis(condition1, condition2, condition1_name, condition2_name, ...
    conditions, animal_names, feature_matrix, cluster_assignments, ...
    visualize, do_export, export_folder, experiment_name)
    
    logger(sprintf('Analyzing within-condition cosine similarity: %s vs %s', condition1_name, condition2_name), 'INFO');
    
    % Extract animal base names (without condition suffix)
    animal_base_names = cellfun(@(x) x(1:find(x=='_',1)-1), animal_names, 'UniformOutput', false);
    
    % Find indices for each condition
    condition1_mask = strcmp(conditions, condition1);
    condition2_mask = strcmp(conditions, condition2);
    
    condition1_animal_names = animal_base_names(condition1_mask);
    condition2_animal_names = animal_base_names(condition2_mask);
    
    % Get unique animals in each condition
    condition1_unique_animals = unique(condition1_animal_names);
    condition2_unique_animals = unique(condition2_animal_names);
    
    num_condition1_animals = length(condition1_unique_animals);
    num_condition2_animals = length(condition2_unique_animals);
    
    logger(sprintf('Found %d animals in condition %s and %d animals in condition %s', ...
        num_condition1_animals, condition1, num_condition2_animals, condition2), 'INFO');
    
    unique_clusters = unique(cluster_assignments);
    num_clusters = length(unique_clusters);
    logger(sprintf('Found %d unique finely parsed behaviors (clusters)', num_clusters), 'INFO');
    
    % Initialize arrays to store mean cosine similarity for each animal in each condition
    condition1_animal_similarities = zeros(num_condition1_animals, 1);
    condition2_animal_similarities = zeros(num_condition2_animals, 1);
    
    % Process each animal in condition 1
    for i = 1:num_condition1_animals
        animal = condition1_unique_animals{i};
        logger(sprintf('Processing animal %s in condition %s', animal, condition1), 'INFO');
        
        % Find indices for this animal in this condition
        animal_indices = find(condition1_mask & strcmp(animal_base_names, animal));
        
        if length(animal_indices) < 2
            logger(sprintf('Insufficient data for animal %s in condition %s, skipping', ...
                animal, condition1), 'WARNING');
            condition1_animal_similarities(i) = NaN;
            continue;
        end
        
        condition1_animal_similarities(i) = calculate_animal_similarity(animal_indices, feature_matrix, ...
            cluster_assignments, unique_clusters, condition1, animal);
    end
    
    % Process each animal in condition 2
    for i = 1:num_condition2_animals
        animal = condition2_unique_animals{i};
        logger(sprintf('Processing animal %s in condition %s', animal, condition2), 'INFO');
        
        % Find indices for this animal in this condition
        animal_indices = find(condition2_mask & strcmp(animal_base_names, animal));
        
        if length(animal_indices) < 2
            logger(sprintf('Insufficient data for animal %s in condition %s, skipping', ...
                animal, condition2), 'WARNING');
            condition2_animal_similarities(i) = NaN;
            continue;
        end
        
        condition2_animal_similarities(i) = calculate_animal_similarity(animal_indices, feature_matrix, ...
            cluster_assignments, unique_clusters, condition2, animal);
    end
    
    % Remove NaN values and create visualization
    create_similarity_visualization(condition1_animal_similarities, condition2_animal_similarities, ...
        condition1, condition2, condition1_name, condition2_name, ...
        visualize, do_export, export_folder, experiment_name);
end

function animal_similarity = calculate_animal_similarity(animal_indices, feature_matrix, ...
    cluster_assignments, unique_clusters, condition, animal)
    
    num_clusters = length(unique_clusters);
    behavior_similarities = zeros(num_clusters, 1);
    valid_behavior_count = 0;
    
    % Process each behavior cluster
    for cluster_idx = 1:num_clusters
        cluster = unique_clusters(cluster_idx);
        
        % Find frames for this cluster in this condition for this animal
        cluster_indices = animal_indices(cluster_assignments(animal_indices) == cluster);
        
        if length(cluster_indices) < 2
            continue;
        end
        
        % Process consecutive time points (average over max 5 consecutive frames)
        processed_features = process_consecutive_frames(cluster_indices, feature_matrix, 5);
        
        if size(processed_features, 1) < 2
            continue;
        end
        
        % Calculate pairwise cosine similarity within this cluster's feature vectors
        cosine_distances = pdist(processed_features, 'cosine');
        cosine_similarities = 1 - cosine_distances;
        
        % Average the similarities for this behavior cluster
        avg_similarity = mean(cosine_similarities);
        
        % Store the result if valid
        if ~isnan(avg_similarity)
            behavior_similarities(cluster_idx) = avg_similarity;
            valid_behavior_count = valid_behavior_count + 1;
            
            logger(sprintf('Condition %s, Animal %s, Cluster %d: Average cosine similarity = %.4f (%d comparisons)', ...
                condition, animal, cluster, avg_similarity, length(cosine_similarities)), 'INFO');
        end
    end
    
    % Calculate the mean similarity across all behaviors for this animal
    if valid_behavior_count > 0
        animal_similarity = mean(behavior_similarities(behavior_similarities ~= 0));
        logger(sprintf('Condition %s, Animal %s: Mean cosine similarity = %.4f across %d behaviors', ...
            condition, animal, animal_similarity, valid_behavior_count), 'INFO');
    else
        animal_similarity = NaN;
        logger(sprintf('Animal %s in condition %s: No valid behaviors for comparison', animal, condition), 'WARNING');
    end
end

function create_similarity_visualization(condition1_animal_similarities, condition2_animal_similarities, ...
    condition1, condition2, condition1_name, condition2_name, ...
    visualize, do_export, export_folder, experiment_name)
    
    % Remove NaN values
    condition1_animal_similarities = condition1_animal_similarities(~isnan(condition1_animal_similarities));
    condition2_animal_similarities = condition2_animal_similarities(~isnan(condition2_animal_similarities));
    
    valid_condition1_count = length(condition1_animal_similarities);
    valid_condition2_count = length(condition2_animal_similarities);
    
    if valid_condition1_count == 0 || valid_condition2_count == 0
        logger('Insufficient valid data for one or both conditions, cannot proceed with comparison', 'ERROR');
        return;
    end
    
    % Calculate statistics for each condition
    condition1_mean = mean(condition1_animal_similarities);
    condition1_sem = std(condition1_animal_similarities) / sqrt(valid_condition1_count);
    
    condition2_mean = mean(condition2_animal_similarities);
    condition2_sem = std(condition2_animal_similarities) / sqrt(valid_condition2_count);
    
    % Perform statistical test (t-test)
    [h, p] = ttest2(condition1_animal_similarities, condition2_animal_similarities);
    significance_label = 'NS'; % Default: Not Significant
    if p < 0.05
        significance_label = '*';
        if p < 0.01
            significance_label = '**';
            if p < 0.001
                significance_label = '***';
            end
        end
    end
    
    logger(sprintf('Statistical comparison: p = %.4f (%s)', p, significance_label), 'INFO');
    
    % Create a bar plot
    fig_cosine = figure('Name', sprintf('Within-condition Cosine Similarity: %s vs %s', condition1_name, condition2_name), ...
        'Color', 'w', 'Position', [500, 500, 400, 400], 'Visible', visualize);
    
    % Define condition names for the plot
    condition_names = {condition1_name, condition2_name};
    
    % Define bar positions
    bar_positions = 1:2;
    bar_width = 0.7;
    
    % Plot bars
    hold on;
    b1 = bar(bar_positions(1), condition1_mean, bar_width, 'FaceColor', [0.8, 0.2, 0.2]);
    b2 = bar(bar_positions(2), condition2_mean, bar_width, 'FaceColor', [0.5, 0.5, 0.5]);
    
    % Add individual data points (scatter)
    scatter(repmat(bar_positions(1), valid_condition1_count, 1) + (rand(valid_condition1_count, 1)-0.5)*0.2, ...
        condition1_animal_similarities, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    
    scatter(repmat(bar_positions(2), valid_condition2_count, 1) + (rand(valid_condition2_count, 1)-0.5)*0.2, ...
        condition2_animal_similarities, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
    
    % Add error bars (SEM)
    errorbar(bar_positions(1), condition1_mean, condition1_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    errorbar(bar_positions(2), condition2_mean, condition2_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    % Add significance bar
    y_max = max([condition1_mean + condition1_sem, condition2_mean + condition2_sem]) * 1.1;
    plot([bar_positions(1), bar_positions(2)], [y_max, y_max], 'k', 'LineWidth', 1.5);
    text(mean(bar_positions), y_max*1.05, significance_label, 'HorizontalAlignment', 'center', 'FontSize', 14);
    
    % Set axis properties
    ylabel('Cosine similarity between feature vectors', 'FontSize', 12);
    title('Within fine behaviors', 'FontSize', 14);
    xticks(bar_positions);
    xticklabels(condition_names);
    xlim([0.5, 2.5]);
    ylim([0, max([condition1_mean, condition2_mean])*1.3]);
    box off;
    
    % Create cleaner y-axis with fewer ticks
    yticks_range = get(gca, 'YTick');
    if length(yticks_range) > 5
        yticks(linspace(0, max(yticks_range), 5));
    end
    
    % Print statistics to console
    logger('Within-condition cosine similarity statistics:', 'INFO');
    logger(sprintf('Condition %s: Mean = %.4f, SEM = %.4f, n = %d', ...
        condition1, condition1_mean, condition1_sem, valid_condition1_count), 'INFO');
    logger(sprintf('Condition %s: Mean = %.4f, SEM = %.4f, n = %d', ...
        condition2, condition2_mean, condition2_sem, valid_condition2_count), 'INFO');
    logger(sprintf('Statistical comparison: p = %.4f (%s)', p, significance_label), 'INFO');
    
    % Export the figure if needed
    if do_export
        logger('Exporting within-condition cosine similarity figure to export folder', 'INFO');
        export_name = sprintf('within_condition_cosine_similarity_%s_%s_vs_%s', ...
            strrep(experiment_name, ' ', '_'), condition1, condition2);
        exportgraphics(fig_cosine, fullfile(export_folder, [export_name '.pdf']), ...
            'ContentType', 'vector', 'BackgroundColor', 'none');
    end
end


%% Helper function to process consecutive frames
function processed_features = process_consecutive_frames(indices, feature_matrix, max_consecutive)
    % This function averages feature vectors for consecutive frames up to max_consecutive
    
    % Sort indices to ensure proper sequence
    indices = sort(indices);
    
    % Initialize processed features list
    processed_features = [];
    
    % Process indices
    i = 1;
    while i <= length(indices)
        % Start a new consecutive sequence
        consec_start = i;
        
        % Find end of consecutive sequence (up to max_consecutive)
        while (i < length(indices)) && (i - consec_start < max_consecutive - 1) && ...
              (indices(i+1) == indices(i) + 1)
            i = i + 1;
        end
        
        % Extract and average the feature vectors in this consecutive sequence
        consec_indices = indices(consec_start:i);
        avg_feature = mean(feature_matrix(consec_indices, :), 1);
        
        % Add to processed features
        processed_features = [processed_features; avg_feature];
        
        % Move to next sequence
        i = i + 1;
    end
end
