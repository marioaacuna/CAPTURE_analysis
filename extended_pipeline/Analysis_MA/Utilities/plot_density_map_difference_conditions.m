
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
% zvals = analysisstruct.zValues;

%% Calculate normalized difference maps between Saline and Formalin conditions
logger('Generating normalized difference maps between Saline and Formalin conditions', 'INFO');
% Use global t-SNE parameters for consistent density estimation
global_zvals = analysisstruct.zValues; % Global coordinates for reference

% First, identify frames for each condition
saline_mask = strcmp(conditions, 'S');
formalin_mask = strcmp(conditions, 'F');

logger(sprintf('Found %d frames for Saline condition', sum(saline_mask)), 'INFO');
logger(sprintf('Found %d frames for Formalin condition', sum(formalin_mask)), 'INFO');

% Extract animal names for each condition
saline_animal_names = animal_names(saline_mask);
formalin_animal_names = animal_names(formalin_mask);

% Extract z-values for each condition
saline_zvals = analysisstruct.zValues(saline_mask, :);
formalin_zvals = analysisstruct.zValues(formalin_mask, :);

% Get unique animals in both conditions
saline_unique_animals = unique(saline_animal_names);
formalin_unique_animals = unique(formalin_animal_names);

% Find common animals between both conditions
common_animals = intersect(saline_unique_animals, formalin_unique_animals);
num_common_animals = length(common_animals);

if num_common_animals == 0
    logger('No common animals found between Saline and Formalin conditions', 'ERROR');
    return;
end

logger(sprintf('Found %d common animals between Saline and Formalin conditions', num_common_animals), 'INFO');

% Calculate subplot layout for the difference maps
n_cols = ceil(sqrt(num_common_animals));
n_rows = ceil(num_common_animals / n_cols);

% Create figure for the difference maps
fig_diff = figure('Name', 'Normalized Difference Maps: (Saline - Formalin) / (Saline + Formalin)', 'Color', 'w', ...
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
    saline_idx = strcmp(saline_animal_names, animal);
    formalin_idx = strcmp(formalin_animal_names, animal);
    
    % Skip if not enough data points in either condition
    if sum(saline_idx) < 10 || sum(formalin_idx) < 10
        logger(sprintf('Skipping animal %s - insufficient data points (S:%d, F:%d)', ...
            animal, sum(saline_idx), sum(formalin_idx)), 'WARNING');
        continue;
    end
    
    % Get z-values for this animal in each condition
    animal_saline_zvals = saline_zvals(saline_idx, :);
    animal_formalin_zvals = formalin_zvals(formalin_idx, :);
    
    % Create subplot for this animal
    subplot(n_rows, n_cols, i);
    h_ax = gca;
    set(h_ax, 'Color', 'w');
    
    % Generate density maps for both conditions
    % Use the same settings for both to ensure comparable results
    density_width = analysisstruct.params.density_width;
    density_max = max(global_zvals(:))*analysisstruct.params.expansion_factor;
    density_res = analysisstruct.params.density_res;
    
    % Generate density maps using your plotdensitymaps function
    [saline_density_maps, xx, yy, ~, ~] = plotdensitymaps({animal_saline_zvals}, 1, h_ax, ...
        density_width, density_max, density_res, 'Blues9');
    
    % We need to clear the axis before plotting the next density map
    cla(h_ax);
    
    [formalin_density_maps, ~, ~, ~, ~] = plotdensitymaps({animal_formalin_zvals}, 1, h_ax, ...
        density_width, density_max, density_res, 'Reds9');
    
    % Clear the axis again for our difference map
    cla(h_ax);
    
    % Extract the actual density maps
    saline_map = saline_density_maps{1};
    formalin_map = formalin_density_maps{1};
    
    % Handle potential NaN or Inf values
    saline_map(isnan(saline_map) | isinf(saline_map)) = 0;
    formalin_map(isnan(formalin_map) | isinf(formalin_map)) = 0;
    
    % Apply small epsilon to avoid division by zero
    epsilon = 1e-10;
    
    % Calculate normalized difference map: (S-F)/(S+F)
    diff_map = (formalin_map - saline_map) ./ (saline_map + formalin_map + epsilon);
    
    % Clip extreme values for better visualization
    diff_map = max(min(diff_map, 1), -1);
    
    % Display the difference map
    imagesc((diff_map));
    colormap(h_ax, diff_cmap);
    caxis([-1, 1]); % Set color scale to -1 to 1
    
    % Add colorbar
    c = colorbar;
    c.Label.String = '(S-F)/(S+F)';
    c.Label.FontSize = 12;
    c.FontSize = 10;
    
    % Add title with animal name
    title(['Animal: ' animal], 'FontSize', 12);
    
    % Turn off axis
    axis off;
    
    logger(sprintf('Completed difference map for animal: %s', animal), 'INFO');
end

% Add overall title
sgtitle('Normalized Difference Maps: (Saline - Formalin) / (Saline + Formalin)', 'FontSize', 14);

% Export the difference maps figure if needed
if do_export
    logger('Exporting normalized difference maps to export folder', 'INFO');
    export_name = 'normalized_difference_maps_S_F';
    exportgraphics(fig_diff, fullfile(export_folder, [export_name '.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
end

%% Create combined difference map (all animals averaged)
logger('Creating combined difference map (average of all animals)', 'INFO');

% Figure for combined difference map
fig_combined_diff = figure('Name', 'Combined Normalized Difference Map (All Animals)', 'Color', 'w', ...
    'Position', [400, 400, 500, 500], 'Visible', visualize);

% Initialize arrays to store accumulated density maps
accumulated_saline_map = zeros(density_res, density_res);
accumulated_formalin_map = zeros(density_res, density_res);
animal_count = 0;

% Loop through common animals to accumulate density maps
for i = 1:num_common_animals
    animal = common_animals{i};
    
    % Get indices for this animal in each condition
    saline_idx = strcmp(saline_animal_names, animal);
    formalin_idx = strcmp(formalin_animal_names, animal);
    
    % Skip if not enough data points in either condition
    if sum(saline_idx) < 10 || sum(formalin_idx) < 10
        continue;
    end
    
    % Get z-values for this animal in each condition
    animal_saline_zvals = saline_zvals(saline_idx, :);
    animal_formalin_zvals = formalin_zvals(formalin_idx, :);
    
    % Temporary invisible axes for generating density maps
    temp_ax = axes('Visible', 'off');
    
    % Generate density maps
    [saline_density_maps, xx, yy, ~, ~] = plotdensitymaps({animal_saline_zvals}, 1, temp_ax, ...
        density_width, density_max, density_res);
    
    cla(temp_ax);
    
    [formalin_density_maps, ~, ~, ~, ~] = plotdensitymaps({animal_formalin_zvals}, 1, temp_ax, ...
        density_width, density_max, density_res);
    
    % Delete temporary axes
    delete(temp_ax);
    
    % Extract density maps
    saline_map = saline_density_maps{1};
    formalin_map = formalin_density_maps{1};
    
    % Handle potential NaN or Inf values
    saline_map(isnan(saline_map) | isinf(saline_map)) = 0;
    formalin_map(isnan(formalin_map) | isinf(formalin_map)) = 0;
    
    % Accumulate density maps
    accumulated_saline_map = accumulated_saline_map + saline_map;
    accumulated_formalin_map = accumulated_formalin_map + formalin_map;
    animal_count = animal_count + 1;
    
    logger(sprintf('Added animal %s to combined map (%d of %d)', animal, animal_count, num_common_animals), 'INFO');
end

% Calculate average density maps
if animal_count > 0
    avg_saline_map = accumulated_saline_map ;
    avg_formalin_map = accumulated_formalin_map ;
    
    % Calculate combined difference map
    epsilon = 1e-10;
    combined_diff_map = (avg_formalin_map - avg_saline_map) ./ (avg_saline_map + avg_formalin_map + epsilon);
    
    % Clip extreme values
    % combined_diff_map = max(min(combined_diff_map, 1), -1);
    
    % Display the combined difference map
    h_ax = axes('Parent', fig_combined_diff);
    imagesc(flipud(combined_diff_map));
    % colormap(h_ax, diff_cmap);
    caxis([-1, 1]);
    
    % Add colorbar
    c = colorbar;
    c.Label.String = '(S-F)/(S+F)';
    c.Label.FontSize = 14;
    c.FontSize = 12;
    
    % Add title
    title('Combined Normalized Difference Map: (Saline - Formalin) / (Saline + Formalin)', 'FontSize', 14);
    subtitle(sprintf('Average of %d animals', animal_count), 'FontSize', 12);
    
    % Turn off axis
    axis off;
    
    % Export the combined difference map if needed
    if do_export
        logger('Exporting combined normalized difference map to export folder', 'INFO');
        export_name = 'combined_normalized_difference_map_S_F';
        exportgraphics(fig_combined_diff, fullfile(export_folder, [export_name '.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
    end
else
    logger('No valid animals to create combined difference map', 'ERROR');
end

logger('Normalized difference map analysis completed', 'INFO');


%% 
%% Calculate cosine similarity between feature vectors for different conditions
logger('Analyzing cosine similarity between feature vectors across conditions', 'INFO');

% We'll analyze the cosine similarity for the feature matrix
% (combined joint features and extra joint features)
feature_matrix = [analysisstruct.jt_features, analysisstruct.extra_jt_features];

% Get condition information from animal identifiers
conditions = cellfun(@(x) x(end), animal_list_used_after_analysis, 'UniformOutput', false);

% Get list of all available conditions
all_conditions = unique(conditions);
logger(sprintf('Available conditions: %s', strjoin(all_conditions, ', ')), 'INFO');


condition1 = 'S'; 
condition2 = 'F';

% Filter frames for each condition
condition1_mask = strcmp(conditions, condition1);
condition2_mask = strcmp(conditions, condition2);

% Extract animal names for each condition
condition1_animal_names = animal_names(condition1_mask);
condition2_animal_names = animal_names(condition2_mask);

% Extract features for each condition
condition1_features = feature_matrix(condition1_mask, :);
condition2_features = feature_matrix(condition2_mask, :);

% Get unique animals in each condition
condition1_unique_animals = unique(condition1_animal_names);
condition2_unique_animals = unique(condition2_animal_names);

num_condition1_animals = length(condition1_unique_animals);
num_condition2_animals = length(condition2_unique_animals);

logger(sprintf('Found %d animals in condition %s and %d animals in condition %s', ...
    num_condition1_animals, condition1, num_condition2_animals, condition2), 'INFO');

% Initialize arrays to store cosine similarity values for each animal
condition1_cosine_sim = zeros(num_condition1_animals, 1);
condition2_cosine_sim = zeros(num_condition2_animals, 1);

% Calculate within-condition cosine similarity for each animal in condition 1
for i = 1:num_condition1_animals
    animal = condition1_unique_animals{i};
    animal_idx = strcmp(condition1_animal_names, animal);
    
    if sum(animal_idx) < 2
        logger(sprintf('Skipping animal %s - insufficient data points (%d)', ...
            animal, sum(animal_idx)), 'WARNING');
        condition1_cosine_sim(i) = NaN;
        continue;
    end
    
    % Get features for this animal
    animal_features = condition1_features(animal_idx, :);
    
    % Calculate pairwise cosine similarity within this animal's behaviors
    n_samples = size(animal_features, 1);
    pair_count = 0;
    total_similarity = 0;
    
    % Compare each pair of behavior samples
    for j = 1:n_samples
        for k = j+1:n_samples
            % Compute cosine similarity between feature vectors
            sim = dot(animal_features(j,:), animal_features(k,:)) / ...
                  (norm(animal_features(j,:)) * norm(animal_features(k,:)));
            total_similarity = total_similarity + sim;
            pair_count = pair_count + 1;
        end
    end
    
    % Calculate average similarity if we have pairs
    if pair_count > 0
        condition1_cosine_sim(i) = total_similarity / pair_count;
    else
        condition1_cosine_sim(i) = NaN;
    end
    
    logger(sprintf('Animal %s (%s): Average cosine similarity = %.4f (%d pairs)', ...
        animal, condition1, condition1_cosine_sim(i), pair_count), 'INFO');
end

% Calculate within-condition cosine similarity for each animal in condition 2
for i = 1:num_condition2_animals
    animal = condition2_unique_animals{i};
    animal_idx = strcmp(condition2_animal_names, animal);
    
    if sum(animal_idx) < 2
        logger(sprintf('Skipping animal %s - insufficient data points (%d)', ...
            animal, sum(animal_idx)), 'WARNING');
        condition2_cosine_sim(i) = NaN;
        continue;
    end
    
    % Get features for this animal
    animal_features = condition2_features(animal_idx, :);
    
    % Calculate pairwise cosine similarity within this animal's behaviors
    n_samples = size(animal_features, 1);
    pair_count = 0;
    total_similarity = 0;
    
    % Compare each pair of behavior samples
    for j = 1:n_samples
        for k = j+1:n_samples
            % Compute cosine similarity between feature vectors
            sim = dot(animal_features(j,:), animal_features(k,:)) / ...
                  (norm(animal_features(j,:)) * norm(animal_features(k,:)));
            total_similarity = total_similarity + sim;
            pair_count = pair_count + 1;
        end
    end
    
    % Calculate average similarity if we have pairs
    if pair_count > 0
        condition2_cosine_sim(i) = total_similarity / pair_count;
    else
        condition2_cosine_sim(i) = NaN;
    end
    
    logger(sprintf('Animal %s (%s): Average cosine similarity = %.4f (%d pairs)', ...
        animal, condition2, condition2_cosine_sim(i), pair_count), 'INFO');
end

% Remove NaN values
condition1_cosine_sim = condition1_cosine_sim(~isnan(condition1_cosine_sim));
condition2_cosine_sim = condition2_cosine_sim(~isnan(condition2_cosine_sim));

% Calculate mean and standard error of the mean (SEM) for each condition
condition1_mean = mean(condition1_cosine_sim);
condition1_sem = std(condition1_cosine_sim) / sqrt(length(condition1_cosine_sim));

condition2_mean = mean(condition2_cosine_sim);
condition2_sem = std(condition2_cosine_sim) / sqrt(length(condition2_cosine_sim));

% Perform statistical test (t-test)
[h, p] = ttest2(condition1_cosine_sim, condition2_cosine_sim);
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

% Create a bar plot similar to the figure
fig_cosine = figure('Name', 'Cosine Similarity Between Feature Vectors', 'Color', 'w', ...
    'Position', [500, 500, 400, 400], 'Visible', visualize);

% Define condition names for the plot
condition_names = {['S' newline 'F'], 'Sham'};

% Define bar positions
bar_positions = 1:2;
bar_width = 0.7;

% Plot bars
hold on;
b1 = bar(bar_positions(1), condition1_mean, bar_width, 'FaceColor', [0.8, 0.2, 0.2]); % Red for DLS lesion
b2 = bar(bar_positions(2), condition2_mean, bar_width, 'FaceColor', [0.5, 0.5, 0.5]); % Gray for Sham

% Add individual data points (scatter)
scatter(repmat(bar_positions(1), length(condition1_cosine_sim), 1) + (rand(length(condition1_cosine_sim), 1)-0.5)*0.2, ...
    condition1_cosine_sim, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.6);

scatter(repmat(bar_positions(2), length(condition2_cosine_sim), 1) + (rand(length(condition2_cosine_sim), 1)-0.5)*0.2, ...
    condition2_cosine_sim, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.6);

% Add error bars (SEM)
errorbar(bar_positions(1), condition1_mean, condition1_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
errorbar(bar_positions(2), condition2_mean, condition2_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Add significance bar if there's a significant difference
if ~strcmp(significance_label, 'NS')
    % Calculate y position for significance bar
    y_max = max([condition1_mean + condition1_sem, condition2_mean + condition2_sem]) * 1.1;
    
    % Plot significance bar
    plot([bar_positions(1), bar_positions(2)], [y_max, y_max], 'k', 'LineWidth', 1.5);
    text(mean(bar_positions), y_max*1.05, significance_label, 'HorizontalAlignment', 'center', 'FontSize', 14);
else
    % Just show NS text above the bars
    y_max = max([condition1_mean + condition1_sem, condition2_mean + condition2_sem]) * 1.1;
    plot([bar_positions(1), bar_positions(2)], [y_max, y_max], 'k', 'LineWidth', 1.5);
    text(mean(bar_positions), y_max*1.05, 'NS', 'HorizontalAlignment', 'center', 'FontSize', 14);
end

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
logger('Cosine similarity statistics:', 'INFO');
logger(sprintf('DLS lesion (%s): Mean = %.4f, SEM = %.4f, n = %d', ...
    condition1, condition1_mean, condition1_sem, length(condition1_cosine_sim)), 'INFO');
logger(sprintf('Sham (%s): Mean = %.4f, SEM = %.4f, n = %d', ...
    condition2, condition2_mean, condition2_sem, length(condition2_cosine_sim)), 'INFO');
logger(sprintf('t-test: p = %.4f (%s)', p, significance_label), 'INFO');

% Export the figure if needed
if do_export
    logger('Exporting cosine similarity figure to export folder', 'INFO');
    export_name = 'cosine_similarity_between_feature_vectors';
    exportgraphics(fig_cosine, fullfile(export_folder, [export_name '.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
end

% %% Alternative: Create a cosine similarity matrix to visualize between all conditions
% logger('Creating cosine similarity matrix between all conditions', 'INFO');
% 
% % Get all unique conditions
% all_conditions = unique(conditions);
% num_conditions = length(all_conditions);
% 
% % Initialize matrix to store mean cosine similarity between conditions
% condition_similarity_matrix = zeros(num_conditions, num_conditions);
% condition_similarity_sem_matrix = zeros(num_conditions, num_conditions);
% condition_sample_counts = zeros(num_conditions, num_conditions);
% 
% % Calculate cosine similarity between all pairs of conditions
% for i = 1:num_conditions
%     cond_i = all_conditions{i};
%     mask_i = strcmp(conditions, cond_i);
%     features_i = feature_matrix(mask_i, :);
%     animals_i = animal_names(mask_i);
% 
%     for j = 1:num_conditions
%         cond_j = all_conditions{j};
%         mask_j = strcmp(conditions, cond_j);
%         features_j = feature_matrix(mask_j, :);
%         animals_j = animal_names(mask_j);
% 
%         % Calculate similarities between all pairs of feature vectors
%         similarities = zeros(1000, 1); % Pre-allocate with reasonable size
%         sim_count = 0;
% 
%         % If comparing the same condition, calculate within-condition similarity
%         if i == j
%             % Sample a subset of comparisons for efficiency if dataset is large
%             n_samples_i = size(features_i, 1);
%             max_comparisons = 5000; % Cap to avoid excessive computation
% 
%             if n_samples_i > 100
%                 % Sample random pairs instead of exhaustive computation
%                 for comp = 1:min(max_comparisons, n_samples_i*(n_samples_i-1)/2)
%                     idx1 = randi(n_samples_i);
%                     idx2 = randi(n_samples_i);
% 
%                     % Ensure we're not comparing a vector with itself
%                     while idx2 == idx1
%                         idx2 = randi(n_samples_i);
%                     end
% 
%                     vec1 = features_i(idx1, :);
%                     vec2 = features_i(idx2, :);
% 
%                     % Calculate cosine similarity
%                     sim = dot(vec1, vec2) / (norm(vec1) * norm(vec2));
% 
%                     sim_count = sim_count + 1;
%                     similarities(sim_count) = sim;
%                 end
%             else
%                 % For smaller datasets, calculate all pairwise comparisons
%                 for idx1 = 1:n_samples_i
%                     for idx2 = (idx1+1):n_samples_i
%                         vec1 = features_i(idx1, :);
%                         vec2 = features_i(idx2, :);
% 
%                         % Calculate cosine similarity
%                         sim = dot(vec1, vec2) / (norm(vec1) * norm(vec2));
% 
%                         sim_count = sim_count + 1;
%                         similarities(sim_count) = sim;
%                     end
%                 end
%             end
%         else
%             % Compare between different conditions
%             % Sample a subset of comparisons for efficiency if dataset is large
%             n_samples_i = size(features_i, 1);
%             n_samples_j = size(features_j, 1);
%             max_comparisons = 5000; % Cap to avoid excessive computation
% 
%             if n_samples_i*n_samples_j > max_comparisons
%                 % Sample random pairs instead of exhaustive computation
%                 for comp = 1:max_comparisons
%                     idx1 = randi(n_samples_i);
%                     idx2 = randi(n_samples_j);
% 
%                     vec1 = features_i(idx1, :);
%                     vec2 = features_j(idx2, :);
% 
%                     % Calculate cosine similarity
%                     sim = dot(vec1, vec2) / (norm(vec1) * norm(vec2));
% 
%                     sim_count = sim_count + 1;
%                     similarities(sim_count) = sim;
%                 end
%             else
%                 % For smaller datasets, calculate all pairwise comparisons
%                 for idx1 = 1:n_samples_i
%                     for idx2 = 1:n_samples_j
%                         vec1 = features_i(idx1, :);
%                         vec2 = features_j(idx2, :);
% 
%                         % Calculate cosine similarity
%                         sim = dot(vec1, vec2) / (norm(vec1) * norm(vec2));
% 
%                         sim_count = sim_count + 1;
%                         similarities(sim_count) = sim;
%                     end
%                 end
%             end
%         end
% 
%         % Trim to actual size
%         similarities = similarities(1:sim_count);
% 
%         % Store results
%         condition_similarity_matrix(i, j) = mean(similarities);
%         condition_similarity_sem_matrix(i, j) = std(similarities) / sqrt(sim_count);
%         condition_sample_counts(i, j) = sim_count;
% 
%         logger(sprintf('Cosine similarity between %s and %s: %.4f ± %.4f (n=%d)', ...
%             cond_i, cond_j, condition_similarity_matrix(i, j), ...
%             condition_similarity_sem_matrix(i, j), sim_count), 'INFO');
%     end
% end
% 
% % Create a heatmap visualization of the similarity matrix
% fig_matrix = figure('Name', 'Cosine Similarity Matrix Between Conditions', 'Color', 'w', ...
%     'Position', [600, 600, 500, 400], 'Visible', visualize);
% 
% % Create the heatmap
% imagesc(condition_similarity_matrix);
% colormap(viridis(256)); % Using viridis colormap for better visualization
% colorbar;
% caxis([0, max(condition_similarity_matrix(:))]); % Set color scale from 0 to max
% 
% % Add condition labels
% xticks(1:num_conditions);
% yticks(1:num_conditions);
% xticklabels(all_conditions);
% yticklabels(all_conditions);
% 
% % Add text annotations with values
% for i = 1:num_conditions
%     for j = 1:num_conditions
%         text(j, i, sprintf('%.3f', condition_similarity_matrix(i, j)), ...
%             'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
%     end
% end
% 
% title('Cosine Similarity Between Conditions', 'FontSize', 14);
% xlabel('Condition');
% ylabel('Condition');
% 
% % Export the matrix figure if needed
% if do_export
%     logger('Exporting cosine similarity matrix to export folder', 'INFO');
%     export_name = 'cosine_similarity_matrix_between_conditions';
%     exportgraphics(fig_matrix, fullfile(export_folder, [export_name '.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
% end
% 
% logger('Cosine similarity analysis completed', 'INFO');

%% BASED ON THE PAPER:


%% New alternative:
%% Analyze within-condition cosine similarity and compare across conditions (S vs F)
logger('Analyzing within-condition cosine similarity and comparing across conditions (S vs F)', 'INFO');

% Combine joint features and extra joint features
feature_matrix = [analysisstruct.jt_features, analysisstruct.extra_jt_features];

% Extract animal identifiers and conditions
animal_list = animal_list_used_after_analysis;
conditions = cellfun(@(x) x(end), animal_list, 'UniformOutput', false);

% Extract animal names (without condition suffix)
animal_base_names = cellfun(@(x) x(1:find(x=='_',1)-1), animal_list, 'UniformOutput', false);

% Define the conditions we're comparing
condition1 = 'S';  % Saline
condition2 = 'F';  % Formalin

logger(sprintf('Comparing within-condition similarity between %s and %s', condition1, condition2), 'INFO');

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

% Get cluster assignments for fine behaviors if available
if isfield(analysisstruct, 'cluster_assignments')
    cluster_assignments = analysisstruct.cluster_assignments;
else
    % If no cluster assignments, create dummy assignment (single cluster)
    logger('No cluster assignments found, treating all data as a single behavior', 'WARNING');
    cluster_assignments = ones(size(feature_matrix, 1), 1);
end

unique_clusters = unique(cluster_assignments);
num_clusters = length(unique_clusters);
logger(sprintf('Found %d unique finely parsed behaviors (clusters)', num_clusters), 'INFO');

% Initialize arrays to store mean cosine similarity for each animal in each condition
condition1_animal_similarities = zeros(num_condition1_animals, 1);
condition2_animal_similarities = zeros(num_condition2_animals, 1);

% Process each animal in condition 1 (Saline)
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
    
    % Initialize array to store similarities for each behavior cluster
    behavior_similarities = zeros(num_clusters, 1);
    valid_behavior_count = 0;
    
    % Process each behavior cluster
    for cluster_idx = 1:num_clusters
        cluster = unique_clusters(cluster_idx);
        
        % Find frames for this cluster in this condition for this animal
        cluster_indices = animal_indices(cluster_assignments(animal_indices) == cluster);
        
        if length(cluster_indices) < 2
            % Need at least 2 frames to compute similarity
            continue;
        end
        
        % Process consecutive time points (average over max 5 consecutive frames)
        processed_features = process_consecutive_frames(cluster_indices, feature_matrix, 5);
        
        if size(processed_features, 1) < 2
            % Need at least 2 processed feature vectors to compute similarity
            continue;
        end
        
        % Calculate pairwise cosine similarity within this cluster's feature vectors
        % Using pdist (not pdist2) since we're comparing within the same set
        cosine_distances = pdist(processed_features, 'cosine');
        cosine_similarities = 1 - cosine_distances;
        
        % Average the similarities for this behavior cluster
        avg_similarity = mean(cosine_similarities);
        
        % Store the result if valid
        if ~isnan(avg_similarity)
            behavior_similarities(cluster_idx) = avg_similarity;
            valid_behavior_count = valid_behavior_count + 1;
            
            logger(sprintf('Condition %s, Animal %s, Cluster %d: Average cosine similarity = %.4f (%d comparisons)', ...
                condition1, animal, cluster, avg_similarity, length(cosine_similarities)), 'INFO');
        end
    end
    
    % Calculate the mean similarity across all behaviors for this animal
    if valid_behavior_count > 0
        condition1_animal_similarities(i) = mean(behavior_similarities(behavior_similarities ~= 0));
        
        logger(sprintf('Condition %s, Animal %s: Mean cosine similarity = %.4f across %d behaviors', ...
            condition1, animal, condition1_animal_similarities(i), valid_behavior_count), 'INFO');
    else
        condition1_animal_similarities(i) = NaN;
        logger(sprintf('Animal %s in condition %s: No valid behaviors for comparison', animal, condition1), 'WARNING');
    end
end

% Process each animal in condition 2 (Formalin)
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
    
    % Initialize array to store similarities for each behavior cluster
    behavior_similarities = zeros(num_clusters, 1);
    valid_behavior_count = 0;
    
    % Process each behavior cluster
    for cluster_idx = 1:num_clusters
        cluster = unique_clusters(cluster_idx);
        
        % Find frames for this cluster in this condition for this animal
        cluster_indices = animal_indices(cluster_assignments(animal_indices) == cluster);
        
        if length(cluster_indices) < 2
            % Need at least 2 frames to compute similarity
            continue;
        end
        
        % Process consecutive time points (average over max 5 consecutive frames)
        processed_features = process_consecutive_frames(cluster_indices, feature_matrix, 5);
        
        if size(processed_features, 1) < 2
            % Need at least 2 processed feature vectors to compute similarity
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
                condition2, animal, cluster, avg_similarity, length(cosine_similarities)), 'INFO');
        end
    end
    
    % Calculate the mean similarity across all behaviors for this animal
    if valid_behavior_count > 0
        condition2_animal_similarities(i) = mean(behavior_similarities(behavior_similarities ~= 0));
        
        logger(sprintf('Condition %s, Animal %s: Mean cosine similarity = %.4f across %d behaviors', ...
            condition2, animal, condition2_animal_similarities(i), valid_behavior_count), 'INFO');
    else
        condition2_animal_similarities(i) = NaN;
        logger(sprintf('Animal %s in condition %s: No valid behaviors for comparison', animal, condition2), 'WARNING');
    end
end

% Remove NaN values (animals with no valid data)
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

% Create a bar plot similar to the figure
fig_cosine = figure('Name', 'Within-condition Cosine Similarity Between Feature Vectors', 'Color', 'w', ...
    'Position', [500, 500, 400, 400], 'Visible', visualize);

% Define condition names for the plot
condition_names = {condition1, condition2};

% Define bar positions
bar_positions = 1:2;
bar_width = 0.7;

% Plot bars
hold on;
b1 = bar(bar_positions(1), condition1_mean, bar_width, 'FaceColor', [0.8, 0.2, 0.2]); % Red for Saline
b2 = bar(bar_positions(2), condition2_mean, bar_width, 'FaceColor', [0.5, 0.5, 0.5]); % Gray for Formalin

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
    export_name = sprintf('within_condition_cosine_similarity_%s_vs_%s', condition1, condition2);
    exportgraphics(fig_cosine, fullfile(export_folder, [export_name '.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
end


logger('Saline vs Formalin cosine similarity analysis completed', 'INFO');


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
