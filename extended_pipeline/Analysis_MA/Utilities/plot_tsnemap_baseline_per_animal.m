%% Initialization
clc;
logger('Starting baseline t-SNE analysis per animal script', 'INFO');
clear;
close all;

GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Configuration for visualization and export
debugging = true;  % Set to true for debugging mode
if debugging
    visualize = 'on';  % Show figures during debugging
    do_export = false;  % Don't export during debugging
else
    visualize = 'off';  % Don't show figures in production mode
    do_export = true;   % Export figures in production mode
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

% Plot t-SNE maps for baseline condition per animal
logger('Plotting t-SNE maps for baseline condition per animal', 'INFO');

%% Extract condition identifiers and animal names
% Extract last character of each identifier to determine condition
conditions = cellfun(@(x) x(end), animal_list_used_after_analysis, 'UniformOutput', false);

% Extract animal names (everything before the underscore)
animal_names = cellfun(@(x) x(1:find(x=='_',1)-1), animal_list_used_after_analysis, 'UniformOutput', false);

% Filter for baseline condition only ('B')
idx_baseline = strcmp(conditions, 'B');
baseline_animal_names = animal_names(idx_baseline);
baseline_zvals = analysisstruct.zValues(idx_baseline, :);

% Get unique animal names for baseline condition
unique_animals = unique(baseline_animal_names);
num_animals = length(unique_animals);

logger(['Found ' num2str(num_animals) ' unique animals in baseline condition'], 'INFO');

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

%% Create density maps for each animal in baseline condition
logger('Creating density maps for each animal in baseline condition', 'INFO');

% Calculate subplot layout
n_cols = ceil(sqrt(num_animals));
n_rows = ceil(num_animals / n_cols);

fig_density = figure('Name', 'Density Maps: Baseline per Animal', 'Color', 'w', ...
    'Position', [100, 100, 300*n_cols, 300*n_rows], 'Visible', visualize);

for i = 1:num_animals
    animal = unique_animals{i};
    idx_animal = strcmp(baseline_animal_names, animal);
    
    subplot(n_rows, n_cols, i);
    h_ax = gca;
    set(h_ax, 'Color', 'w');
    
    if sum(idx_animal) > 0  % Check if animal has data points
        plotdensitymaps({baseline_zvals(idx_animal,:)}, 1, h_ax, analysisstruct.params.density_width, ...
            max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor, analysisstruct.params.density_res);
    end
    
    title(['Animal: ' animal]);
    axis square
    
    logger(['Processed density map for animal: ' animal ' (' num2str(sum(idx_animal)) ' frames)'], 'INFO');
end

% Export the density figure if needed
if do_export
    logger(['Exporting baseline density maps per animal to: ' export_folder], 'INFO');
    exportgraphics(fig_density, [export_folder '/density_map_baseline_per_animal.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end

%% Create combined scatter plot with all animals in baseline condition
logger('Creating combined scatter plot for all animals in baseline condition', 'INFO');

fig_scatter = figure('Name', 'Scatter Plot: Baseline All Animals', 'Color', 'w', 'Visible', visualize);
hold on;

% Plot each animal with its specific color
for i = 1:num_animals
    animal = unique_animals{i};
    idx_animal = strcmp(baseline_animal_names, animal);
    
    if sum(idx_animal) > 0  % Check if animal has data points
        scatter(baseline_zvals(idx_animal,1), baseline_zvals(idx_animal,2), 10, ...
            animal_color_map(animal), 'Marker', '.', 'DisplayName', animal);
        
        logger(['Added scatter points for animal: ' animal ' (' num2str(sum(idx_animal)) ' frames)'], 'INFO');
    end
end

hold off;
legend('Location', 'best');
title('Baseline Condition t-SNE Map - All Animals');
xlabel('t-SNE Dimension 1');
ylabel('t-SNE Dimension 2');
axis equal tight;

% Export the scatter figure if needed
if do_export
    logger(['Exporting baseline scatter plot per animal to: ' export_folder], 'INFO');
    exportgraphics(fig_scatter, [export_folder '/scatter_plot_baseline_per_animal.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end

%% Create individual scatter plots for each animal (optional)
logger('Creating individual scatter plots for each animal', 'INFO');

% Calculate subplot layout for individual plots
fig_individual = figure('Name', 'Individual Scatter Plots: Baseline per Animal', 'Color', 'w', ...
    'Position', [200, 200, 300*n_cols, 300*n_rows], 'Visible', visualize);

for i = 1:num_animals
    animal = unique_animals{i};
    idx_animal = strcmp(baseline_animal_names, animal);
    
    subplot(n_rows, n_cols, i);
    
    if sum(idx_animal) > 0  % Check if animal has data points
        scatter(baseline_zvals(idx_animal,1), baseline_zvals(idx_animal,2), 10, ...
            animal_color_map(animal), 'Marker', '.');
    end
    
    title(['Animal: ' animal]);
    xlabel('t-SNE Dimension 1');
    ylabel('t-SNE Dimension 2');
    axis equal tight;
end

% Export the individual plots figure if needed
if do_export
    logger(['Exporting individual baseline scatter plots per animal to: ' export_folder], 'INFO');
    exportgraphics(fig_individual, [export_folder '/scatter_plots_baseline_individual_animals.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end

%% Display summary statistics
logger('=== BASELINE CONDITION SUMMARY ===', 'INFO');
logger(['Total animals in baseline condition: ' num2str(num_animals)], 'INFO');
logger(['Total frames in baseline condition: ' num2str(sum(idx_baseline))], 'INFO');

for i = 1:num_animals
    animal = unique_animals{i};
    idx_animal = strcmp(baseline_animal_names, animal);
    frames_count = sum(idx_animal);
    percentage = (frames_count / sum(idx_baseline)) * 100;
    logger(['Animal ' animal ': ' num2str(frames_count) ' frames (' num2str(percentage, '%.1f') '%)'], 'INFO');
end

logger('Baseline per animal analysis complete', 'INFO');
