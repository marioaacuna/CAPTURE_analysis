%% Initialization
logger('Starting behavioral cluster analysis script', 'INFO');
clear;
close all;
clc;
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
% export_folder = '~/Desktop/figs_presentation_painAI';


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
animal_list_used_after_analysis =  long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});



% Plot t-SNE maps per condition, using scatter
logger('Plotting t-SNE maps per condition', 'INFO');


%% Extract condition identifiers and set up colors
% Extract last character of each identifier to determine condition
conditions = cellfun(@(x) x(end), animal_list_used_after_analysis, 'UniformOutput', false);

% Define colors for each condition
color_map = containers.Map();
color_map('B') = [0, 0.4470, 0.7410];     % Blue
color_map('F') = [0.8500, 0.3250, 0.0980]; % Red
color_map('H') = [0.9290, 0.6940, 0.1250]; % Yellow
color_map('N') = [0.4940, 0.1840, 0.5560]; % Purple
color_map('S') = [0.4660, 0.6740, 0.1880]; % Green

% Store all z-values for easier access
zvals = analysisstruct.zValues;

%% Create density maps for each condition comparison
% 1. S vs F comparison
logger('Plotting S vs F comparison', 'INFO');
fig_SF_density = figure('Name', 'Density Maps: S vs F', 'Color', 'w', 'Position', [100, 100, 800, 400], 'Visible', visualize);

% S condition
subplot(1, 2, 1);
idx_S = strcmp(conditions, 'S');
h_S = gca;
set(h_S, 'Color', 'w');
plotdensitymaps({zvals(idx_S,:)}, 1, h_S, analysisstruct.params.density_width, ...
    max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor, analysisstruct.params.density_res);
title('Condition S');
axis square

% F condition
subplot(1, 2, 2);
idx_F = strcmp(conditions, 'F');
h_F = gca;
set(h_F, 'Color', 'w');
plotdensitymaps({zvals(idx_F,:)}, 1, h_F, analysisstruct.params.density_width, ...
    max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor, analysisstruct.params.density_res);
title('Condition F');
axis square

% Export the figure if needed
if do_export
    logger(['Exporting S vs F density map to: ' export_folder], 'INFO');
    exportgraphics(fig_SF_density, [export_folder '/density_map_SF_comparison.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end

% 2. H vs N comparison
logger('Plotting H vs N comparison', 'INFO');
fig_HN_density = figure('Name', 'Density Maps: H vs N', 'Color', 'w', 'Position', [100, 100, 800, 400], 'Visible', visualize);

% H condition
subplot(1, 2, 1);
idx_H = strcmp(conditions, 'H');
h_H = gca;
set(h_H, 'Color', 'w');
plotdensitymaps({zvals(idx_H,:)}, 1, h_H, analysisstruct.params.density_width, ...
    max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor, analysisstruct.params.density_res);
title('Condition H');
axis square

% N condition
subplot(1, 2, 2);
idx_N = strcmp(conditions, 'N');
h_N = gca;
set(h_N, 'Color', 'w');
plotdensitymaps({zvals(idx_N,:)}, 1, h_N, analysisstruct.params.density_width, ...
    max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor, analysisstruct.params.density_res);
title('Condition N');
axis square

% Export the figure if needed
if do_export
    logger(['Exporting H vs N density map to: ' export_folder], 'INFO');
    exportgraphics(fig_HN_density, [export_folder '/density_map_HN_comparison.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end

% 3. B vs S vs H comparison
logger('Plotting B vs S vs H comparison', 'INFO');
fig_BSH_density = figure('Name', 'Density Maps: B vs S vs H', 'Color', 'w', 'Position', [100, 100, 1200, 400], 'Visible', visualize);

% B condition
subplot(1, 3, 1);
idx_B = strcmp(conditions, 'B');
h_B = gca;
set(h_B, 'Color', 'w');
plotdensitymaps({zvals(idx_B,:)}, 1, h_B, analysisstruct.params.density_width, ...
    max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor, analysisstruct.params.density_res);
title('Condition B');
axis square

% S condition (reusing idx_S from above)
subplot(1, 3, 2);
h_S2 = gca;
set(h_S2, 'Color', 'w');
plotdensitymaps({zvals(idx_S,:)}, 1, h_S2, analysisstruct.params.density_width, ...
    max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor, analysisstruct.params.density_res);
title('Condition S');
axis square

% H condition (reusing idx_H from above)
subplot(1, 3, 3);
h_H2 = gca;
set(h_H2, 'Color', 'w');
plotdensitymaps({zvals(idx_H,:)}, 1, h_H2, analysisstruct.params.density_width, ...
    max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor, analysisstruct.params.density_res);
title('Condition H');
axis square

% Export the figure if needed
if do_export
    logger(['Exporting B vs S vs H density map to: ' export_folder], 'INFO');
    exportgraphics(fig_BSH_density, [export_folder '/density_map_BSH_comparison.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end

%% Create scatter plots for each comparison
% 1. S vs F scatter plot
logger('Creating scatter plot for S vs F', 'INFO');
fig_SF_scatter = figure('Name', 'Scatter Plot: S vs F', 'Color', 'w', 'Visible', visualize);
hold on;
scatter(zvals(idx_S,1), zvals(idx_S,2), 10, color_map('S'), 'Marker', '.', 'DisplayName', 'S');
scatter(zvals(idx_F,1), zvals(idx_F,2), 10, color_map('F'), 'Marker', '.', 'DisplayName', 'F');
hold off;
legend('Location', 'best');
title('S vs F t-SNE Map');
xlabel('t-SNE Dimension 1');
ylabel('t-SNE Dimension 2');
axis equal tight;

% Export the figure if needed
if do_export
    logger(['Exporting S vs F scatter plot to: ' export_folder], 'INFO');
    exportgraphics(fig_SF_scatter, [export_folder '/scatter_plot_SF_comparison.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end

% 2. H vs N scatter plot
logger('Creating scatter plot for H vs N', 'INFO');
fig_HN_scatter = figure('Name', 'Scatter Plot: H vs N', 'Color', 'w', 'Visible', visualize);
hold on;
scatter(zvals(idx_H,1), zvals(idx_H,2), 10, color_map('H'), 'Marker', '.', 'DisplayName', 'H');
scatter(zvals(idx_N,1), zvals(idx_N,2), 10, color_map('N'), 'Marker', '.', 'DisplayName', 'N');
hold off;
legend('Location', 'best');
title('H vs N t-SNE Map');
xlabel('t-SNE Dimension 1');
ylabel('t-SNE Dimension 2');
axis equal tight;

% Export the figure if needed
if do_export
    logger(['Exporting H vs N scatter plot to: ' export_folder], 'INFO');
    exportgraphics(fig_HN_scatter, [export_folder '/scatter_plot_HN_comparison.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end

% 3. B vs S vs H scatter plot
logger('Creating scatter plot for B vs S vs H', 'INFO');
fig_BSH_scatter = figure('Name', 'Scatter Plot: B vs S vs H', 'Color', 'w', 'Visible', visualize);
hold on;
scatter(zvals(idx_B,1), zvals(idx_B,2), 10, color_map('B'), 'Marker', '.', 'DisplayName', 'B');
scatter(zvals(idx_S,1), zvals(idx_S,2), 10, color_map('S'), 'Marker', '.', 'DisplayName', 'S');
scatter(zvals(idx_H,1), zvals(idx_H,2), 10, color_map('H'), 'Marker', '.', 'DisplayName', 'H');
hold off;
legend('Location', 'best');
title('B vs S vs H t-SNE Map');
xlabel('t-SNE Dimension 1');
ylabel('t-SNE Dimension 2');
axis equal tight;

% Export the figure if needed
if do_export
    logger(['Exporting B vs S vs H scatter plot to: ' export_folder], 'INFO');
    exportgraphics(fig_BSH_scatter, [export_folder '/scatter_plot_BSH_comparison.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
end

logger('Plotting complete', 'INFO');