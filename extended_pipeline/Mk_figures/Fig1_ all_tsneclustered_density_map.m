%% Initialization
logger('Starting behavioral cluster analysis script', 'INFO');
clear;
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;


% Export folder
export_folder = fullfile(GC.temp_root, 'figs');
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


%% Plots
sizes = [100 100 1100 1100];
%% plot a tsne map -- see plotting script for parameter definitions
h1=figure(609);
clf;
params.nameplot=1;
params.density_plot =0;
params.watershed = 1;
params.sorted = 1;
params.markersize = 1;
params.coarseboundary =0;
params.do_coarse = 0;
% plot tsne
plot_clustercolored_tsne(analysisstruct,1,params.watershed,h1,params)
set(h1,'Position',(sizes))
  
% bird specific axes
axisparams.zlim = ([200 300]);
axisparams.xlim = ([-400 400]);
axisparams.ylim = ([-400 400]);

cluster_figure_filename = fullfile(GC.figure_folder, 'Tsne_clusters.pdf');

% Export
disp('Exporting tsne-clustered map...')
exportgraphics(h1, [export_folder '/all_Tsne_clustered.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
exportgraphics(h1, [export_folder '/all_Tsne_clustered.jpeg'], 'ContentType', 'vector', 'BackgroundColor', 'none');
disp('done')

%% density
% All
fig_SF_density = figure('Name', 'Density Map', 'Color', 'w', 'Position', sizes, 'Visible', visualize);

% S condition
idx_S = strcmp(conditions, 'S');
h_S = gca;
set(h_S, 'Color', 'w');
plotdensitymaps({zvals}, 1, h_S, analysisstruct.params.density_width, ...
    max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor, analysisstruct.params.density_res);
title('All');
axis square

% Export
disp('Exporting density map...')
exportgraphics(fig_SF_density, [export_folder '/all_density_map.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
disp('done')