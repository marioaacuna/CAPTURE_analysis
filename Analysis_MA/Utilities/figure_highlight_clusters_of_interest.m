% From script_02__02_behavioral_cluster_anaylsis.m save the cluster of interest (predominant clusters in F)
% and run this script to highlight the clusters of interest in the figure


% Get cthe analysis struct
%% Initialization
clear;
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Get clusters of interest
load(fullfile(rootpath, 'cluster_analysis_results.mat'), 'results');
clusters_of_interest = results.predominantF;



%% Load Data
% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');

% get Parameters
params.nameplot=1;
params.density_plot =0;
params.watershed = 1;
params.sorted = 1;
params.markersize = 1;
params.coarseboundary =0;
params.do_coarse = 0;

% Make figure
iter = 1;
Fig = figure();
plot_clustercolored_tsne_highlight(analysisstruct,1,params.watershed,Fig,params, clusters_of_interest(:)); 

