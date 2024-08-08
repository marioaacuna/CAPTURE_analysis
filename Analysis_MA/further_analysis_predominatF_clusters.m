%% Clean up the predominat frames.
% in order to understand the clusters that span what we call pain, they
% need to be present in at least 75% of the total number of animals under
% the pain condition (in this case F). Later on we can look at the clusters
% that are present only on individual animals, to see animal specific
% poses. 
% this script is intended to analyse thse kid of features.


%% Init
clear;
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% load clusters

load(fullfile(rootpath, 'cluster_analysis_results.mat'), 'results');
predominantF = results.predominantF;
predominantS = results.predominantS;


% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');

% Take the proportions and see which contain zeros
all_F_clusters = find(sum(~ismember(results.cluster_proportions_F, 0))>=2);


pain_clusters = predominantF(ismember(predominantF,all_F_clusters));


plot_identified_clusters(pain_clusters, analysisstruct)