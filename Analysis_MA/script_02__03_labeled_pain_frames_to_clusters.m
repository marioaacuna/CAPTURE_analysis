% Read lalbeled facial pain frames and see what clusters they belong to
% You have to go animal by animal
% you have to create fist the cluster vector, that indicates per frame, to
% what cluster it belongs

%% INIT

clear; close all; clc;
global GC
% load cluster results
rootpath = GC.preprocessing_rootpath;
clusters_struct_file = fullfile(rootpath, 'clusters_struct.mat');

load(fullfile(rootpath, 'cluster_analysis_results.mat'), 'results');

% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');

%% choose animal  
animal_ID = 'AK_553_F';

%% load cluster vector
load(clusters_struct_file, 'clusters_struct')

%% match labelebed frames to clusters
% load pain frames
facial_rootpath = fullfile (rootpath ,'facial_labels');
labels_file = fullfile(facial_rootpath, [animal_ID '_labels.csv']);
T =  readtable(labels_file);
pain_frames = find(T.pain == 1);

% clusters of interest
coi = clusters_struct.(animal_ID);
coi_pain = coi(pain_frames);

pain_idx = ismember(coi_pain, results.predominantF);

unique()
%%
figure, 
histogram(coi_pain(pain_idx))