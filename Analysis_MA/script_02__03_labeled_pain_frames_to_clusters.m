% Read lalbeled facial pain frames and see what clusters they belong to
% You have to go animal by animal
% you have to create fist the cluster vector, that indicates per frame, to
% what cluster it belongs

%% INIT
clear; close all; clc;

% load cluster results
rootpath = GC.preprocessing_rootpath;
load(fullfile(rootpath, 'cluster_analysis_results.mat'), 'results');
%% choose animal  
animal_ID = 'AK_553_F';

%% load cluster vector
load(clusters_struct_file, 'clusters_struct')

%% match labelebed frames to clusters
% load pain frames
facial_rootpath = 'C:\Users\acuna\OneDrive - Universitaet Bern\Spontaneous_pain_kinematics\data\0_preprocessing\facial_labels';
T =  readtable(fullfile(facial_rootpath, [animal_ID '_labels.csv']));