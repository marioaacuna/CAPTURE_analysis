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

% Load cluster vectors for all the animals

load(clusters_struct_file, 'clusters_struct')

%% iterate therough animals
animal_ids = fieldnames(clusters_struct);
animal_ids(ismember(animal_ids, 'conditions')) = [];

predominantF = results.predominantF;
predominantS = results.predominantS;
n_predominantF = length(predominantF);
n_predominantS = length(predominantS);

% << fill up here >>
% Create empty cell array to fill up the clusters per condition per animal
% clusters_pain_no_pain = cell(0,0); proportions: ... etc
% << until here >>
for iid = 1: length(animal_ids)
    % animal_ID = animal_ids{iid}; ...

    animal_ID = 'AK_553_F'; % example

    % calibrate condition
    if ~endsWith(animal_ID, 'F')
        continue
    end

    %% match labelebed frames to clusters
    % load pain frames
    facial_rootpath = fullfile (rootpath ,'facial_labels');
    labels_file = fullfile(facial_rootpath, [animal_ID '_labels.csv']);
    T =  readtable(labels_file);
    pain_frames = find(T.pain == 1);

    % clusters of interest
    coi = clusters_struct.(animal_ID);
    coi_pain = coi(pain_frames);
    coi_pain_inter = intersect(coi_pain, predominantF); % this depicts how many identified pain frames are present in predominant pain clusters
    % << fill up here
    % proportion: 100*numel(coi_pain_inter)/numel(coi_pain); indicates the
    % proportion of clusters in labeled frames that are present in the
    % automatic detection (predominantF)
    % get coi from no pain (S group)
    coi_no_pain = intersect(coi_pain_inter, predominantS); % I hypothetize that there's going to be empty, meaning no overlap between clusters of pain frames and predominat S clusters
    
    % << until here >>
    
    % << fill up here >>
    % Gather the clusters of interest (coi_pain and coi_no_pain), coi_pain
    % Gather the proportion of intersect pain frames relative to the
    % coi_pain number
    % in column index 1, coi_no_pain in column index 2
    % << until here >>

end % end for loop

% Do statistics and plots
% << fill up here >>

% << until here >>