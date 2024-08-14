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

% Create empty cell array to fill up the clusters per condition per animal
clusters_pain_no_pain = cell(length(animal_ids), 2); % 2 columns for pain and no pain proportions
for iid = 1: length(animal_ids)
    % animal_ID = animal_ids{iid}; ...

    animal_ID = animal_ids{iid}; % use the animal ID from the list

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
    % proportion: 100*numel(coi_pain_inter)/numel(coi_pain); indicates the
    % proportion of clusters in labeled frames that are present in the
    % automatic detection (predominantF)
    proportion_pain = 100 * numel(coi_pain_inter) / numel(coi_pain);

    % get coi from no pain (S group)
    coi_no_pain = intersect(coi_pain_inter, predominantS); % I hypothesize that there's going to be empty, meaning no overlap between clusters of pain frames and predominant S clusters

    % proportion: 100*numel(coi_no_pain)/numel(coi_pain); indicates the
    % proportion of clusters in labeled frames that are present in the
    % automatic detection (predominantS)
    proportion_no_pain = 100 * numel(coi_no_pain) / numel(coi_pain);
    
    % Gather the clusters of interest (coi_pain and coi_no_pain), coi_pain
    % Gather the proportion of intersect pain frames relative to the
    % coi_pain number
    % in column index 1, coi_no_pain in column index 2
    clusters_pain_no_pain{iid, 1} = proportion_pain;
    clusters_pain_no_pain{iid, 2} = proportion_no_pain;

end % end for loop

% Do statistics and plots
% Convert cell array to matrix for plotting
proportions_matrix = cell2mat(clusters_pain_no_pain);

% Generate bar plots
figure;
bar(proportions_matrix);
legend('Pain Proportion', 'No Pain Proportion');
title('Proportions of Pain and No Pain Clusters');
xlabel('Animal ID');
ylabel('Proportion (%)');
