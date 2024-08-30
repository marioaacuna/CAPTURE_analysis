% Parameters experiments from Oct 2023
% The data was done by Jun between X-Y 2023.
% This is a cross sectional (xs) study, where some animals had saline (S)
% injection and other formalin (F)
%   - version: 30 Aug 2024
%% 1. Get data from server
global GC
xs_animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
xs_conditions = {'F', 'F', 'S', 'S', 'N', 'S', 'F', 'F', 'F'};

% data_in_server_path  = '/Volumes/Mario/Results_Capture'; % this needs to be adjusted to be either in PC or Mac
data_in_server_path = fullfile(GC.preprocessing_rootpath, 'miniscope_batch_behavior_data');
%% Load data
load(fullfile(data_in_server_path, "analysisstruct_clusters.mat"), 'analysisstruct');

% Load predictions
load(fullfile(data_in_server_path, "agg_predictions.mat"), 'predictions', 'animal_condition_identifier');

% Load ratception structure -  not sure where it is yet.
%load(fullfile(data_in_server_path, , 'ratception_struct');

if ~exist('animal_condition_identifier', 'var')
    do_identifier = 1;
else
    do_identifier = 0;
end
if do_identifier
    % animal_frames_identifier contains a vector sting with animal IDs per frame. The animal list and each condition are given by xs_animal_list and xs_conditions
    % Add '_F' or '_S' to the animal IDs in the animal_frames_identifier to create the animal_condition_identifier (if not loaded from file)
    % let's do that in a for loop
    animal_condition_identifier = cell(1, length(animal_frames_identifier));
    for i = 1:length(animal_frames_identifier)
        this_animal = animal_frames_identifier{i};
        if contains(this_animal, xs_animal_list)
            animal_condition_identifier{i} = strcat(this_animal, '_', xs_conditions{strcmp(this_animal, xs_animal_list)});
        else
            animal_condition_identifier{i} = this_animal;
        end
    end

    % Save the animal_condition_identifier in the predictions.mat file
    save(fullfile(data_in_server_path, "agg_predictions.mat"), 'predictions', 'animal_condition_identifier');

end

%% Plot some figures

% remove animal 'N'
to_remove= find(contains(xs_conditions, 'N'));
long_animal_frames_identifier = repelem(animal_condition_identifier,GC.repfactor);

animal_list_used_after_analysis = long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});


% for inter animal analysis:
animal_list = unique(animal_list_used_after_analysis, 'stable');
cond_inds = zeros(1,length(analysisstruct.condition_inds)); % sorting per animal
for iid = 1:length(animal_list)
    animal_ID = animal_list{iid};
    if endsWith(animal_ID, '_N')
        continue
    end
    idx = ismember(animal_list_used_after_analysis, animal_ID);
    cond_inds(idx) = iid;
end


% analysisstruct.condition_inds = cond_inds;

% conditions ie. actua conditions, concatenating animals in same cond
condition_inds = cond_inds; % sorting per condition
for iff = 1:length(animal_list_used_after_analysis)
     animal_ID = animal_list_used_after_analysis{iff};
     if endsWith(animal_ID, '_N')
         continue
     end
    condition_inds(iff) = double(endsWith(animal_ID, '_F') +1);
   
end

conds = unique(condition_inds);

%% Plot the density maps
for icond = 2:3 % cond2: sal, cond 3: F
    fighand_in =figure(444 + icond-1);
    set(fighand_in,'Color','w')    
    
    title(num2str(this_cond))
    this_cond = conds(icond);
    Zvals =analysisstruct.zValues(ismember(condition_inds, this_cond),:);
    % Zvals =analysisstruct.zValues(contains(condition_inds, this_cond),:);
    % subplot(1,2,icond)
     plotdensitymaps({Zvals},1,fighand_in,analysisstruct.params.density_width,...
        max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor,analysisstruct.params.density_res);
    % this will normalize it to all values
end

%% get Predominant clusters
input_params.repfactor = GC.repfactor;
% input_params = analysisstruct.input_params;
upsampled_identifiers = repelem(animal_condition_identifier, input_params.repfactor);

% Get the frames with good tracking
good_frames = analysisstruct.frames_with_good_tracking{1, 1};

% Extract cluster IDs and corresponding identifiers for good frames
cluster_ids = analysisstruct.annot_reordered{end};
frame_identifiers = upsampled_identifiers(good_frames);

% Extract unique animal IDs and conditions
[animal_ids, ~, animal_indices] = unique(cellfun(@(x) x(1:end-2), frame_identifiers, 'UniformOutput', false));
conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);
% get the clusters
unique_clusters = unique(cluster_ids);
num_clusters = length(unique_clusters);
num_animals = length(animal_ids);


% Initialize storage for the density of each cluster per condition
cluster_density_condition_1 = zeros(size(unique_clusters));
cluster_density_condition_0 = zeros(size(unique_clusters));

% Initialize matrices to store results
clusterComposition = zeros(num_clusters, 2);  % [S_count, F_count]

% animal_frames_ids = frame_identifiers(good_frames);
    
for c = 1:length(unique_clusters)
    cluster_index = unique_clusters(c);
    
    % Find all frames belonging to the current cluster
    frames_in_cluster = find(cluster_ids == cluster_index);
    
    % Identify the animal-condition combinations for these frames
    animal_conditions_in_cluster = frame_identifiers(frames_in_cluster);
    
    % Count frames for each condition
    cluster_density_condition_1(c) = sum(cellfun(@(x) endsWith(x, '_F'), animal_conditions_in_cluster));
    cluster_density_condition_0(c) = sum(cellfun(@(x) endsWith(x, '_S'), animal_conditions_in_cluster));
    
    
    clusterComposition(c, 1) =  sum(cellfun(@(x) endsWith(x, '_F'), animal_conditions_in_cluster));
    clusterComposition(c, 2) =  sum(cellfun(@(x) endsWith(x, '_S'), animal_conditions_in_cluster));
end

% Calculate the difference in number of frames between condition 1 and 0 for each cluster
difference_frames = cluster_density_condition_1 - cluster_density_condition_0;

% Create a bar graph
figure;
bar(difference_frames);
title('Difference in Frame Counts: Condition 1 vs. Condition 0');
xlabel('Cluster');
ylabel('Difference in Frame Counts');
set(gca, 'XTick', 0:10:length(unique_clusters), 'XTickLabel', arrayfun(@num2str, 0:10:length(unique_clusters), 'UniformOutput', false));
box off

totalFrames = sum(clusterComposition, 2);
clusterProportions = clusterComposition ./ totalFrames;

% Identify predominantly associated clusters
threshold = 0.99;  % Define threshold for "predominant" association
predominantF = find(clusterProportions(:, 1) >= threshold);
predominantS = find(clusterProportions(:, 2) >= threshold);

to_take = predominantF;
fig_predominant = figure('pos', [10,10,2056,1350], 'color','w');
n_rows = ceil(sqrt(numel(to_take)));
n_cols = ceil(sqrt(numel(to_take)));

for ic = 1:numel(to_take)
    subplot(n_rows, n_cols, ic)
    this_cls = to_take(ic);
    fprintf('ic = %i - \n', this_cls)
    plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1},...
        find(analysisstruct.annot_reordered{end}==this_cls),['cl nr :  ', num2str(this_cls)]);
    title(this_cls)
end
predominat_figure_name = fullfile(GC.figure_folder, 'predominant_F.pdf');
%export_fig(predominat_figure_name, '-pdf', fig_predominant)



