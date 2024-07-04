close all
clc
clear
%% INIT

% List of animals

overwrite_pred_concat = 0; % in case you want to add more animals, set this to 1
overwrite_ratception = 0; % will do or not mocap
overwrite_MLmatobjfile = 0; % Overwrite extarcted features
overwrite_coefficient = 0; % Overwrite dim red coeffs 
overwrite_zvals = 0; % Overwrite tsne zvals
do_extra_features = 0;

init_frame_rate = 100; % effective frame rate of videos

if ispc
    % so far for these test, we have the data on local drive D
    rootpath = 'D:\CAPTURE';
else
    rootpath = '/Users/mario/Library/Mobile Documents/com~apple~CloudDocs/_todo/SpontaneousPainProject/240704_data'; % I have updated where the files are to gather them more quickly
end

if ~exist("rootpath", 'dir'), mkdir(rootpath); end
 
% animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
% conditions = {'F', 'F', 'S', 'S', 'N', 'S', 'F', 'F', 'F'};

% Define the animals and conditions
animals = {'AK_552', 'AK_553', 'AK_665', 'AK_667'}; % animal ids
conditions = {'0', '1'}; % 0 for condition 'S', and 1 for condition 'F'
dannce_path = 'D:\DANNCE';
% Initialize structures
agg_predictions = struct();
animal_condition_identifier = {};


analysis_filename = fullfile(rootpath,'raw_concat_analysis.mat' );
filename_predictions = fullfile(rootpath, 'agg_predictions.mat');

if ~any(exist(filename_predictions, 'file'))
    run_pred_concat = 1 ;
else
    run_pred_concat = 0;
end

if run_pred_concat || overwrite_pred_concat
    disp('%% Running concatenation')
    % Initialize the aggregate predictions structure
    agg_predictions = struct();
    animal_condition_identifier = {};
    % initialize frames id
    % agg_predictions.frames_id = [];
    % Current offset for tracking frames
    offset = 0;

    % % Loop through each animal to aggregate predictions
    % for i = 1:length(animal_list)
    %     animal_id = animal_list{i};
    % 
    %     % Load predictions for the current animal
    %     load_path = fullfile(rootpath, animal_id, 'predictions.mat');
    %     load(load_path); % Assuming the loaded structure is called 'predictions'
    % 
    %     body_parts = fieldnames(predictions); % Get the list of body parts
    % 
    %     for j = 1:length(body_parts)
    %         part_name = body_parts{j};
    % 
    %         if isfield(agg_predictions, part_name)
    %             if ~strcmp(part_name, 'sampleID')
    %                 % Concatenate if the field already exists in agg_predictions
    %                 agg_predictions.(part_name) = cat(1, agg_predictions.(part_name), predictions.(part_name));
    %             else
    %                 agg_predictions.(part_name) = cat(2, agg_predictions.(part_name), predictions.(part_name) + 8 + agg_predictions.(part_name)(end));
    %             end
    %         else
    %             % Otherwise, initialize the field in agg_predictions
    %             agg_predictions.(part_name) = predictions.(part_name);
    %         end
    %     end
    %     % store the frames id
    %     % agg_predictions.frames_id = [agg_predictions.frames_id; [1:1:length(agg_predictions.sampleID)]'];
    %     % Update the animal_frames_identifier
    %     num_frames_for_current_animal = size(predictions.(body_parts{1}), 1);
    %     animal_frames_identifier = [animal_frames_identifier; repmat({animal_id}, num_frames_for_current_animal, 1)];
    % end

    % Loop through each animal
    for i = 1:length(animals)
        animal_id = animals{i};

        % Loop through each condition for the current animal
        for j = 1:length(conditions)
            condition = conditions{j};

            % Construct the path to the predictions file
            load_path = fullfile(dannce_path, animal_id, condition,'DANNCE\predict_results', 'predictions.mat');

            % Check if the file exists
            if ~exist(load_path, 'file')
                warning('File not found: %s', load_path);
                continue;
            end

            % Load predictions for the current animal and condition
            load(load_path); % Assuming the loaded structure is called 'predictions'

            body_parts = fieldnames(predictions); % Get the list of body parts

            % Process each body part
            for k = 1:length(body_parts)
                part_name = body_parts{k};

                if isfield(agg_predictions, part_name)
                    if ~strcmp(part_name, 'sampleID')
                        % Concatenate if the field already exists in agg_predictions
                        agg_predictions.(part_name) = cat(1, agg_predictions.(part_name), predictions.(part_name));
                    else
                        % For sampleID, adjust the values before concatenating
                        last_sample_id = 0;
                        if ~isempty(agg_predictions.(part_name))
                            last_sample_id = agg_predictions.(part_name)(end);
                        end
                        agg_predictions.(part_name) = cat(2, agg_predictions.(part_name), predictions.(part_name) + last_sample_id);
                    end
                else
                    % Otherwise, initialize the field in agg_predictions
                    agg_predictions.(part_name) = predictions.(part_name);
                end
            end

            % Update the animal_condition_identifier
            num_frames_for_current_combo = size(predictions.(body_parts{1}), 1);
            condition_letter = condition_to_letter(condition);
            animal_condition_identifier = [animal_condition_identifier; repmat({sprintf('%s_%s', animal_id, condition_letter)}, num_frames_for_current_combo, 1)];
        end
    end
    % At this point, agg_predictions contains concatenated body part data for all animals
    % and animal_frames_identifier indicates the animal for each frame.


    % save the data
    % predictions = struct();
    predictions = agg_predictions;
    % predictions.predictions = agg_predictions;
    save(filename_predictions, 'predictions', 'animal_condition_identifier')
else
    disp('predictions previously concatenated, now loading them')
    load(filename_predictions)
end

%%
long_animal_frames_identifier = repelem(animal_condition_identifier,3);
%% INIT ratception procedure
filename_ratception = fullfile(rootpath, "ratception_prediction.mat");

animal_name = 'mario_mouse22';

input_params = struct();
input_params.SpineF_marker = 'SpineF';
input_params.SpineM_marker = 'SpineM';
% input_params.repfactor = 300/30;
input_params.repfactor = round(300/init_frame_rate);
input_params.conversion_factor = 1;

% Run prepro if it doesn't exist
if ~exist(filename_ratception, "file") || overwrite_ratception
    disp('%% Running Pre-Pro %% ')
    ratception_struct = preprocess_dannce(filename_predictions,filename_ratception,animal_name,input_params);
else
    disp('Loading previously analysied prepro data')
   load(filename_ratception);
end

%%
roothpath_CAPTURE = fullfile(rootpath, 'CAPTURE');

coefficient_file = fullfile(roothpath_CAPTURE,'coefficients.mat');
linkname = 'mario_mouse22';
ratname ='myrat';% 'test_mouse';

%% Load Mocapstruct
mocapstruct = ratception_struct;
if do_extra_features
    % In case you want to do some extra features
    savefilename_extra = fullfile(rootpath, 'CAPTURE', 'myextratsnefeature', 'extraMLFeatures.mat');
    eigenposture_save_filder = fullfile(rootpath, 'CAPTURE', 'myextratsnefeature');
end
% savefilename ='myextratsnefeature';
% directory_here = pwd;
%feature filename and whether or not to overwrite

savefilename_features = fullfile(roothpath_CAPTURE);
MLmatobjfile = fullfile(savefilename_features, 'myMLfeatures.mat');

%visualize the mocap data
% animate_markers_nonaligned_fullmovie_demo(mocapstruct,1:10:88000, []);

%% Create behavioral features
%this determines the set of frames to use -- in general if the animal is
%resting for too long it will cause errors
mocapstruct.modular_cluster_properties.clipped_index{8} = 1:size(mocapstruct.aligned_mean_position,1 );
% to control the wavelet parameters, you can change the properties in the
% compute_wl_transform_features file

if ~exist(fullfile(rootpath, 'CAPTURE'),'dir' )
    mkdir(fullfile(rootpath, 'CAPTURE'))
end

if ~exist(MLmatobjfile,'file') || overwrite_MLmatobjfile
    MLmatobj = create_behavioral_features(mocapstruct,coefficient_file,overwrite_coefficient,linkname);
    save(MLmatobjfile, 'MLmatobj', '-v7.3')
else
    disp('Loading ML features')
    MLmatobj = load(MLmatobjfile, 'MLmatobj');
    MLmatobj = MLmatobj.MLmatobj;
end

%%
analysisparams.tsnegranularity = 25;% 50:default ;120; 250(1.5, corr_thr) ; 100(1.5)- very detailed; 100(2.5)-very good
% Mario: it seems that high number of frames (between 100 - 120) yield better results

%subselect a particular set of features
analysisstruct = compute_tsne_features(MLmatobj,mocapstruct,analysisparams);

animal_list_used_after_analysis =  long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});
% analysisstruct.condition_inds = [];

% for inter animal analysis:
animal_list = unique(animal_list_used_after_analysis, 'stable');
cond_inds = zeros(1,length(analysisstruct.condition_inds)); % sorting per animal
for iid = 1:length(animal_list)
    animal_ID = animal_list{iid};
    idx = ismember(animal_list_used_after_analysis, animal_ID);
    cond_inds(idx) = iid;
end

% analysisstruct.condition_inds = cond_inds;

% conditions ie. actua conditions, concatenating animals in same cond
condition_inds = cond_inds; % sorting per condition
for iff = 1:length(animal_list_used_after_analysis)
     animal_ID = animal_list_used_after_analysis{iff};
    condition_inds(iff) = double(endsWith(animal_ID, '_F') +1);
   
end
% analysisstruct.condition_inds = condition_inds;
% this does not work, see  compute_tsne_features.m line 144
%% Plot dim red
zvals_filename = fullfile(roothpath_CAPTURE, 'zvals.mat');
perplexity = 200; % 75;

if ~exist(zvals_filename, 'file') || overwrite_zvals

    %run tsne
    disp('%% Running TSNE %%')
    zvals = tsne(analysisstruct.jt_features, "Perplexity",perplexity, 'Exaggeration', 20,'verbose',1,'LearnRate', 1200); %perplexity 90 works well too (less nr of clusters), but maybe not recommended due to few nr of frames (see length(analysisstruct.jt_features))
    % save zvals to then read later if necessary
    save(zvals_filename, 'zvals','-mat')
else
    disp(' Loading TSNE zvals')
    load(zvals_filename)
end



disp('Done')

figure(1)
gscatter(zvals(:,1), zvals(:,2), cond_inds)
% plot(zvals(:,1),zvals(:,2),'ob','MarkerFaceColor','b', 'MarkerSize',2)
title({['Granu: ',num2str(analysisparams.tsnegranularity)], ['Perp: ', num2str(perplexity)]})
set(gcf,'Position',([100 100 1100 1100]))
set(gcf, 'color', 'w')

%% clustering parameters

disp('%% INIT clustering %%')
analysisstruct.zValues = zvals;
analysisstruct.params.density_res = 1001; %resolution of the map
analysisstruct.params.density_width = 2;%default:3; 2; %density kernel in tsne space 0.5 for pca; 2.5 tsne .  !! 5 so far works well
analysisstruct.params.expansion_factor = 1.1; %add a little room to the map after kernel smoothing
analysisstruct.params.density_threshold = 1*10^(-5); %remove regions in plots with low density
analysisstruct.matchedconds = {[unique(analysisstruct.condition_inds)]}; %if running over multiple conditions
analysisstruct.conditions_to_run = [unique(analysisstruct.condition_inds)];
analysisstruct.tsnegranularity = analysisparams.tsnegranularity;

params.reorder=1;
analysisstruct = compute_analysis_clusters_demo(analysisstruct,params); % check line 248, cluster_tsne_map.m
disp('%% Done clustering %%')

%% behavior plots and movies
analysisstruct.conditionnames = ratname;
analysisstruct.ratnames = ratname;
analysisstruct.filesizes = {size(mocapstruct.aligned_mean_position,1 );};

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
set(h1,'Position',([100 100 1100 1100]))
  
% bird specific axes
axisparams.zlim = ([200 300]);
axisparams.xlim = ([-400 400]);
axisparams.ylim = ([-400 400]);

cluster_figure_filename = fullfile(roothpath_CAPTURE, 'Tsne_clusters.pdf');
export_fig(cluster_figure_filename, '-pdf', h1)

%% run sequence and state analysis
params.do_show_pdistmatrix =1;
params.decimation_factor = 1; %downsample if needed to save on memory
% params.decimation_factor = 5; %downsample if needed to save on memory
params.doclustering = 1;

%clustering parameters
params.corr_threshold = 0.2;% 0.2
params.clustercutoff =0.65;%0.12; % 0.65
analysisstruct.plotdirectory = '';
%timescale to use, in seconds
% params.timescales = [1./4 2];  % this is in minutes, check the find_sequences_state_demo code 
params.timescales = [1./4];  % this is in minutes, check the find_sequences_state_demo code 
% paper 15 s and 12 s -> 1/4, 2 
% params.timescales = [15, 120]; 
% params.timescales = [0.1];%[1./4 2]; 
% params.timescales = [0.05];%[1./4 2]; 


analysisstruct.conditionnames = {'test'};
analysisstruct.ratname = {ratname};
condition =1;

% Uncomment this to to sequence analysis

% hierarchystruct=   find_sequences_states_demo(analysisstruct,condition,params);
% % identofy the unique cls
% [seq_cls, seq_c_idx, seq_r] = unique(hierarchystruct.clustered_behavior{1}, 'stable');

% h=figure(370);
% CL = cell(length(seq_c_idx),1);
% for seq_ic = 1:numel(seq_cls)
%     this_cls = seq_cls(seq_ic);    fprintf('ic = %i - ', this_cls)
%      CL(seq_ic) =  {find(hierarchystruct.clustered_behavior{1}==this_cls)};
%     if this_cls==0,  fprintf('\n'),continue, end
%     animate_markers_nonaligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},...
%         find(hierarchystruct.clustered_behavior{1}==this_cls), h, [], ['ic =  ',num2str(this_cls)]);
% 
% end
%%
%%
[cls, c_idx, r] = unique(analysisstruct.annot_reordered{end}, 'stable');

plot_poses = 1;
if plot_poses
    % h= figure(370);
    % clf;

    fig_poses = figure('pos', [10,300,1500,1900]);
    nclus = numel(cls);
    n_rows = ceil(sqrt(nclus));
    n_cols = ceil(sqrt(nclus));
    for ic = 1:numel(cls)
        subplot(n_rows, n_cols, ic)
        this_cls = cls(ic);
        fprintf('ic = %i - \n', this_cls)
        plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1},...
            find(analysisstruct.annot_reordered{end}==this_cls),['cl nr :  ', num2str(this_cls)]);
        title(this_cls)
    end
end
% save cluster plot
cluster_poses_figure_filename = fullfile(roothpath_CAPTURE, 'Poses_clusters.pdf');
export_fig(cluster_poses_figure_filename, '-pdf', fig_poses)
keyboard
%% save analysis
save(analysis_filename, 'analysisstruct' , '-v7.3')

%% Do statistics -  this needs to be modified
% first upsample the animal_frames_identifier to match the  data. take in
% consideration the factor at the beginning (input.repfactor)

up_sampled_identifier = repelem(animal_condition_identifier, input_params.repfactor);
good_frames_animal_identifier = up_sampled_identifier(analysisstruct.frames_with_good_tracking{1, 1});

 
clusters = analysisstruct.annot_reordered{end};

% Identify the unique clusters
unique_clusters = unique(clusters);
numClusters = numel(unique_clusters);
% Initialize storage for the density of each cluster per condition
cluster_density_condition_1 = zeros(size(unique_clusters));
cluster_density_condition_0 = zeros(size(unique_clusters));

% Initialize matrices to store results
clusterComposition = zeros(numClusters, 2);  % [S_count, F_count]

   
    
for c = 1:length(unique_clusters)
    cluster_index = unique_clusters(c);
    
    % Find all frames belonging to the current cluster
    frames_in_cluster = find(clusters == cluster_index);
    
    % Identify the animal-condition combinations for these frames
    animal_conditions_in_cluster = good_frames_animal_identifier(frames_in_cluster);
    
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

% Check what are the clusters only in condition 1
only_condition_1 = find(difference_frames == max(difference_frames));
to_take = only_condition_1;

figure('pos', [10,300,1500,1900])
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


%% new

totalFrames = sum(clusterComposition, 2);
clusterProportions = clusterComposition ./ totalFrames;

% Identify predominantly associated clusters
threshold = 0.75;  % Define threshold for "predominant" association
predominantS = find(clusterProportions(:, 2) >= threshold);
predominantF = find(clusterProportions(:, 1) >= threshold);

to_take = predominantF;
fig_predominant = figure('pos', [10,300,1500,1900]);
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
% title('Predominant pain frames')
pain_frames_fig = fullfile(roothpath_CAPTURE,'Predominant_pain_frames.pdf');
export_fig(pain_frames_fig, '-pdf', fig_predominant)


%%

% % Assuming clusters are stored in analysisstruct.annot_reordered{end}
% clusters = analysisstruct.annot_reordered{end};
%
% % Create a mapping between animal ID and condition
% animal_to_condition = containers.Map(animal_list, conditions); 
% 
% % Identify the unique clusters
% unique_clusters = unique(clusters, 'stable');
% 
% % Initialize storage for the density of each cluster per condition
% cluster_density_formalin = zeros(size(unique_clusters));
% cluster_density_saline = zeros(size(unique_clusters));
% 
% for c = 1:length(unique_clusters)
%     % current_cluster = unique_clusters(c);
%     % cluster_index = find(unique_clusters == current_cluster);
%      cluster_index = unique_clusters(c);
%     % Find all frames belonging to the current cluster
%     frames_in_cluster = find(cluster_index == clusters);
% 
%     % Identify the animals for these frames
%     animals_in_cluster = good_frames_animal_identifier(frames_in_cluster);
% 
%     % Identify the conditions for these animals
%     conditions_in_cluster = values(animal_to_condition, animals_in_cluster);
% 
%     % Compute density for each condition
%     cluster_density_formalin(cluster_index) = sum(strcmp(conditions_in_cluster, 'F')) ;%/ length(conditions_in_cluster);
%     cluster_density_saline(cluster_index) = sum(strcmp(conditions_in_cluster, 'S')) ;%/ length(conditions_in_cluster);
% end
% 
% % At this point, cluster_density_formalin and cluster_density_saline
% % contain the densities of each cluster for the respective conditions.
% 
% % Calculate the difference in number of frames between formalin and saline for each cluster
% difference_frames = cluster_density_formalin - cluster_density_saline;
% 
% % Create a bar graph
% figure;
% bar(difference_frames);
% title('Difference in Frame Counts: Formalin vs. Saline');
% xlabel('Cluster');
% ylabel('Difference in Frame Counts');
% set(gca, 'XTick', 0:10:length(unique_clusters), 'XTickLabel', arrayfun(@num2str, 0:10:length(unique_clusters), 'UniformOutput', false));
% box off


%% if density, check what are the clusters only in formalin
only_formalin = find(difference_frames ==1);
to_take = only_formalin;
figure('pos', [10,300,1500,1900])
nclus = numel(31);
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


%% animal per animal stats
% Number of animals
num_animals = length(animal_list);

% Unique clusters
% unique_clusters = unique(analysisstruct.annot_reordered);

% Initialize the storage
density_per_animal = zeros(num_animals, length(unique_clusters));

for a = 1:num_animals
    animal_id = animal_list{a};
    frames_of_this_animal = find(strcmp(good_frames_animal_identifier, animal_id));
    
    for c = 1:length(unique_clusters)
        cluster_id = unique_clusters(c);
        cluster_index = find(unique_clusters == cluster_id);

        
        % Count the number of frames this animal contributes to this cluster
        density_per_animal(a, cluster_id) = sum(analysisstruct.annot_reordered{1}(frames_of_this_animal) == cluster_id);
    end
end

% Indices for formalin and saline treated animals
formalin_indices = find(strcmp(conditions, 'F'));
saline_indices = find(strcmp(conditions, 'S'));

p_values = zeros(1, length(unique_clusters));

for c = 1:length(unique_clusters)
    this_c = unique_clusters(c);
    % Perform a two-sample t-test comparing formalin and saline densities for this cluster
    [~, p_values(this_c)] = ttest2(density_per_animal(formalin_indices, this_c), density_per_animal(saline_indices, this_c));
end

% Correct for multiple comparisons (optional but recommended)
corrected_p_values = mafdr(p_values, 'BHFDR', true);
sig_cls = find(p_values<0.05);
hold on
plot(sig_cls, repmat(0, size(sig_cls)), '*')
hold off

%% Plot significant cls
figure('pos', [10,300,1500,1900])
nclus = numel(sig_cls);
n_rows = ceil(sqrt(numel(sig_cls)));
n_cols = ceil(sqrt(numel(sig_cls)));
for ic = 1:numel(sig_cls)
    subplot(n_rows, n_cols, ic)
    this_cls = sig_cls(ic);
    fprintf('ic = %i - \n', this_cls)
    plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1},...
        find(analysisstruct.annot_reordered{end}==this_cls),['cl nr :  ', num2str(this_cls)]);
    title(this_cls)
end

%% Plot cls after thr
to_take = find(difference_frames >200);
figure('pos', [10,300,1500,1900])
nclus = numel(31);
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

%% get the videoframes for the clusters
clusters_to_frames_map = containers.Map('KeyType','double','ValueType','any');

for c = 1:length(unique_clusters)
    cluster_id = unique_clusters(c);
    
    % Find indices in upsampled data that belong to this cluster
    indices_in_this_cluster = find(analysisstruct.annot_reordered{1,1} == cluster_id);
    
    % Map these indices to the original frames in mocapstruct.move_frames
    indices_in_original = min(indices_in_this_cluster * 50, length(mocapstruct.move_frames));
    mocap_frames_for_this_cluster = mocapstruct.move_frames(indices_in_original);

    % mocap_frames_for_this_cluster = mocapstruct.move_frames(indices_in_this_cluster * 50);
    
    % Store in the map
    clusters_to_frames_map(cluster_id) = mocap_frames_for_this_cluster;
end

% Adjusted Mapping mocapstruct.move_frames to specific animals:
cluster_to_animal_map = containers.Map('KeyType','double','ValueType','any');

for c = 1:length(unique_clusters)
    cluster_id = unique_clusters(c);
    
    frames_for_this_cluster = clusters_to_frames_map(cluster_id);
    
    % Determine the indices in animal_frames_identifier that these frames correspond to
    indices_in_animal_identifier = frames_for_this_cluster;  % This gives the correct index in animal_frames_identifier
    
    % Ensure indices do not exceed the length of animal_frames_identifier
    valid_indices = min(indices_in_animal_identifier, length(animal_frames_identifier));
    
    % Getting the animal ids for these indices
    animals_for_this_cluster = animal_condition_identifier(valid_indices);
    
    cluster_to_animal_map(cluster_id) = animals_for_this_cluster;
end

%% get the frames per cluster per animal
% Initialize a structure to store cluster-animal-frame relationship
cluster_animal_frames = struct();

for c = 1:length(unique_clusters)
    cluster_id = unique_clusters(c);
    
    % Step 1: Get frames associated with this cluster
    frames_for_this_cluster = clusters_to_frames_map(cluster_id);
    
    % Step 2: Get animals associated with these frames
    animals_for_this_cluster = cluster_to_animal_map(cluster_id);
    
    % Store each animal's frames separately
    for a = 1:length(animals_for_this_cluster)
        animal_id = animals_for_this_cluster{a};
        
        % Use logical indexing to get frames specifically for this animal
        specific_frames_for_animal = frames_for_this_cluster(strcmp(animals_for_this_cluster, animal_id));
        
        % Check if this cluster already exists in the structure
        if ~isfield(cluster_animal_frames, ['cluster_' num2str(cluster_id)])
            cluster_animal_frames.(['cluster_' num2str(cluster_id)]) = struct();
        end
        
        % Step 3: Store frames for this animal under this cluster
        field_name_for_animal = ['animal_' animal_id];
        cluster_animal_frames.(['cluster_' num2str(cluster_id)]).(field_name_for_animal) = specific_frames_for_animal;
    end
end



%% get the frames
%example
frames_for_animal_326_in_cluster_5 = cluster_animal_frames.cluster_5.('animal_334');

%% Solve problem with initial frame ids
% Assuming animal_list contains the IDs of animals in the order they were concatenated
animal_list = {'326', '327', '328', '330', '332_training', '332', '334'};

% Calculate the cumulative frames for each animal
cumulative_frames = zeros(1, length(animal_list));
for i = 2:length(animal_list)
    animal_id = animal_list{i-1};
    cumulative_frames(i) = cumulative_frames(i-1) + length(predictions.(animal_id)); %check this
end

% Adjust the frame indices for each animal and each cluster
adjusted_cluster_animal_frames = struct();

for cluster_id = unique_clusters
    for animal_id = animal_list
        % Extract current frames for this cluster and animal
        current_frames = cluster_animal_frames.(['cluster_' num2str(cluster_id)]).(['animal_' animal_id]);
        
        % Adjust these frames to get original frame indices for the animal
        animal_idx = find(strcmp(animal_list, animal_id));
        adjusted_frames = current_frames - cumulative_frames(animal_idx);
        
        % Store these adjusted frames
        adjusted_cluster_animal_frames.(['cluster_' num2str(cluster_id)]).(['animal_' animal_id]) = adjusted_frames;
    end
end



%% Functions
function letter = condition_to_letter(condition)
    if strcmp(condition, '0')
        letter = 'S';
    elseif strcmp(condition, '1')
        letter = 'F';
    else
        error('Unknown condition: %s', condition);
    end
end