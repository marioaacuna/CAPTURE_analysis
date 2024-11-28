% So far this script does not work well as it is not finished. Do not use.
% load data
% so far for these test, we have the data on local drive D.
clear
clc
iid = 1;
rootpath = 'D:\test_CAPTURE';
animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
conditions = {'F', 'F', 'S', 'S', 'N', 'S', 'F', 'F', 'F'};
% load mocap
filename_ratception = fullfile(rootpath, animal_list{iid},'ratception_prediction.mat');
load(filename_ratception);
mocapstruct = ratception_struct;
mocapstruct.modular_cluster_properties.clipped_index{8} = 1:size(mocapstruct.aligned_mean_position,1 );

% load MLmatob
MLmatobjfile =  fullfile(rootpath, animal_list{iid}, 'CAPTURE', 'myMLfeatures.mat');
MLmatobj = load(MLmatobjfile, 'MLmatobj');
MLmatobj = MLmatobj.MLmatobj;

%% Run analysis struct if not done already

%if ~exist(...)%%
analysisparams.tsnegranularity = 50;% 120; 250(1.5, corr_thr) ; 100(1.5)- very detailed; 100(2.5)-very good
% Mario: it seems that high number of frames (between 100 - 120) yield better results

%subselect a particular set of features
analysisstruct = compute_tsne_features(MLmatobj,mocapstruct,analysisparams);

%run tsne
disp('%% Running TSNE %%')
zvals = tsne(analysisstruct.jt_features, "Perplexity",30, 'verbose',1); %perplexity 90 works well too (less nr of clusters), but maybe not recommended due to few nr of frames (see length(analysisstruct.jt_features))

disp('%% INIT clustering %%')
analysisstruct.zValues = zvals;
analysisstruct.params.density_res = 1001; %resolution of the map
analysisstruct.params.density_width = 3;%2; %density kernel in tsne space 0.5 for pca; 2.5 tsne .  !! 5 so far works well
analysisstruct.params.expansion_factor = 1.1; %add a little room to the map after kernel smoothing
analysisstruct.params.density_threshold = 1*10^(-5); %remove regions in plots with low density
analysisstruct.matchedconds = {[unique(analysisstruct.condition_inds)]}; %if running over multiple conditions
analysisstruct.conditions_to_run = [unique(analysisstruct.condition_inds)];
analysisstruct.tsnegranularity = analysisparams.tsnegranularity;

params.reorder=1;
analysisstruct = compute_analysis_clusters_demo(analysisstruct,params);
disp('%% Done clustering %%')

analysisstruct.conditionnames = ratname;
analysisstruct.ratnames = ratname;
analysisstruct.filesizes = {size(mocapstruct.aligned_mean_position,1 );};



%% Run hierarchical clsutering

params.do_show_pdistmatrix =1;
params.decimation_factor = 1; %downsample if needed to save on memory
% params.decimation_factor = 5; %downsample if needed to save on memory
params.doclustering = 1;

%clustering parameters
params.corr_threshold = 0.2;% 0.2
params.clustercutoff =0.65;%0.12; % 0.65
analysisstruct.plotdirectory = '';
%timescale to use, in seconds
params.timescales = [1./4 2];  % this is in minutes, check the find_sequences_state_demo code 
% paper 15 s and 12 s -> 1/4, 2 
% params.timescales = [15, 120]; 
% params.timescales = [0.1];%[1./4 2]; 
% params.timescales = [0.05];%[1./4 2]; 


analysisstruct.conditionnames = {'test'};
analysisstruct.ratname = {ratname};
condition =1;
hierarchystruct=   find_sequences_states_demo(analysisstruct,condition,params);
% identofy the unique cls
[seq_cls, seq_c_idx, seq_r] = unique(hierarchystruct.clustered_behavior{1}, 'stable');

%% Assuming you have a list 'cluster_sequence' that denotes the cluster (state) at each frame
% For example: cluster_sequence = [1, 2, 2, 1, 3, 2, ...];
% choose the timescale (1 or 2)
cluster_sequence = hierarchystruct.clustered_behavior{1};
cluster_sequence = cluster_sequence +1; % nomenclature starts from 0

unique_clusters = unique(cluster_sequence);
n_unique_clusters = length(unique_clusters);
cluster_to_index = containers.Map('KeyType', 'double', 'ValueType', 'double');

for i = 1:n_unique_clusters
    cluster_to_index(unique_clusters(i)) = i;
end

transition_matrix = zeros(n_unique_clusters, n_unique_clusters);

for i = 1:length(cluster_sequence) - 1
    current_cluster = cluster_sequence(i);
    next_cluster = cluster_sequence(i+1);
    
    current_index = cluster_to_index(current_cluster);
    next_index = cluster_to_index(next_cluster);
    
    transition_matrix(current_index, next_index) = transition_matrix(current_index, next_index) + 1;
end

row_sums = sum(transition_matrix, 2);
transition_matrix = bsxfun(@rdivide, transition_matrix, row_sums);


% 
% unique_clusters = unique(cluster_sequence);  % Find the unique clusters (states)
% n_clusters = length(unique_clusters);  % Number of unique clusters
% % n_clusters = 1:max(unique_clusters);
% % Initialize transition matrix
% transition_matrix = zeros(n_clusters, n_clusters);
% 
% % Calculate transition probabilities
% for i = 1:length(cluster_sequence) - 1
% 
%     current_cluster = cluster_sequence(i);
%     next_cluster = cluster_sequence(i+1);
% 
%     transition_matrix(current_cluster, next_cluster) = transition_matrix(current_cluster, next_cluster) + 1;
% end
% 
% 
% 
% % Normalize each row to sum to 1
% row_sums = sum(transition_matrix, 2);
% transition_matrix = bsxfun(@rdivide, transition_matrix, row_sums);

% Display the transition matrix
disp('Transition Matrix:');
% disp(transition_matrix);
figure
% If you want to visualize this matrix
imagesc(transition_matrix);
colorbar;
xlabel('Next Cluster');
ylabel('Current Cluster');
title('Markovian Transition Matrix');
%%
keyboard

% Assume your original transition matrix is called transition_matrix

expected_transition_matrix_15s = transition_matrix ^ 15;
expected_transition_matrix_120s = transition_matrix ^ 120;

observed_transition_matrix_15s = zeros(n_clusters, n_clusters);
observed_transition_matrix_120s = zeros(n_clusters, n_clusters);

% 15s timescale
for i = 1:15:length(cluster_sequence) - 15
    current_cluster = cluster_sequence(i);
    next_cluster = cluster_sequence(i + 15);
    
    observed_transition_matrix_15s(current_cluster, next_cluster) = observed_transition_matrix_15s(current_cluster, next_cluster) + 1;
end

% 120s timescale
for i = 1:120:length(cluster_sequence) - 120
    current_cluster = cluster_sequence(i);
    next_cluster = cluster_sequence(i + 120);
    
    observed_transition_matrix_120s(current_cluster, next_cluster) = observed_transition_matrix_120s(current_cluster, next_cluster) + 1;
end

% Normalize the matrices
row_sums_15s = sum(observed_transition_matrix_15s, 2);
observed_transition_matrix_15s = bsxfun(@rdivide, observed_transition_matrix_15s, row_sums_15s);

row_sums_120s = sum(observed_transition_matrix_120s, 2);
observed_transition_matrix_120s = bsxfun(@rdivide, observed_transition_matrix_120s, row_sums_120s);


