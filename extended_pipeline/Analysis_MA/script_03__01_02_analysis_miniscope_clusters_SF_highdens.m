%% script_03__01_analysis_miniscope_clusters_SF
% Preamble
% Similar to the original script but analyzing S and F conditions
% S condition uses session 2, F condition uses session 3
clear, close all, clc
global GC

% Modified inputs for S/F analysis
ROI_traces_path = GC.traces_folder;
clusters_path = '';
suffix = '_raw_deltaF_over_F.mat';

cluster_folder = fullfile(GC.preprocessing_rootpath);

% Cluster data structure
clusters_struct_file = fullfile(cluster_folder, 'clusters_struct_high_density.mat');
clusters_struct = load(clusters_struct_file);
clusters_struct = clusters_struct.clusters_struct;

% Animal conditions
animals_of_interest = fieldnames(clusters_struct);
animals_of_interest(ismember(animals_of_interest, 'conditions')) = [];
animal_conditions = (clusters_struct.conditions);

ds_factor = GC.frame_rate/5; % downsampling factor

%%  Pre-allocate structures to store data for each condition
data_S = struct();
data_F = struct();

% Loop through animals
for animal = 1:length(animals_of_interest)
    animal_ID = animals_of_interest{animal};
    animal_condition = animal_conditions{animal};
    
    % Skip if not S or F condition
    if sum(~strcmp(animal_condition, {'S', 'F'})) == 2
        continue
    end
    
    this_cluster_vector = clusters_struct.(animal_ID);

    % Load fluorescence data
    ROI_traces_filename = fullfile(ROI_traces_path, [animal_ID(1:end-2), suffix]);
    
    try
        data = load(ROI_traces_filename);
    catch ME
        disp(ME.identifier)
        continue
    end


    % Animal: ID_1386 has only 2 sessions, so select 1 and 2
    % Select appropriate session based on condition
    if strcmp(animal_condition, 'S')
        session_to_use = 2;
    else % F condition
        session_to_use = 3;
    end
    
    if startsWith(animal_ID, 'ID_1386')
        session_to_use = session_to_use - 1;
    end
    traces = data.dFF(:,session_to_use);
    traces = cell2mat(traces);

    % Process clusters and traces
    cluster_vector_ds = downsample_vector(this_cluster_vector, (traces));
    
    % Ensure lengths are consistent after downsampling
    min_length = min(length(cluster_vector_ds), length(traces));
    cluster_vector_ds = cluster_vector_ds(1:min_length);
    traces_interpolated = traces(:, 1:min_length);

    % Calculate metrics
    [max_amplitude, unique_clusters] = calculate_max_amplitude(traces_interpolated, cluster_vector_ds);
    ensemble_activity = analyze_neural_ensembles_poses(traces_interpolated, cluster_vector_ds);
    
    % Store the data depending on the animal condition
    switch animal_condition
        case 'S'
            data_S.(animal_ID).max_amplitude = max_amplitude;
            data_S.(animal_ID).unique_clusters = unique_clusters;
            data_S.(animal_ID).ensemble_activity = ensemble_activity;
            data_S.(animal_ID).traces = traces_interpolated;
            data_S.(animal_ID).cluster_vector = cluster_vector_ds;
        case 'F'
            data_F.(animal_ID).max_amplitude = max_amplitude;
            data_F.(animal_ID).unique_clusters = unique_clusters;
            data_F.(animal_ID).ensemble_activity = ensemble_activity;
            data_F.(animal_ID).traces = traces_interpolated;
            data_F.(animal_ID).cluster_vector = cluster_vector_ds;
        otherwise
            warning(['Unknown condition for animal: ' animal_ID]);
    end
end

%% 
% Call the function to analyze and plot the calcium metrics comparison
[A,Aor] = analyze_calcium_metrics_comparison(data_S, data_F, animals_of_interest, 'max_amplitude', 'S_v_F');
[P, Por] = analyze_calcium_metrics_comparison(data_S, data_F, animals_of_interest, 'peaks', 'S_v_F');
[F,Ford] = analyze_calcium_metrics_comparison(data_S, data_F, animals_of_interest, 'freqs', 'S_v_F');


%% Plot poses that are gained and lost in activity
% for now let's take only one. Peak amplitude (P)
% load analysis struct
logger('Loading analysisstrcut', 'INFO');
% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');
%%




cls = F.increased;
% cls = [5,6,2];
plot_poses = 1;
if plot_poses
    % h= figure(370);
    % clf;

    fig_i = figure('pos', [10,300,1500,1900]);
    nclus = numel(cls);
    n_rows = ceil(sqrt(nclus));
    n_cols = ceil(sqrt(nclus));
    for ic = 1:numel(cls)
        subplot(n_rows, n_cols, ic)
        this_cls = cls(ic);
        fprintf('ic = %i - \n', this_cls)
        plot_mean_cluster_aligned(analysisstruct.highdensity_analysisstruct.mocapstruct_reduced_agg{1},...
            find(analysisstruct.highdensity_analysisstruct.annot_reordered{end}==this_cls),['cl nr :  ', num2str(this_cls)]);
        title(this_cls)
    end
end


cls = F.decreased;
plot_poses = 1;
if plot_poses
    % h= figure(370);
    % clf;

    fig_d = figure('pos', [10,300,1500,1900]);
    nclus = numel(cls);
    n_rows = ceil(sqrt(nclus));
    n_cols = ceil(sqrt(nclus));
    for ic = 1:numel(cls)
        subplot(n_rows, n_cols, ic)
        this_cls = cls(ic);
        fprintf('ic = %i - \n', this_cls)
        plot_mean_cluster_aligned(analysisstruct.highdensity_analysisstruct.mocapstruct_reduced_agg{1},...
            find(analysisstruct.highdensity_analysisstruct.annot_reordered{end}==this_cls),['cl nr :  ', num2str(this_cls)]);
        title(this_cls)
    end
end

%% analysis neuronal ensembles pca
% names = {'Saline', 'Formalin'};
% script_03__01_01_test_analysis_neuronalensembles(data_S, data_F, names)