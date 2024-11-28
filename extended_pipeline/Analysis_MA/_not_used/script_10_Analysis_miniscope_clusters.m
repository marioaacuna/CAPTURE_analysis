%% script_10_Analysis_miniscope_clusters
% Preamble
% This script will analyse the calcium acitvity of neuronal recordings and
% categorize them into the different clusters, based on clusters ids ,
% This script will read the ROI fluorescence and the cluster vector as inputs.
clear, close all, clc
% inputs

session_to_take = 2; % this session contains the Sal or PFA data
ROI_traces_path = 'V:\Ca_imaging_pain\4_fluorescence_traces'; % for now we leave them blank, fix later
clusters_path = '';
suffix = '_raw_deltaF_over_F.mat';

cluster_folder = "H:\Mario\Results_Capture\clusters";

% Cluster data structure (this contains clusters_struct.(animal_ID_in_struct) = clusters; per animal)
clusters_struct_file = fullfile(cluster_folder, 'clusters_struct.mat');
clusters_struct = load(clusters_struct_file);
clusters_struct = clusters_struct.clusters_struct;
% Animal conditions
animals_of_interest = fieldnames (clusters_struct);
animals_of_interest(ismember(animals_of_interest, 'conditions')) = []; % remove the conditions field
animal_conditions = (clusters_struct.conditions); 
% animal_of_interest = {'326', '327', '328', '330', '332', '334', '335', '336'};
%conditions = {'F', 'F', 'S', 'S', 'S', 'F', 'F', 'F'};
ds_factor = 120/5; % downsampling factor


% Pre-allocate structures to store data for each condition
data_S = struct(); 
data_F = struct();
% Loop through animals
for animal = 1:length(animals_of_interest)
    animal_ID = animals_of_interest{animal};
    animal_condition = animal_conditions{animal};
    this_cluster_vector = clusters_struct.(animal_ID);

    % Load fluorescence data
    ROI_traces_filename = fullfile(ROI_traces_path,  [animal_ID, suffix]); % animal Id is already JH_XXX
    % load data
    data = load(ROI_traces_filename);
    traces = data.dFF(:,session_to_take); % Assuming that the fluorescence data is stored in a cell named 'session_data'
    traces = cell2mat(traces);

    % downsample cluster vector to match fluorescence sampling rate
    cluster_vector_ds = downsample_vector(this_cluster_vector, (traces));
    
    % Ensure lengths are consistent after downsampling
    min_length = min(length(cluster_vector_ds), length(traces));
    cluster_vector_ds = cluster_vector_ds(1:min_length);
    traces_interpolated = traces( :, 1:min_length);

    % calculate the max amplitude per cluster
    [max_amplitude, unique_clusters] = calculate_max_amplitude(traces_interpolated, cluster_vector_ds);
    
    % Store the data depending on the animal condition (either S or F)
    switch animal_condition
        case 'S'
            data_S.(animal_ID).max_amplitude = max_amplitude;
            data_S.(animal_ID).unique_clusters = unique_clusters;
        case 'F'
            data_F.(animal_ID).max_amplitude = max_amplitude;
            data_F.(animal_ID).unique_clusters = unique_clusters;
        otherwise
            warning(['Unknown condition for animal: ' animal_ID]);
    end

end



%% V2


global_clusters = [];
for animal = 1:length(animals_of_interest)
    animal_ID = animals_of_interest{animal};

    % Check for the animal in 'S' condition
    if isfield(data_S, animal_ID)
        global_clusters = union(global_clusters, data_S.(animal_ID).unique_clusters);
    end

    % Check for the animal in 'F' condition
    if isfield(data_F, animal_ID)
        global_clusters = union(global_clusters, data_F.(animal_ID).unique_clusters);
    end
end


% Concatenate data for each cluster in each group (either S or F)
concatenated_S = cell(1, max(global_clusters) + 1); % +1 since it starts from 0
concatenated_F = cell(1, max(global_clusters) + 1);

for cluster_idx = 0:max(global_clusters) % Starting from 0 as mentioned
    % For 'S' condition
    concatenated_S{cluster_idx + 1} = {}; % Initializing as a cell
    for animal = 1:length(animals_of_interest)
        animal_ID = animals_of_interest{animal};
        if isfield(data_S, animal_ID) % Check if the animal exists in the 'S' structure
            if ismember(cluster_idx, data_S.(animal_ID).unique_clusters) % Check if this cluster exists for the animal
                concatenated_S{cluster_idx + 1}{end+1} = data_S.(animal_ID).max_amplitude(:, data_S.(animal_ID).unique_clusters == cluster_idx); 
            end
        end
    end
    
    % For 'F' condition
    concatenated_F{cluster_idx + 1} = {}; % Initializing as a cell
    for animal = 1:length(animals_of_interest)
        animal_ID = animals_of_interest{animal};
        if isfield(data_F, animal_ID) % Check if the animal exists in the 'F' structure
            if ismember(cluster_idx, data_F.(animal_ID).unique_clusters) % Check if this cluster exists for the animal
                concatenated_F{cluster_idx + 1}{end+1} = data_F.(animal_ID).max_amplitude(:, data_F.(animal_ID).unique_clusters == cluster_idx);
            end
        end
    end
end

% Now, concatenated_S and concatenated_F are cell arrays containing cells for each cluster.
% Each of these inner cells contains matrices from individual animals.
% Step 1: Concatenate all ROIs for each cluster

% Initializing structures to store concatenated data for all ROIs per cluster
all_ROIs_per_cluster_S = struct();
all_ROIs_per_cluster_F = struct();

for cluster_idx = 1:length(global_clusters)
    
    % For 'S' condition
    all_ROIs_S = [];
    if ~isempty(concatenated_S{cluster_idx})
        for animal = 1:length(concatenated_S{cluster_idx})
            all_ROIs_S = [all_ROIs_S; concatenated_S{cluster_idx}{animal}];
        end
    end
    cluster_str = ['cluster_' num2str(global_clusters(cluster_idx))];
    all_ROIs_per_cluster_S.(cluster_str) = all_ROIs_S;
    
    % For 'F' condition
    all_ROIs_F = [];
    if ~isempty(concatenated_F{cluster_idx})
        for animal = 1:length(concatenated_F{cluster_idx})
            all_ROIs_F = [all_ROIs_F; concatenated_F{cluster_idx}{animal}];
        end
    end
    all_ROIs_per_cluster_F.(cluster_str) = all_ROIs_F;
    
end

% Step 2: Compare
% At this point, all_ROIs_per_cluster_S and all_ROIs_per_cluster_F contain all ROI data per cluster
% for S and F conditions respectively. Now, you can run the appropriate statistical tests 
% (e.g., t-tests) to compare all_ROIs_per_cluster_S.(cluster) vs all_ROIs_per_cluster_F.(cluster).
% Initializing structures to store p-values and test statistics

%% stats and plotting
Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]); % Create a rectangular figure

% Define some constants for plotting aesthetics
barWidth = 0.75; % Width of the bars
gapWidth = 1; % Gap between each cluster
currentX = 1; % Starting x position for the first bar


p_values = struct();
t_statistics = struct();

clusters = fieldnames(all_ROIs_per_cluster_S); % Extracting cluster names

for i = 1:length(clusters)
    cluster_name = clusters{i};
    
    data_S = all_ROIs_per_cluster_S.(cluster_name);
    data_F = all_ROIs_per_cluster_F.(cluster_name);
    
    % Perform ttest2
    try
    [h, p, ci, stats] = ttest2(data_S, data_F);
    catch
        p = 0;
        stats.tstat = [];
    end
    % Store p-value and t-statistic
    p_values.(cluster_name) = p;
    t_statistics.(cluster_name) = stats.tstat;


    % Calculate mean and SEM for S
    mean_S = mean(data_S);
    SEM_S = std(data_S)/sqrt(length(data_S));
    
    % Calculate mean and SEM for F
    mean_F = mean(data_F);
    SEM_F = std(data_F)/sqrt(length(data_F));
    
    % Plot bar for S
    bar(currentX, mean_S, barWidth, 'b'); 
    hold on;
    % Plot error bar for S
    errorbar(currentX, mean_S, SEM_S, 'k', 'LineStyle', 'none');
    
    % Plot bar for F right next to S
    bar(currentX + barWidth, mean_F, barWidth, 'r'); 
    % Plot error bar for F
    errorbar(currentX + barWidth, mean_F, SEM_F, 'k', 'LineStyle', 'none');
    
    if p < 0.05 && p > 0.01
        text(currentX , max(mean_F, mean_S) + max(mean_F, mean_S) * 0.20, '*')
    elseif p < 0.01 && p > 0.001
        text(currentX - barWidth, max(mean_F, mean_S) + max(mean_F, mean_S) * 0.20, '**')
    elseif p < 0.001 && p > 0
        text(currentX - barWidth, max(mean_F, mean_S) + max(mean_F, mean_S) * 0.20, '***')
    end
        

    % Update x position for the next cluster
    currentX = currentX + barWidth * 2 + gapWidth;
end


% Some aesthetic properties for the plot
% legend('Condition S', 'Condition F');

% Create dummy bars for legend
hS = bar(NaN, NaN, 'b'); % dummy bar for S
hF = bar(NaN, NaN, 'r'); % dummy bar for F

legend([hS, hF], {'Saline', 'Formalin'});


ylabel('Mean Amplitude');
xlabel('Clusters');
title('Mean Amplitude of ROIs per Cluster for Conditions S and F');
set(gca, 'XTick', 1.5:barWidth*2+gapWidth:length(clusters)*(barWidth*2+gapWidth));
set(gca, 'XTickLabel', clusters, 'TickLabelInterpreter', 'none');
box off
set(gca, 'TickDir', 'out');


%% save figure
keyboard



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions
% Fluorescece   data is sampled at 5Hz
% Cluster data is sampled at 120Hz
% Create a function for downsample vector clusters to match the fluorescence data
function cluster_vector_ds = downsample_vector(cluster_vector, traces)
    %% Downsample cluster vector to match fluorescence sampling rate
    % Considering missing frames in fluorescence, let's linearly interpolate
    %missing_frames_ratio = length(cluster_vector)/length_traces;
    %ds_factor_adjusted = ds_factor * missing_frames_ratio;
    %cluster_vector_ds = downsample(cluster_vector, round(ds_factor));

    % %%
    % Linearly interpolate the cluster_vector to match the fluorescence traces length
    cluster_vector_interpolated = interp1(1:length(cluster_vector), cluster_vector, linspace(1, length(cluster_vector), length(traces)), 'nearest');


    % cluster_vector_ds = downsample(cluster_vector, ds_factor);
    % Ensure lengths are consistent after downsampling
    min_length = min(length(cluster_vector_interpolated), length(traces));
    cluster_vector_ds = cluster_vector_interpolated(1:min_length);
end

function [max_amplitudes, unique_clusters] = calculate_max_amplitude(traces_interpolated, cluster_vector_ds)
    %% Collect max amplitude of all cells in this session for all clusters
    unique_clusters = unique(cluster_vector_ds);
    n_clusters = length(unique_clusters);
    n_rois = size(traces_interpolated, 1);
    max_amplitudes = zeros(n_rois, n_clusters);

    for i = 1:n_clusters
        current_cluster = unique_clusters(i);
        idx = cluster_vector_ds == current_cluster; % indices where current cluster is present
        for j = 1:n_rois
            roi_trace = traces_interpolated(j, idx);
            max_roi = max(roi_trace);
            if max_roi<0
                max_roi = 0;
            end
            max_amplitudes(j, i) =max_roi;
        end
    end
end