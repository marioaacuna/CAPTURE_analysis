%% script_10_Analysis_miniscope_clusters
% Preamble
% This script will analyse the calcium acitvity of neuronal recordings and
% categorize them into the different clusters, based on clusters ids ,
% This script will read the ROI fluorescence and the cluster vector as inputs.
clear, close all, clc
global GC
% inputs

session_to_take = 2; % this session contains the Sal or PFA data
ROI_traces_path = 'V:\Ca_imaging_pain\4_fluorescence_traces'; % for now we leave them blank, fix later
clusters_path = '';
suffix = '_raw_deltaF_over_F.mat';

cluster_folder = fullfile(GC.preprocessing_rootpath) ;

% Cluster data structure (this contains clusters_struct.(animal_ID_in_struct) = clusters; per animal)
clusters_struct_file = fullfile(cluster_folder, 'clusters_struct.mat');
clusters_struct = load(clusters_struct_file);
clusters_struct = clusters_struct.clusters_struct;
% Animal conditions
animals_of_interest = fieldnames (clusters_struct);
animals_of_interest(ismember(animals_of_interest, 'conditions')) = []; % remove the conditions field
animal_conditions = (clusters_struct.conditions); 

ds_factor = GC.frame_rate/5; % downsampling factor

%%  Pre-allocate structures to store data for each condition
data_H = struct(); 
data_N = struct();
% Loop through animals
for animal = 1:length(animals_of_interest)
    animal_ID = animals_of_interest{animal};
    animal_condition = animal_conditions{animal};
    if sum(~strcmp(animal_condition, {'H', 'N'})) == 2
        continue
    end
    this_cluster_vector = clusters_struct.(animal_ID);

    % Load fluorescence data
    ROI_traces_filename = fullfile(ROI_traces_path,  [animal_ID(1:end-2), suffix]); % animal Id is already JH_XXX
    % load data
    try
    data = load(ROI_traces_filename);
    catch ME
        disp(ME.identifier)
        continue
    end
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

     % Calculate ensemble activity
    ensemble_activity = analyze_neural_ensembles_poses(traces_interpolated, cluster_vector_ds);
    
    
    
    % Store the data depending on the animal condition
    switch animal_condition
        case 'H'
            data_H.(animal_ID).max_amplitude = max_amplitude;
            data_H.(animal_ID).unique_clusters = unique_clusters;
            data_H.(animal_ID).ensemble_activity = ensemble_activity;  % Add ensemble activity
            data_H.(animal_ID).traces = traces_interpolated;  % Store traces for later analysis
            data_H.(animal_ID).cluster_vector = cluster_vector_ds;  % Store cluster vector
        case 'N'
            data_N.(animal_ID).max_amplitude = max_amplitude;
            data_N.(animal_ID).unique_clusters = unique_clusters;
            data_N.(animal_ID).ensemble_activity = ensemble_activity;  % Add ensemble activity
            data_N.(animal_ID).traces = traces_interpolated;  % Store traces for later analysis
            data_N.(animal_ID).cluster_vector = cluster_vector_ds;  % Store cluster vector
        otherwise
            warning(['Unknown condition for animal: ' animal_ID]);
    end

end



%% V2


global_clusters = [];
for animal = 1:length(animals_of_interest)
    animal_ID = animals_of_interest{animal};

    % Check for the animal in 'H' condition
    if isfield(data_H, animal_ID)
        global_clusters = union(global_clusters, data_H.(animal_ID).unique_clusters);
    end

    % Check for the animal in 'N' condition
    if isfield(data_N, animal_ID)
        global_clusters = union(global_clusters, data_N.(animal_ID).unique_clusters);
    end
end


% Concatenate data for each cluster in each group (either H or N)
concatenated_H = cell(1, max(global_clusters) + 1); % +1 since it starts from 0
concatenated_N = cell(1, max(global_clusters) + 1);

for cluster_idx = 0:max(global_clusters) % Starting from 0 as mentioned
    % For 'H' condition
    concatenated_H{cluster_idx + 1} = {}; % Initializing as a cell
    for animal = 1:length(animals_of_interest)
        animal_ID = animals_of_interest{animal};
        if isfield(data_H, animal_ID) % Check if the animal exists in the 'H' structure
            if ismember(cluster_idx, data_H.(animal_ID).unique_clusters) % Check if this cluster exists for the animal
                concatenated_H{cluster_idx + 1}{end+1} = data_H.(animal_ID).max_amplitude(:, data_H.(animal_ID).unique_clusters == cluster_idx); 
            end
        end
    end
    
    % For 'N' condition
    concatenated_N{cluster_idx + 1} = {}; % Initializing as a cell
    for animal = 1:length(animals_of_interest)
        animal_ID = animals_of_interest{animal};
        if isfield(data_N, animal_ID) % Check if the animal exists in the 'N' structure
            if ismember(cluster_idx, data_N.(animal_ID).unique_clusters) % Check if this cluster exists for the animal
                concatenated_N{cluster_idx + 1}{end+1} = data_N.(animal_ID).max_amplitude(:, data_N.(animal_ID).unique_clusters == cluster_idx);
            end
        end
    end
end

% Now, concatenated_H and concatenated_N are cell arrays containing cells for each cluster.
% Each of these inner cells contains matrices from individual animals.
% Step 1: Concatenate all ROIs for each cluster

% Initializing structures to store concatenated data for all ROIs per cluster
all_ROIs_per_cluster_H = struct();
all_ROIs_per_cluster_N = struct();

for cluster_idx = 1:length(global_clusters)
    
    % For 'H' condition
    all_ROIs_H = [];
    if ~isempty(concatenated_H{cluster_idx})
        for animal = 1:length(concatenated_H{cluster_idx})
            all_ROIs_H = [all_ROIs_H; concatenated_H{cluster_idx}{animal}];
        end
    end
    cluster_str = ['cluster_' num2str(global_clusters(cluster_idx))];
    all_ROIs_per_cluster_H.(cluster_str) = all_ROIs_H;
    
    % For 'N' condition
    all_ROIs_N = [];
    if ~isempty(concatenated_N{cluster_idx})
        for animal = 1:length(concatenated_N{cluster_idx})
            all_ROIs_N = [all_ROIs_N; concatenated_N{cluster_idx}{animal}];
        end
    end
    all_ROIs_per_cluster_N.(cluster_str) = all_ROIs_N;
    
end

% Step 2: Compare
% At this point, all_ROIs_per_cluster_H and all_ROIs_per_cluster_N contain all ROI data per cluster
% for H and N conditions respectively. Now, you can run the appropriate statistical tests 
% (e.g., t-tests) to compare all_ROIs_per_cluster_H.(cluster) vs all_ROIs_per_cluster_N.(cluster).
% Initializing structures to store p-values and test statistics

%% stats and plotting
Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]); % Create a rectangular figure

% Define some constants for plotting aesthetics
barWidth = 0.75; % Width of the bars
gapWidth = 1; % Gap between each cluster
currentX = 1; % Starting x position for the first bar


p_values = struct();
t_statistics = struct();

clusters = fieldnames(all_ROIs_per_cluster_H); % Extracting cluster names

for i = 1:length(clusters)
    cluster_name = clusters{i};
    
    data_H = all_ROIs_per_cluster_H.(cluster_name);
    data_N = all_ROIs_per_cluster_N.(cluster_name);
    
    % Perform ttest2
    try
    [h, p, ci, stats] = ttest2(data_H, data_N);
    catch
        p = 0;
        stats.tstat = [];
    end
    % Store p-value and t-statistic
    p_values.(cluster_name) = p;
    t_statistics.(cluster_name) = stats.tstat;


    % Calculate mean and SEM for H
    mean_H = mean(data_H);
    SEM_H = std(data_H)/sqrt(length(data_H));
    
    % Calculate mean and SEM for N
    mean_N = mean(data_N);
    SEM_N = std(data_N)/sqrt(length(data_N));
    
    % Plot bar for H
    bar(currentX, mean_H, barWidth, 'b'); 
    hold on;
    % Plot error bar for H
    errorbar(currentX, mean_H, SEM_H, 'k', 'LineStyle', 'none');
    
    % Plot bar for N right next to H
    bar(currentX + barWidth, mean_N, barWidth, 'r'); 
    % Plot error bar for N
    errorbar(currentX + barWidth, mean_N, SEM_N, 'k', 'LineStyle', 'none');
    try
        if p < 0.05 && p > 0.01
            text(currentX , max(mean_N, mean_H) + max(mean_N, mean_H) * 0.20, '*')
        elseif p < 0.01 && p > 0.001
            text(currentX - barWidth, max(mean_N, mean_H) + max(mean_N, mean_H) * 0.20, '**')
        elseif p < 0.001 && p > 0
            text(currentX - barWidth, max(mean_N, mean_H) + max(mean_N, mean_H) * 0.20, '***')
        end
    end
        

    % Update x position for the next cluster
    currentX = currentX + barWidth * 2 + gapWidth;
end


% Some aesthetic properties for the plot
% legend('Condition H', 'Condition N');

% Create dummy bars for legend
hH = bar(NaN, NaN, 'b'); % dummy bar for H
hN = bar(NaN, NaN, 'r'); % dummy bar for N

legend([hH, hN], {'Sham', 'Neuropathic'});


ylabel('Mean Amplitude');
xlabel('Clusters');
title('Mean Amplitude of ROIs per Cluster for Conditions H and N');
set(gca, 'XTick', 1.5:barWidth*2+gapWidth:length(clusters)*(barWidth*2+gapWidth));
set(gca, 'XTickLabel', clusters, 'TickLabelInterpreter', 'none');
box off
set(gca, 'TickDir', 'out');


%% save figure
keyboard



%% Difference plot
%% Difference Plot
Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]);

barWidth = 0.85; % Wider bars since we only have one per cluster
currentX = 1;
clusters = fieldnames(all_ROIs_per_cluster_H);

% Preallocate arrays for differences and errors
differences = zeros(1, length(clusters));
errors = zeros(1, length(clusters));
p_values_array = zeros(1, length(clusters));

for i = 1:length(clusters)
    cluster_name = clusters{i};
    
    data_H = all_ROIs_per_cluster_H.(cluster_name);
    data_N = all_ROIs_per_cluster_N.(cluster_name);
    
    % Calculate difference and pooled SEM
    mean_diff = mean(data_N) - mean(data_H);
    pooled_SEM = sqrt((std(data_N)^2/length(data_N)) + (std(data_H)^2/length(data_H)));
    
    differences(i) = mean_diff;
    errors(i) = pooled_SEM;
    
    % Store p-value
    try
        [~, p] = ttest2(data_H, data_N);
        p_values_array(i) = p;
    catch
        p_values_array(i) = 1;
    end
end

% Create color array based on p-values
colors = zeros(length(clusters), 3);
colors(p_values_array < 0.05 & p_values_array >= 0.01, :) = repmat([0.8 0.4 0], sum(p_values_array < 0.05 & p_values_array >= 0.01), 1);
colors(p_values_array < 0.01 & p_values_array >= 0.001, :) = repmat([0.9 0.2 0], sum(p_values_array < 0.01 & p_values_array >= 0.001), 1);
colors(p_values_array < 0.001, :) = repmat([1 0 0], sum(p_values_array < 0.001), 1);
colors(p_values_array >= 0.05, :) = repmat([0.7 0.7 0.7], sum(p_values_array >= 0.05), 1);

% Plot
b = bar(differences, barWidth);
b.FaceColor = 'flat';
b.CData = colors;

% Add error bars
hold on;
errorbar(1:length(clusters), differences, errors, 'k.', 'LineStyle', 'none');

% Customize plot
ylabel('Difference in Mean Amplitude (Neuropathic - Sham)');
xlabel('Clusters');
title('Difference in Neuronal Activity between Neuropathic and Sham Conditions');
set(gca, 'XTick', 1:length(clusters));
set(gca, 'XTickLabel', clusters, 'TickLabelInterpreter', 'none');
box off;
set(gca, 'TickDir', 'out');

% Add zero line
yline(0, 'k--', 'Alpha', 0.3);

% Add legend
legend_elements = [
    patch([0 0], [0 0], [0.7 0.7 0.7], 'DisplayName', 'n.s.');
    patch([0 0], [0 0], [0.8 0.4 0], 'DisplayName', 'p < 0.05');
    patch([0 0], [0 0], [0.9 0.2 0], 'DisplayName', 'p < 0.01');
    patch([0 0], [0 0], [1 0 0], 'DisplayName', 'p < 0.001')
];
legend(legend_elements, 'Location', 'northeast');

%% heatmap plot
%% Heat Map
Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]);

clusters = fieldnames(all_ROIs_per_cluster_H);
num_clusters = length(clusters);

% Preallocate matrices for means and p-values
means_matrix = zeros(2, num_clusters);
p_values_array = zeros(1, num_clusters);

% Calculate means and p-values
for i = 1:num_clusters
    cluster_name = clusters{i};
    
    data_H = all_ROIs_per_cluster_H.(cluster_name);
    data_N = all_ROIs_per_cluster_N.(cluster_name);
    
    means_matrix(1,i) = mean(data_H);
    means_matrix(2,i) = mean(data_N);
    
    try
        [~, p] = ttest2(data_H, data_N);
        p_values_array(i) = p;
    catch
        p_values_array(i) = 1;
    end
end

% Create heatmap
imagesc(means_matrix);

% Custom colormap (blue to red)
colormap('jet');
colorbar;

% Add significance markers
hold on;
for i = 1:num_clusters
    if p_values_array(i) < 0.001
        text(i, 1.5, '***', 'HorizontalAlignment', 'center', 'Color', 'k');
    elseif p_values_array(i) < 0.01
        text(i, 1.5, '**', 'HorizontalAlignment', 'center', 'Color', 'k');
    elseif p_values_array(i) < 0.05
        text(i, 1.5, '*', 'HorizontalAlignment', 'center', 'Color', 'k');
    end
end

% Customize plot
ylabel('Condition');
xlabel('Clusters');
title('Mean Amplitude Heatmap of ROIs per Cluster');
set(gca, 'YTick', 1:2);
set(gca, 'YTickLabel', {'Sham', 'Neuropathic'});
set(gca, 'XTick', 1:num_clusters);
set(gca, 'XTickLabel', clusters, 'TickLabelInterpreter', 'none');

%% %%%% Scatter Plot with Connected Pairs
Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]);

clusters = fieldnames(all_ROIs_per_cluster_H);
num_clusters = length(clusters);

% Preallocate arrays
means_H = zeros(1, num_clusters);
means_N = zeros(1, num_clusters);
p_values_array = zeros(1, num_clusters);

% Calculate means and p-values
for i = 1:num_clusters
    cluster_name = clusters{i};
    
    data_H = all_ROIs_per_cluster_H.(cluster_name);
    data_N = all_ROIs_per_cluster_N.(cluster_name);
    
    means_H(i) = mean(data_H);
    means_N(i) = mean(data_N);
    
    try
        [~, p] = ttest2(data_H, data_N);
        p_values_array(i) = p;
    catch
        p_values_array(i) = 1;
    end
end

% Create color array based on p-values
colors = zeros(num_clusters, 3);
colors(p_values_array < 0.05 & p_values_array >= 0.01, :) = repmat([0.8 0.4 0], sum(p_values_array < 0.05 & p_values_array >= 0.01), 1);
colors(p_values_array < 0.01 & p_values_array >= 0.001, :) = repmat([0.9 0.2 0], sum(p_values_array < 0.01 & p_values_array >= 0.001), 1);
colors(p_values_array < 0.001, :) = repmat([1 0 0], sum(p_values_array < 0.001), 1);
colors(p_values_array >= 0.05, :) = repmat([0.7 0.7 0.7], sum(p_values_array >= 0.05), 1);

% Plot connecting lines
for i = 1:num_clusters
    line([means_H(i) means_N(i)], [i i], 'Color', colors(i,:));
end

% Plot scatter points
hold on;
scatter(means_H, 1:num_clusters, 50, 'b', 'filled');
scatter(means_N, 1:num_clusters, 50, 'r', 'filled');

% Customize plot
ylabel('Clusters');
xlabel('Mean Amplitude');
title('Comparison of Mean Amplitude between Conditions');
set(gca, 'YTick', 1:num_clusters);
set(gca, 'YTickLabel', clusters, 'TickLabelInterpreter', 'none');
box off;
set(gca, 'TickDir', 'out');

% Add legend
legend_elements = [
    scatter(NaN, NaN, 50, 'b', 'filled', 'DisplayName', 'Sham');
    scatter(NaN, NaN, 50, 'r', 'filled', 'DisplayName', 'Neuropathic');
    line([NaN NaN], [NaN NaN], 'Color', [0.7 0.7 0.7], 'DisplayName', 'n.s.');
    line([NaN NaN], [NaN NaN], 'Color', [0.8 0.4 0], 'DisplayName', 'p < 0.05');
    line([NaN NaN], [NaN NaN], 'Color', [0.9 0.2 0], 'DisplayName', 'p < 0.01');
    line([NaN NaN], [NaN NaN], 'Color', [1 0 0], 'DisplayName', 'p < 0.001')
];
legend(legend_elements, 'Location', 'northeast');
%% analysis neuronal ensembles pca
script_03__01_01_test_analysis_neuronalensembles()




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions
% Fluorescece   data is sampled at 5Hz
% Cluster data is sampled at 120Hz
% Create a function for downsample vector clusters to match the fluorescence data

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