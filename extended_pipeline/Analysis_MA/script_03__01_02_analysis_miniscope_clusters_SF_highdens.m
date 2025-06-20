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




cls = A.increased;
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


cls = A.decreased;
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


%%



%% analysis neuronal ensembles pca
% names = {'Saline', 'Formalin'};
% script_03__01_01_test_analysis_neuronalensembles(data_S, data_F, names)

%% Plot t-SNE map with increased and decreased clusters (generic)
logger('Creating t-SNE visualization for cluster analysis', 'INFO');

% Create t-SNE difference visualization for A analysis (max amplitude)
plot_tsne_clusters_SF_difference(analysisstruct, A, 'max_amplitude');

% Uncomment below to plot other analyses
% plot_tsne_clusters_SF_difference(analysisstruct, P, 'peaks');
% plot_tsne_clusters_SF_difference(analysisstruct, F, 'frequency');

%% Function to plot t-SNE map with cluster differences
function plot_tsne_clusters_SF_difference(analysisstruct, data_struct, analysis_type)
    % Generic function to plot t-SNE cluster differences
    % Inputs:
    %   - analysisstruct: analysis structure containing t-SNE data
    %   - data_struct: structure containing .increased and .decreased fields
    %   - analysis_type: string describing the analysis ('max_amplitude', 'peaks', 'frequency')
    
    % Define analysis type mapping for titles and filenames
    analysis_mapping = containers.Map(...
        {'max_amplitude', 'peaks', 'frequency'}, ...
        {'Max Amplitude', 'Peak Count', 'Frequency'});
    
    % Get display name for the analysis type
    if isKey(analysis_mapping, analysis_type)
        display_name = analysis_mapping(analysis_type);
    else
        display_name = strrep(analysis_type, '_', ' ');
        display_name = [upper(display_name(1)), display_name(2:end)]; % Capitalize first letter
    end
    % Create figure with cluster differences visualization
    fig = figure('Position', [100, 100, 1400, 600]);
    set(fig, 'Color', 'w');
    
    % Define colors for different cluster states
    increased_color = [1, 0, 0]; % Red for increased
    decreased_color = [0, 0, 1]; % Blue for decreased
    neutral_color = [0.7, 0.7, 0.7]; % Gray for unchanged
    
    % Get the high-density analysis data
    if isfield(analysisstruct, 'highdensity_analysisstruct')
        zValues = analysisstruct.highdensity_analysisstruct.zValues;
        cluster_assignments = analysisstruct.highdensity_analysisstruct.annot_reordered{end};
    else
        zValues = analysisstruct.zValues;
        cluster_assignments = analysisstruct.annot_reordered{end};
    end
      % Create masks for different cluster types
    increased_mask = false(size(zValues, 1), 1);
    decreased_mask = false(size(zValues, 1), 1);
    
    % Create masks for increased clusters
    for i = 1:length(data_struct.increased)
        cluster_id = data_struct.increased(i);
        increased_mask = increased_mask | (cluster_assignments == cluster_id)';
    end
    
    % Create masks for decreased clusters
    for i = 1:length(data_struct.decreased)
        cluster_id = data_struct.decreased(i);
        decreased_mask = decreased_mask | (cluster_assignments == cluster_id)';
    end
    
    % % Plot with cluster boundaries if available
    % subplot(1, 2, 1);
    % 
    % % Plot all points first (neutral/unchanged)
    % plot(zValues(:,1), zValues(:,2), '.', ...
    %     'Color', neutral_color, 'MarkerSize', 1);
    % hold on;
    % 
    % % Overlay increased clusters (red)
    % if sum(increased_mask) > 0
    %     plot(zValues(increased_mask,1), zValues(increased_mask,2), '.', ...
    %         'Color', increased_color, 'MarkerSize', 3);
    % end
    % 
    % % Overlay decreased clusters (blue)
    % if sum(decreased_mask) > 0
    %     plot(zValues(decreased_mask,1), zValues(decreased_mask,2), '.', ...
    %         'Color', decreased_color, 'MarkerSize', 3);
    % end
    % 
    % % Add watershed boundaries if available
    % if isfield(analysisstruct, 'highdensity_analysisstruct')
    %     analysis_struct_to_use = analysisstruct.highdensity_analysisstruct;
    % else
    %     analysis_struct_to_use = analysisstruct;
    % end
    % 
    % if isfield(analysis_struct_to_use, 'sorted_watershed') && ...
    %    isfield(analysis_struct_to_use, 'xx') && ...
    %    isfield(analysis_struct_to_use, 'yy')
    %     nnn = analysis_struct_to_use.sorted_watershed;
    %     nnn(nnn > 0) = 1;
    %     B = bwboundaries(nnn);
    % 
    %     for kk = 1:numel(B)
    %         if size(B{kk}, 1) > 0
    %             plot(analysis_struct_to_use.xx(B{kk}(:,2)), ...
    %                  analysis_struct_to_use.yy(B{kk}(:,1)), ...
    %                  'k-', 'LineWidth', 0.5);
    %         end
    %     end
    % end
    % 
    % title(sprintf('S vs F Cluster Differences (%s Analysis)', display_name));
    % xlabel('t-SNE 1');
    % ylabel('t-SNE 2');
    % axis equal;
    % legend({'Unchanged', 'Increased in F', 'Decreased in F', 'Cluster Boundaries'}, ...
    %        'Location', 'best');
    
    % Create separate plots for increased and decreased clusters
    % subplot(1, 2, 2);
    
    % Plot base map
    plot(zValues(:,1), zValues(:,2), '.', ...
        'Color', [0.9, 0.9, 0.9], 'MarkerSize', 1);
    hold on;
    
    % Plot only the significantly changed clusters
    if sum(increased_mask) > 0
        plot(zValues(increased_mask,1), zValues(increased_mask,2), '.', ...
            'Color', increased_color, 'MarkerSize', 4);
    end
    
    if sum(decreased_mask) > 0
        plot(zValues(decreased_mask,1), zValues(decreased_mask,2), '.', ...
            'Color', decreased_color, 'MarkerSize', 4);
    end
    
    title('Significantly Changed Clusters Only');
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    axis square;
      % Add cluster ID annotations for the changed clusters
    if sum(increased_mask) > 0 || sum(decreased_mask) > 0        % Label increased clusters
        for i = 1:length(data_struct.increased)
            cluster_id = data_struct.increased(i);
            cluster_points = cluster_assignments == cluster_id;
            if sum(cluster_points) > 0
                centroid_x = mean(zValues(cluster_points, 1));
                centroid_y = mean(zValues(cluster_points, 2));
                text(centroid_x, centroid_y, num2str(cluster_id), ...
                     'FontSize', 8, 'FontWeight', 'bold', ...
                     'Color', increased_color, 'HorizontalAlignment', 'center');
            end
        end
        
        % Label decreased clusters
        for i = 1:length(data_struct.decreased)
            cluster_id = data_struct.decreased(i);
            cluster_points = cluster_assignments == cluster_id;
            if sum(cluster_points) > 0
                centroid_x = mean(zValues(cluster_points, 1));
                centroid_y = mean(zValues(cluster_points, 2));
                text(centroid_x, centroid_y, num2str(cluster_id), ...
                     'FontSize', 8, 'FontWeight', 'bold', ...
                     'Color', decreased_color, 'HorizontalAlignment', 'center');
            end
        end
    end
    
    legend({'Background', 'Increased in F', 'Decreased in F'}, 'Location', 'best');    % Add main title with summary statistics
    n_increased = length(data_struct.increased);
    n_decreased = length(data_struct.decreased);
    sgtitle(sprintf('%s Analysis: %d Increased, %d Decreased Clusters', ...
            display_name, n_increased, n_decreased), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    % Log cluster information
    fprintf('t-SNE Cluster Analysis Summary (%s):\n', display_name);
    fprintf('Increased clusters (n=%d): %s\n', n_increased, mat2str(data_struct.increased));
    fprintf('Decreased clusters (n=%d): %s\n', n_decreased, mat2str(data_struct.decreased));    % Save figure using exportgraphics as PDF
    global GC
    if isfield(GC, 'temp_root')
        export_folder = fullfile(GC.temp_root, 'figs_presentation_painAI');
        if exist(export_folder, 'dir')
            fig_filename = fullfile(export_folder, sprintf('tsne_clusters_SF_%s_S_F.pdf', analysis_type));
            exportgraphics(fig, fig_filename, 'ContentType', 'vector');
            fprintf('Saved t-SNE figure as PDF: %s\n', fig_filename);
        end
    end
end

%%