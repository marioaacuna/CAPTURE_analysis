%% Multi-Experiment Behavioral Cluster Analysis Script
%% Initialization
%logger('Starting multi-experiment behavioral cluster analysis script', 'INFO');
clear;
close all;
clc;
GC = general_configs;

% Configuration for visualization and export
debugging = false;  % Set to true for debugging mode
if debugging
    visualize = 'on';  % Show figures during debugging
    do_export = false;  % Don't export during debugging
else
    visualize = 'off';  % Don't show figures in production mode
    do_export = true;   % Export figures in production mode
end

% Export folder
export_folder = fullfile(GC.temp_root, 'figs_many_zvals');
if ~exist(export_folder, 'dir')
    mkdir(export_folder);
end


% Define, zvals to plot
zvals_to_plot =  {{'zValues_all_features'      }
    {'zValues_important_features'}
    {'zValues_jt_features_only'  }
    {'zValues_hand_made'         }
    {'zValues_only_important'    }
};

% <TODO> modify the code such it iterates throught the zvals


%% Define experiments and their conditions
experiments = struct();
experiments(1).name = 'BSFC_300hz';
experiments(1).conditions = {'B', 'S', 'F', 'C'};
experiments(1).folder = '0_preprocessing_BSFC_300hz';

experiments(2).name = 'BHNG_300hz';
experiments(2).conditions = {'B', 'H', 'N', 'G'};
experiments(2).folder = '0_preprocessing_BHNG_300hz';

% Define colors for all conditions
color_map = containers.Map();
color_map('B') = [0.5, 0.5, 0.5];         % Gray
color_map('S') = [0.4660, 0.6740, 0.1880]; % Green
color_map('F') = [0.8500, 0.3250, 0.0980]; % Red
color_map('C') = [0, 0.4470, 0.7410];     % Blue
color_map('H') = [0.4660, 0.6740, 0.1880]; % Green
color_map('N') = [0.8500, 0.3250, 0.0980]; % Red
color_map('G') = [0, 0.4470, 0.7410];     % Blue

%% Process each experiment
for exp_idx = 1:length(experiments)
    experiment = experiments(exp_idx);
    logger(['Processing experiment: ' experiment.name], 'INFO');
    
    % Define paths for this experiment
    exp_data_folder = fullfile(GC.project_path, 'data', experiment.folder);
    
    % Check if experiment folder exists
    if ~exist(exp_data_folder, 'dir')
        logger(['Warning: Experiment folder does not exist: ' exp_data_folder], 'WARN');
        continue;
    end
    
    % Load data for this experiment
    exp_data = load_experiment_data(exp_data_folder, GC);
    
    if isempty(exp_data)
        logger(['Warning: Could not load data for experiment: ' experiment.name], 'WARN');
        continue;
    end
    
    % Create experiment-specific export folder
    exp_export_folder = fullfile(export_folder, experiment.name);
    if ~exist(exp_export_folder, 'dir')
        mkdir(exp_export_folder);
    end
    
    % Loop over z-values configurations and plot analysis for this experiment
    for zv_idx = 1:numel(zvals_to_plot)
        zv_entry = zvals_to_plot{zv_idx};
        if iscell(zv_entry)
            zv_name = zv_entry{1};
        else
            zv_name = zv_entry;
        end

        % Check that the requested zvals field exists
        if ~isfield(exp_data.analysisstruct, zv_name)
            logger(['Warning: z-values field not found in analysisstruct: ' zv_name], 'WARN');
            continue;
        end

        % Update zvals to the selected field
        exp_data.zvals = exp_data.analysisstruct.(zv_name);

        % Create subfolder for this zvals selection
        zv_export_folder = fullfile(exp_export_folder, zv_name);
        if ~exist(zv_export_folder, 'dir')
            mkdir(zv_export_folder);
        end

        % Plot analysis for this experiment and zvals selection
        plot_experiment_analysis(exp_data, experiment, color_map, zv_export_folder, visualize, do_export);
    end
    
    logger(['Completed analysis for experiment: ' experiment.name], 'INFO');
end

logger('Multi-experiment analysis complete', 'INFO');

%% Helper Functions

function exp_data = load_experiment_data(data_folder, GC)
    % Load data for a specific experiment
    exp_data = struct();
    
    try
        % Define file paths
        filename_analysis = fullfile(data_folder, 'raw_concat_analysis.mat');
        filename_predictions = fullfile(data_folder, 'agg_predictions.mat');
        filename_ratception = fullfile(data_folder, 'ratception_prediction.mat');
        
        % Check if files exist
        if ~exist(filename_analysis, 'file')
            logger(['Analysis file not found: ' filename_analysis], 'ERROR');
            exp_data = [];
            return;
        end
        
        if ~exist(filename_predictions, 'file')
            logger(['Predictions file not found: ' filename_predictions], 'ERROR');
            exp_data = [];
            return;
        end
        
        if ~exist(filename_ratception, 'file')
            logger(['Ratception file not found: ' filename_ratception], 'ERROR');
            exp_data = [];
            return;
        end
        
        % Load data
        logger('Loading analysis structure', 'INFO');
        load(filename_analysis, 'analysisstruct');
        
        logger('Loading predictions', 'INFO');
        load(filename_predictions, 'predictions', 'animal_condition_identifier');
        
        logger('Loading ratception structure', 'INFO');
        load(filename_ratception, 'ratception_struct');
        
        % Process data
        upsampling_factor = GC.repfactor;
        long_animal_frames_identifier = repelem(animal_condition_identifier, upsampling_factor);
        animal_list_used_after_analysis = long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});
        
        % Extract condition identifiers
        conditions = cellfun(@(x) x(end), animal_list_used_after_analysis, 'UniformOutput', false);
        
        % Store in output structure
        exp_data.analysisstruct = analysisstruct;
        exp_data.predictions = predictions;
        exp_data.ratception_struct = ratception_struct;
        exp_data.animal_condition_identifier = animal_condition_identifier;
        exp_data.conditions = conditions;
        exp_data.zvals = analysisstruct.zValues;
        
        logger('Data loaded successfully', 'INFO');
        
    catch ME
        logger(['Error loading experiment data: ' ME.message], 'ERROR');
        exp_data = [];
    end
end

function plot_experiment_analysis(exp_data, experiment, color_map, export_folder, visualize, do_export)
    % Plot analysis for a specific experiment
    
    conditions = exp_data.conditions;
    zvals = exp_data.zvals;
    analysisstruct = exp_data.analysisstruct;
    
    % Define condition comparisons based on experiment
    if strcmp(experiment.name, 'BSFC_300hz')
        comparisons = struct();
        comparisons(1).name = 'SF';
        comparisons(1).conditions = {'S', 'F'};
        comparisons(1).title = 'S vs F';
        
        comparisons(2).name = 'BC';
        comparisons(2).conditions = {'B', 'C'};
        comparisons(2).title = 'B vs C';
        
        comparisons(3).name = 'BS';
        comparisons(3).conditions = {'B', 'S'};
        comparisons(3).title = 'B vs S';
        
        comparisons(4).name = 'FC';
        comparisons(4).conditions = {'F', 'C'};
        comparisons(4).title = 'F vs C';
        
        comparisons(5).name = 'SC';
        comparisons(5).conditions = {'S', 'C'};
        comparisons(5).title = 'S vs C';
        
        comparisons(6).name = 'BSF';
        comparisons(6).conditions = {'B', 'S', 'F'};
        comparisons(6).title = 'B vs S vs F';
        
        comparisons(7).name = 'ALL';
        comparisons(7).conditions = {'B', 'S', 'F', 'C'};
        comparisons(7).title = 'All Conditions (B, S, F, C)';
        
    elseif strcmp(experiment.name, 'BHNG_300hz')
        comparisons = struct();
        comparisons(1).name = 'HN';
        comparisons(1).conditions = {'H', 'N'};
        comparisons(1).title = 'H vs N';
        
        comparisons(2).name = 'BG';
        comparisons(2).conditions = {'B', 'G'};
        comparisons(2).title = 'B vs G';
        
        comparisons(3).name = 'NG';
        comparisons(3).conditions = {'N', 'G'};
        comparisons(3).title = 'N vs G';
        
        comparisons(4).name = 'BHN';
        comparisons(4).conditions = {'B', 'H', 'N'};
        comparisons(4).title = 'B vs H vs N';
        
        comparisons(5).name = 'HNG';
        comparisons(5).conditions = {'H', 'N', 'G'};
        comparisons(5).title = 'H vs N vs G';
        
        comparisons(6).name = 'ALL';
        comparisons(6).conditions = {'B', 'H', 'N', 'G'};
        comparisons(6).title = 'All Conditions (B, H, N, G)';
    else
        logger(['Unknown experiment: ' experiment.name], 'ERROR');
        return;
    end
    
    % Create plots for each comparison
    for comp_idx = 1:length(comparisons)
        comparison = comparisons(comp_idx);
        
        % Create density maps
        create_density_maps(comparison, conditions, zvals, analysisstruct, experiment, ...
                          export_folder, visualize, do_export);
        
        % Create scatter plots
        create_scatter_plots(comparison, conditions, zvals, color_map, experiment, ...
                           export_folder, visualize, do_export);
    end
end

function create_density_maps(comparison, conditions, zvals, analysisstruct, experiment, export_folder, visualize, do_export)
    % Create density maps for a specific comparison
    
    logger(['Creating density maps for ' comparison.title ' (' experiment.name ')'], 'INFO');
    
    n_conditions = length(comparison.conditions);
    
    % Calculate figure width based on number of conditions
    fig_width = 400 * n_conditions;
    
    fig_name = ['Density Maps: ' comparison.title ' (' experiment.name ')'];
    fig = figure('Name', fig_name, 'Color', 'w', 'Position', [100, 100, fig_width, 400], 'Visible', visualize);
    
    for i = 1:n_conditions
        cond = comparison.conditions{i};
        
        % Find indices for this condition
        idx = strcmp(conditions, cond);
        
        if sum(idx) == 0
            logger(['Warning: No data found for condition ' cond], 'WARN');
            continue;
        end
        
        subplot(1, n_conditions, i);
        h = gca;
        set(h, 'Color', 'w');
        
        plotdensitymaps({zvals(idx,:)}, 1, h, analysisstruct.params.density_width, ...
            max(zvals(:))*analysisstruct.params.expansion_factor, ...
            analysisstruct.params.density_res);
        
        title(['Condition ' cond]);
        axis square;
    end
    
    % Export the figure if needed
    if do_export
        filename = sprintf('density_map_%s_%s.pdf', experiment.name, comparison.name);
        filepath = fullfile(export_folder, filename);
        logger(['Exporting density map to: ' filepath], 'INFO');
        exportgraphics(fig, filepath, 'ContentType', 'vector', 'BackgroundColor', 'none');
    end
end

function create_scatter_plots(comparison, conditions, zvals, color_map, experiment, export_folder, visualize, do_export)
    % Create scatter plots for a specific comparison
    
    logger(['Creating scatter plot for ' comparison.title ' (' experiment.name ')'], 'INFO');
    
    fig_name = ['Scatter Plot: ' comparison.title ' (' experiment.name ')'];
    fig = figure('Name', fig_name, 'Color', 'w', 'Visible', visualize);
    
    % Set axes background to white
    ax = gca;
    set(ax, 'Color', 'w');
    
    hold on;
    
    % Plot each condition
    for i = 1:length(comparison.conditions)
        cond = comparison.conditions{i};
        
        % Find indices for this condition
        idx = strcmp(conditions, cond);
        
        if sum(idx) == 0
            logger(['Warning: No data found for condition ' cond], 'WARN');
            continue;
        end
        
        % Get color for this condition
        if isKey(color_map, cond)
            color = color_map(cond);
        else
            % Default color if not in map
            color = [0.5, 0.5, 0.5];
            logger(['Warning: No color defined for condition ' cond ', using default'], 'WARN');
        end
        
        scatter(zvals(idx,1), zvals(idx,2), 1, color, 'Marker', '.', 'DisplayName', cond);
    end
    
    hold off;
    legend('Location', 'best');
    title([comparison.title ' t-SNE Map (' experiment.name ')']);
    xlabel('t-SNE Dimension 1');
    ylabel('t-SNE Dimension 2');
    axis equal tight;
    
    % Ensure white background for export
    set(fig, 'Color', 'w');
    set(ax, 'Color', 'w');
    
    % Export the figure if needed
    if do_export
        filename = sprintf('scatter_plot_%s_%s.pdf', experiment.name, comparison.name);
        filepath = fullfile(export_folder, filename);
        logger(['Exporting scatter plot to: ' filepath], 'INFO');
        exportgraphics(fig, filepath, 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
end