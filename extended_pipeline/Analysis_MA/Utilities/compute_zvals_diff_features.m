% Compute t-SNE using different feature sets and plot results
%% Initialization
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
export_folder = fullfile(GC.temp_root, 'figs_tsne_features');
if ~exist(export_folder, 'dir')
    mkdir(export_folder);
end

%% Define experiments and their conditions
experiments = struct();
experiments(1).name = 'BSFC_300hz';
experiments(1).conditions = {'B', 'S', 'F', 'C'};
experiments(1).folder = '0_preprocessing_BSFC_300hz';
experiments(1).plot_conditions = {'S', 'F', 'C'};  % SFC for plotting

experiments(2).name = 'BHNG_300hz';
experiments(2).conditions = {'B', 'H', 'N', 'G'};
experiments(2).folder = '0_preprocessing_BHNG_300hz';
experiments(2).plot_conditions = {'H', 'N', 'G'};  % HNG for plotting

%% Define feature sets
% Define extra features with their corresponding numbers
tsne_features = {'ja_dyadic_spectrograms','appearance_features_agg_score_whitened','pose_score',...
    'spectrogram_pcs_wl_head_angle','spectrogram_pcs_wl_trunk_angle',...
    'ja_eig_spectrogrampcs',...
    'ja_eig_wlpcs',...
    'ja_velocityhead_angle_difforder_100',...
    'ja_velocityhead_angle_difforder_300',...
    'ja_velocitytrunk_angle_difforder_100', ...
    'ja_velocitytrunk_angle_difforder_300',...
    'high_rear',...
    'ext_left_paw',...
    'RGroom',...
    'LGroom',...
    'face_groom_R',...
    'face_groom_L',...
    'guard_left_paw',...
    'lick_bite_left',...
    'weight_asymmetry',...
    'hunch_ratio',...
    'lateral_shift', ...
    'paw_clustering', ...
    };

num_feat = [10,6,10,15,15, 5, 5, 5, 5, 5];  % Updated to include single features

% Generate feature names for extra features
tsnefeatname = cell(0,1);
for ll = 1:numel(tsne_features)
    if numel(num_feat)>=ll
        for mm = 1:num_feat(ll)
            tsnefeatname{numel(tsnefeatname)+1} = strcat(tsne_features{ll},'_',num2str(mm));
        end
    else
        tsnefeatname{numel(tsnefeatname)+1} = tsne_features{ll};
    end
end

% Define important single features for feature set 2
important_single_features = {'guard_left_paw', 'lick_bite_left', 'weight_asymmetry', ...
                           'hunch_ratio', 'lateral_shift', 'paw_clustering'};

% Define feature sets
feature_sets = struct();
feature_sets(1).name = 'all_features';
feature_sets(1).description = 'All jt_features + extra_jt_features';

feature_sets(2).name = 'important_features';
feature_sets(2).description = 'jt_features + important single features (66 total)';

feature_sets(3).name = 'jt_features_only';
feature_sets(3).description = 'jt_features only';

% New feature sets that use ONLY extra features (no concatenation with jt_features)
feature_sets(4).name = 'hand_made';
feature_sets(4).description = 'Only hand-made single features from extra features';

feature_sets(5).name = 'only_important';
feature_sets(5).description = 'Only important single features from extra features';

% Define colors for conditions
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

    % Load analysis structure
    filename_analysis = fullfile(exp_data_folder, 'raw_concat_analysis.mat');
    if ~exist(filename_analysis, 'file')
        logger(['Analysis file not found: ' filename_analysis], 'ERROR');
        continue;
    end

    logger('Loading analysis structure', 'INFO');
    load(filename_analysis, 'analysisstruct');

    % Load predictions for condition identifiers
    filename_predictions = fullfile(exp_data_folder, 'agg_predictions.mat');
    if ~exist(filename_predictions, 'file')
        logger(['Predictions file not found: ' filename_predictions], 'ERROR');
        continue;
    end

    logger('Loading predictions', 'INFO');
    load(filename_predictions, 'predictions', 'animal_condition_identifier');

    % Process data for condition identification
    upsampling_factor = GC.repfactor;
    long_animal_frames_identifier = repelem(animal_condition_identifier, upsampling_factor);
    animal_list_used_after_analysis = long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});
    conditions = cellfun(@(x) x(end), animal_list_used_after_analysis, 'UniformOutput', false);

    % Create experiment-specific export folder
    exp_export_folder = fullfile(export_folder, experiment.name);
    if ~exist(exp_export_folder, 'dir')
        mkdir(exp_export_folder);
    end

    %% Process each feature set
    for feat_idx = 4%:length(feature_sets)
        feature_set = feature_sets(feat_idx);
        logger(['Computing t-SNE for feature set: ' feature_set.name ' (' experiment.name ')'], 'INFO');

        % Extract features based on feature set
        X_features = extract_features(analysisstruct, feature_set, tsnefeatname, important_single_features);

        if isempty(X_features)
            logger(['Warning: No features extracted for ' feature_set.name], 'WARN');
            continue;
        end

        logger(['Feature set size: ' num2str(size(X_features, 2)) ' features'], 'INFO');

        % Compute t-SNE
        perplexity = GC.perplexity;

        X_features = fillmissing(X_features, 'linear');

        if strcmp(feature_set.name, 'all_features')
            exaggeration = 5;
        elseif any(strcmp(feature_set.name, {'hand_made','only_important'}))
            exaggeration = 15;
        else
            exaggeration = 15;
        end

        logger(['Computing t-SNE with perplexity: ' num2str(perplexity), ' and Exaggeration: ', num2str(exaggeration)], 'INFO');

        zvals = tsne(X_features, "Perplexity", perplexity, 'Exaggeration', exaggeration, 'verbose', 1);

        % Save zvals to analysisstruct with appropriate name
        field_name = ['zValues_' feature_set.name];
        analysisstruct.(field_name) = zvals;

        % Plot results for this feature set
        plot_feature_set_results(zvals, conditions, experiment, feature_set, color_map, ...
                                exp_export_folder, visualize, do_export);
    end

    % Save updated analysisstruct
    logger(['Saving updated analysisstruct for ' experiment.name], 'INFO');
    save(filename_analysis, 'analysisstruct', '-append');

    logger(['Completed processing for experiment: ' experiment.name], 'INFO');
end

logger('Feature analysis complete', 'INFO');

%% Helper Functions

function X_features = extract_features(analysisstruct, feature_set, tsnefeatname, important_single_features)
    % Extract features based on the specified feature set

    % Start with jt_features (base features) - this is a matrix
    X_features = analysisstruct.jt_features;

    switch feature_set.name
        case 'all_features'
            % Add all extra features by concatenating the entire extra_jt_features matrix
            if isfield(analysisstruct, 'extra_jt_features') && ~isempty(analysisstruct.extra_jt_features)
                X_features = [X_features, analysisstruct.extra_jt_features];
            else
                logger('Warning: extra_jt_features not found or empty', 'WARN');
            end

        case 'important_features'
            % Add only important single features by finding their column indices
            if isfield(analysisstruct, 'extra_jt_features') && ~isempty(analysisstruct.extra_jt_features)
                % Find indices of important features in tsnefeatname
                important_indices = [];
                for i = 1:length(important_single_features)
                    feat_name = important_single_features{i};
                    % Find the index in tsnefeatname that matches this feature
                    feat_idx = find(strcmp(tsnefeatname, feat_name));
                    if ~isempty(feat_idx)
                        important_indices = [important_indices, feat_idx];
                    else
                        logger(['Warning: Important feature ' feat_name ' not found in feature list'], 'WARN');
                    end
                end

                % Extract the columns corresponding to important features
                if ~isempty(important_indices)
                    X_features = [X_features, analysisstruct.extra_jt_features(:, important_indices)];
                end
            else
                logger('Warning: extra_jt_features not found or empty', 'WARN');
            end

        case 'jt_features_only'
            % Use only jt_features (already assigned above)
            % X_features is already set to analysisstruct.jt_features

        case 'hand_made'
            % Use only the specified hand-made single features from extra_jt_features
            hand_made_features = {'high_rear', 'ext_left_paw', 'RGroom', 'LGroom', ...
                                  'face_groom_R', 'face_groom_L', 'guard_left_paw', 'lick_bite_left', ...
                                  'weight_asymmetry', 'hunch_ratio', 'lateral_shift', 'paw_clustering'};
            if isfield(analysisstruct, 'extra_jt_features') && ~isempty(analysisstruct.extra_jt_features)
                % Map names to indices in tsnefeatname
                idx = [];
                for i = 1:numel(hand_made_features)
                    feat_name = hand_made_features{i};
                    feat_idx = find(strcmp(tsnefeatname, feat_name));
                    if ~isempty(feat_idx)
                        idx = [idx, feat_idx];
                    else
                        logger(['Warning: Hand-made feature ' feat_name ' not found in feature list'], 'WARN');
                    end
                end
                if ~isempty(idx)
                    X_features = analysisstruct.extra_jt_features(:, idx); % ONLY extra features
                else
                    X_features = [];
                end
            else
                logger('Warning: extra_jt_features not found or empty', 'WARN');
                X_features = [];
            end

        case 'only_important'
            % Use only the important single features from extra_jt_features
            if isfield(analysisstruct, 'extra_jt_features') && ~isempty(analysisstruct.extra_jt_features)
                idx = [];
                for i = 1:length(important_single_features)
                    feat_name = important_single_features{i};
                    feat_idx = find(strcmp(tsnefeatname, feat_name));
                    if ~isempty(feat_idx)
                        idx = [idx, feat_idx];
                    else
                        logger(['Warning: Important feature ' feat_name ' not found in feature list'], 'WARN');
                    end
                end
                if ~isempty(idx)
                    X_features = analysisstruct.extra_jt_features(:, idx); % ONLY extra features
                else
                    X_features = [];
                end
            else
                logger('Warning: extra_jt_features not found or empty', 'WARN');
                X_features = [];
            end

        otherwise
            logger(['Unknown feature set: ' feature_set.name], 'ERROR');
            X_features = [];
    end
end

function plot_feature_set_results(zvals, conditions, experiment, feature_set, color_map, export_folder, visualize, do_export)
    % Plot t-SNE results for a specific feature set

    logger(['Plotting results for feature set: ' feature_set.name ' (' experiment.name ')'], 'INFO');

    % Create scatter plot for the specified conditions
    fig_name = ['t-SNE: ' feature_set.description ' (' experiment.name ')'];
    fig = figure('Name', fig_name, 'Color', 'w', 'Visible', visualize, 'Position', [100, 100, 800, 600]);

    % Set axes background to white
    ax = gca;
    set(ax, 'Color', 'w');

    hold on;

    % Plot each condition in plot_conditions
    legend_entries = {};
    for i = 1:length(experiment.plot_conditions)
        cond = experiment.plot_conditions{i};

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
            color = [0.5, 0.5, 0.5];
            logger(['Warning: No color defined for condition ' cond ', using default'], 'WARN');
        end

        scatter(zvals(idx,1), zvals(idx,2), 10, color, 'Marker', '.', 'DisplayName', cond);
        legend_entries{end+1} = cond;
    end

    hold off;
    legend(legend_entries, 'Location', 'best');
    title([feature_set.description ' - ' strjoin(experiment.plot_conditions, '') ' (' experiment.name ')']);
    xlabel('t-SNE Dimension 1');
    ylabel('t-SNE Dimension 2');
    axis equal tight;

    % Ensure white background for export
    set(fig, 'Color', 'w');
    set(ax, 'Color', 'w');

    % Export the figure if needed
    if do_export
        filename = sprintf('tsne_%s_%s_%s.pdf', experiment.name, feature_set.name, strjoin(experiment.plot_conditions, ''));
        filepath = fullfile(export_folder, filename);
        logger(['Exporting plot to: ' filepath], 'INFO');
        exportgraphics(fig, filepath, 'ContentType', 'vector', 'BackgroundColor', 'white');
    end
end
