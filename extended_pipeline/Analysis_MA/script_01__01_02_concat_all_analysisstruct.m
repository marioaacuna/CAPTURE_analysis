%% Preamble
% This script concatenates all analysisstruct data from individual 
% animal-condition analyses into a single dataset for cross-animal analysis.
% It combines jt_features and extra_jt_features for each animal-condition
% and maintains tracking of which frames belong to which animal-condition.

%% preINIT
close all
clc
clear
GC = general_configs();

%% INIT
rootpath = GC.preprocessing_rootpath;

% Check if the per-animal analysis results exist
results_file = fullfile(rootpath, 'per_animal_analysis_results.mat');
if ~exist(results_file, 'file')
    error('Per-animal analysis results not found. Please run script_01__01_run_analysisstruct_per_animal.m first.');
end

% Load the per-animal analysis results
fprintf('Loading per-animal analysis results...\n');
load(results_file, 'all_analysisstruct', 'all_jt_features_extra', 'all_animal_identifiers');

fprintf('Found %d animal-condition combinations to concatenate.\n', length(all_analysisstruct));

%% Process and concatenate all data
% Initialize the concatenated data matrix
D = [];
% Initialize frame tracking for animal-condition identification
frame_to_animal_condition = {};
animal_condition_frame_counts = [];

fprintf('\nConcatenating data from all animal-conditions:\n');

for i = 1:length(all_analysisstruct)
    animal_identifier = all_animal_identifiers{i};
    analysisstruct = all_analysisstruct{i};
    jt_features_extra = all_jt_features_extra{i};

    fprintf('  Processing %s...\n', animal_identifier);

    % Check if we have extra features for this animal-condition
    if isempty(jt_features_extra)
        % If no extra features, just use the regular jt_features
        d = analysisstruct.jt_features;
        fprintf('    Using only jt_features (no extra features available)\n');
    else
        % Concatenate jt_features and extra_jt_features horizontally
        d = [analysisstruct.jt_features, jt_features_extra];
        fprintf('    Concatenated jt_features (%d) + extra_jt_features (%d) = %d features\n', ...
            size(analysisstruct.jt_features, 2), size(jt_features_extra, 2), size(d, 2));
    end

    % Get number of frames for this animal-condition
    num_frames = size(d, 1);
    animal_condition_frame_counts(i) = num_frames;

    % Create frame identifiers for this animal-condition
    animal_frames = repmat({animal_identifier}, num_frames, 1);

    % Concatenate data vertically
    if isempty(D)
        D = d;
        frame_to_animal_condition = animal_frames;
    else
        D = [D; d];
        frame_to_animal_condition = [frame_to_animal_condition; animal_frames];
    end

    fprintf('    Added %d frames (total frames: %d)\n', num_frames, size(D, 1));
end

%% Create comprehensive concatenated analysisstruct
fprintf('\n=== Creating concatenated analysisstruct ===\n');

% Initialize the concatenated analysisstruct using the first one as template
concatenated_analysisstruct = struct();

% Get the template field names from the first analysisstruct
template_struct = all_analysisstruct{1};
field_names = fieldnames(template_struct);

% Initialize fields that will be concatenated
concatenated_analysisstruct.frames_with_good_tracking = {};
concatenated_analysisstruct.frames_tracking_appendages = [];
concatenated_analysisstruct.subset_of_points_to_plot_tsne_capped = {};
concatenated_analysisstruct.subset_of_points_to_plot_tsne_move = {};
concatenated_analysisstruct.condition_inds = {};
concatenated_analysisstruct.jt_features = [];
concatenated_analysisstruct.jt_features_raw = [];
concatenated_analysisstruct.jt_features_mean = [];
concatenated_analysisstruct.jt_features_std = [];
concatenated_analysisstruct.file_sizes = {};
concatenated_analysisstruct.mocapstruct_reduced_agg = {};
concatenated_analysisstruct.extra_jt_features = [];

% Take tsnefeat_name from the first analysisstruct (should be the same for all)
concatenated_analysisstruct.tsnefeat_name = template_struct.tsnefeat_name;

% Initialize tracking variables
current_frame_offset = 0;

fprintf('Processing %d animal-condition combinations:\n', length(all_analysisstruct));

for i = 1:length(all_analysisstruct)
    animal_identifier = all_animal_identifiers{i};
    analysisstruct = all_analysisstruct{i};
    jt_features_extra = all_jt_features_extra{i};
    
    fprintf('  Processing %s (animal-condition %d)...\n', animal_identifier, i);
    
    % Get number of frames for this animal-condition
    num_frames = size(analysisstruct.jt_features, 1);
    
    % 1. frames_with_good_tracking - add as new cell in second dimension
    concatenated_analysisstruct.frames_with_good_tracking{i} = analysisstruct.frames_with_good_tracking{1};
    
    % 2. frames_tracking_appendages - concatenate sequentially with offset
    new_frame_indices = (current_frame_offset + 1):(current_frame_offset + num_frames);
    concatenated_analysisstruct.frames_tracking_appendages = [concatenated_analysisstruct.frames_tracking_appendages; new_frame_indices'];
    
    % 3. subset_of_points_to_plot_tsne_capped - add as new cell
    if isfield(analysisstruct, 'subset_of_points_to_plot_tsne_capped') && ~isempty(analysisstruct.subset_of_points_to_plot_tsne_capped)
        concatenated_analysisstruct.subset_of_points_to_plot_tsne_capped{i} = analysisstruct.subset_of_points_to_plot_tsne_capped{1};
    else
        error('subset_of_points_to_plot_tsne_capped')
        break
        
        % concatenated_analysisstruct.subset_of_points_to_plot_tsne_capped{i} = new_frame_indices;
    end
    
    % 4. subset_of_points_to_plot_tsne_move - add as new cell
    if isfield(analysisstruct, 'subset_of_points_to_plot_tsne_move') && ~isempty(analysisstruct.subset_of_points_to_plot_tsne_move)
        concatenated_analysisstruct.subset_of_points_to_plot_tsne_move{i} = analysisstruct.subset_of_points_to_plot_tsne_move{1} + current_frame_offset;
    else
        concatenated_analysisstruct.subset_of_points_to_plot_tsne_move{i} = new_frame_indices;
    end
    
    % 5. condition_inds - create cell with animal-condition ID repeated
    condition_id_vector = repmat(i, 1, num_frames);
    concatenated_analysisstruct.condition_inds{i} = condition_id_vector;
    
    % 6. jt_features - concatenate vertically
    concatenated_analysisstruct.jt_features = [concatenated_analysisstruct.jt_features; analysisstruct.jt_features];
    
    % 7. jt_features_raw - concatenate vertically
    if isfield(analysisstruct, 'jt_features_raw') && ~isempty(analysisstruct.jt_features_raw)
        concatenated_analysisstruct.jt_features_raw = [concatenated_analysisstruct.jt_features_raw; analysisstruct.jt_features_raw];
    else
        concatenated_analysisstruct.jt_features_raw = [concatenated_analysisstruct.jt_features_raw; analysisstruct.jt_features];
    end
    
    % 8. jt_features_mean - concatenate in second dimension (horizontally)
    if isfield(analysisstruct, 'jt_features_mean') && ~isempty(analysisstruct.jt_features_mean)
        if isempty(concatenated_analysisstruct.jt_features_mean)
            concatenated_analysisstruct.jt_features_mean = analysisstruct.jt_features_mean;
        else
            concatenated_analysisstruct.jt_features_mean = [concatenated_analysisstruct.jt_features_mean; analysisstruct.jt_features_mean];
        end
    end
    
    % 9. jt_features_std - concatenate in second dimension (horizontally)
    if isfield(analysisstruct, 'jt_features_std') && ~isempty(analysisstruct.jt_features_std)
        if isempty(concatenated_analysisstruct.jt_features_std)
            concatenated_analysisstruct.jt_features_std = analysisstruct.jt_features_std;
        else
            concatenated_analysisstruct.jt_features_std = [concatenated_analysisstruct.jt_features_std; analysisstruct.jt_features_std];
        end
    end
    
    % 10. file_sizes - add as new cell with value 1
    concatenated_analysisstruct.file_sizes{i} = 1;
    
    % 11. mocapstruct_reduced_agg - add as new cell
    if isfield(analysisstruct, 'mocapstruct_reduced_agg') && ~isempty(analysisstruct.mocapstruct_reduced_agg)
        concatenated_analysisstruct.mocapstruct_reduced_agg{i} = analysisstruct.mocapstruct_reduced_agg{1};
    end
    
    % 12. extra_jt_features - concatenate vertically
    if ~isempty(jt_features_extra)
        concatenated_analysisstruct.extra_jt_features = [concatenated_analysisstruct.extra_jt_features; jt_features_extra];
    end
    
    % Update frame offset for next iteration
    current_frame_offset = current_frame_offset + num_frames;
    
    fprintf('    Added %d frames (total: %d)\n', num_frames, current_frame_offset);
end

% Add metadata about the concatenation
concatenated_analysisstruct.concatenation_info.total_frames = size(concatenated_analysisstruct.jt_features, 1);
concatenated_analysisstruct.concatenation_info.total_features = size(concatenated_analysisstruct.jt_features, 2);
concatenated_analysisstruct.concatenation_info.num_animal_conditions = length(all_animal_identifiers);
concatenated_analysisstruct.concatenation_info.animal_identifiers = all_animal_identifiers;
concatenated_analysisstruct.concatenation_info.creation_date = datestr(now);

% Copy other fields that might exist but weren't specifically mentioned
other_fields = setdiff(field_names, {'tsnefeat_name', 'frames_with_good_tracking', 'frames_tracking_appendages', ...
    'subset_of_points_to_plot_tsne_capped', 'subset_of_points_to_plot_tsne_move', 'condition_inds', ...
    'jt_features', 'jt_features_raw', 'jt_features_mean', 'jt_features_std', 'file_sizes', ...
    'mocapstruct_reduced_agg', 'extra_jt_features'});

for field_idx = 1:length(other_fields)
    field_name = other_fields{field_idx};
    if isfield(template_struct, field_name)
        concatenated_analysisstruct.(field_name) = template_struct.(field_name);
        fprintf('  Copied field: %s\n', field_name);
    end
end



%% Save concatenated results
fprintf('\n=== Saving Results ===\n');
concatenated_file = fullfile(rootpath, 'concatenated_analysis_results.mat');
save(concatenated_file, 'D','frame_to_animal_condition', 'animal_condition_frame_counts', ...
     'all_animal_identifiers', 'concatenated_analysisstruct', '-v7.3');

fprintf('Concatenated results saved to: %s\n', concatenated_file);

%% Create condition indices for analysis (based on condition letters)
% Extract condition information from animal identifiers
unique_conditions = {};
condition_indices = zeros(size(frame_to_animal_condition));

fprintf('\nCreating condition indices based on condition letters:\n');
for i = 1:length(frame_to_animal_condition)
    animal_condition = frame_to_animal_condition{i};
    
    % Extract condition letter (last character after underscore)
    parts = split(animal_condition, '_');
    condition_letter = parts{end};
    
    % Find or create condition index
    condition_idx = find(strcmp(unique_conditions, condition_letter));
    if isempty(condition_idx)
        unique_conditions{end+1} = condition_letter;
        condition_idx = length(unique_conditions);
    end
    
    condition_indices(i) = condition_idx;
end

% Add condition information to the concatenated structure
concatenated_analysisstruct.condition_indices_by_letter = condition_indices;
concatenated_analysisstruct.unique_conditions = unique_conditions;

fprintf('Found %d unique conditions: %s\n', length(unique_conditions), strjoin(unique_conditions, ', '));

% Display condition distribution
for i = 1:length(unique_conditions)
    condition_count = sum(condition_indices == i);
    fprintf('  Condition %s: %d frames\n', unique_conditions{i}, condition_count);
end

%% Final save with all information
save(concatenated_file, 'D', 'frame_to_animal_condition', 'animal_condition_frame_counts', ...
     'all_animal_identifiers', 'concatenated_analysisstruct', 'condition_indices', ...
     'unique_conditions', '-v7.3');

fprintf('\n=== Concatenation Complete ===\n');
fprintf('Data matrix D: %d frames × %d features\n', size(D, 1), size(D, 2));
fprintf('Frame tracking: %d entries\n', length(frame_to_animal_condition));
fprintf('Concatenated analysisstruct created with proper field structure\n');
fprintf('Ready for cross-animal analysis (t-SNE, clustering, etc.)\n');


%% get Zvals and do the rest of CAPTURE
