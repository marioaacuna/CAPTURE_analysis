%% Preamble
% This script processes each animal-condition combination separately.
% Unlike the concatenated version, this script creates individual 
% analysisstruct for each animal in each condition, which can later
% be concatenated for cross-animal analysis. Each animal-condition
% gets its own dedicated folder structure for organized analysis.

%% preINIT
close all
clc
clear
GC = general_configs();

%% INIT
% Get prediction concatenation settings
settings = get_settings_concat_preds();

% Unpack settings into variables
overwrite_ratception = settings.overwrite_ratception;
overwrite_MLmatobjfile = settings.overwrite_MLmatobjfile;
overwrite_coefficient = settings.overwrite_coefficient;
do_extra_features = settings.do_extra_features;

init_frame_rate = GC.frame_rate; % effective frame rate of videos
rootpath = GC.preprocessing_rootpath;
if ~exist(rootpath, 'dir'), mkdir(rootpath); end

% Read metadata from YAML
metadata = readExperimentMetadata();

animal_name = GC.ratception_name;

input_params = struct();
input_params.SpineF_marker = 'SpineF';
input_params.SpineM_marker = 'SpineM';
input_params.repfactor = GC.repfactor; % round(300/init_frame_rate);
input_params.conversion_factor = 1;

linkname = GC.linkname;
ratname = 'myrat'; % 'test_mouse';

analysisparams.tsnegranularity = GC.tsnegranularity; % 25:default
% Mario: it seems that high number of frames (between 100 - 120) yield better results

% Get all conditions from metadata
conditions = fieldnames(metadata.conditions);

% Initialize storage for final concatenation
all_analysisstruct = {};
all_jt_features_extra = {};
all_animal_identifiers = {};

%% Main loop: Process each condition and animal separately
for condition_idx = 1:length(conditions)
    condition_name = conditions{condition_idx};
    condition_data = metadata.conditions.(condition_name);
    condition_letter = condition_to_letter(condition_name);
    
    fprintf('Processing condition: %s\n', condition_name);
    
    % Loop through animals in this condition
    try
        for animal_idx = 1:length(condition_data.animals)
            animal = condition_data.animals(animal_idx);
            animal_id = animal.id;
            
            % Check if animal has 6cam recording
            if ~animal.has_6cam
                warning('No 6cam recording found for animal %s in condition %s', animal_id, condition_name);
                continue;
            end
            
            % Create animal-condition identifier
            animal_identifier = sprintf('%s_%s', animal_id, condition_letter);
            fprintf('  Processing animal: %s\n', animal_identifier);
            
            % Create specialized folder for this animal-condition
            animal_condition_folder = fullfile(rootpath, animal_identifier);
            if ~exist(animal_condition_folder, 'dir')
                mkdir(animal_condition_folder);
            end
            
            % First construct path to animal folder
            animal_path = fullfile(metadata.data_root_dir, ...
                condition_data.data_path, ...
                animal.path);
            
            % Get the date folder (assuming it's the only folder in there that's a date)
            date_folders = dir(animal_path);
            date_folders = date_folders([date_folders.isdir]); % Only get directories
            date_folders = date_folders(~ismember({date_folders.name}, {'.', '..'})); % Remove . and ..
            
            if isempty(date_folders)
                warning('No date folder found for animal %s in condition %s', animal_id, condition_name);
                continue;
            end
            
            date_folder = date_folders(1).name;
            % Construct full path to predictions
            filename_predictions = fullfile(animal_path, ...
                date_folder, ...
                'DANNCE', 'predict_results', ...
                'predictions.mat');
            
            % Check if file exists
            if ~exist(filename_predictions, 'file')
                warning('File not found: %s', filename_predictions);
                continue;
            end
            
            % Set up paths for this animal-condition
            filename_ratception = fullfile(animal_condition_folder, 'ratception_predictions.mat');
            roothpath_CAPTURE = animal_condition_folder;
            
            %% RATCEPTION preprocessing
            % Run prepro if it doesn't exist
            if ~exist(filename_ratception, "file") || overwrite_ratception
                fprintf('    Running Pre-Pro for %s\n', animal_identifier);
                ratception_struct = preprocess_dannce(filename_predictions, filename_ratception, animal_name, input_params);
            else
                fprintf('    Loading previously analyzed prepro data for %s\n', animal_identifier);
                load(filename_ratception);
            end
            
            %% Load Mocapstruct and set up paths
            mocapstruct = ratception_struct;
            clear ratception_struct % clear memory
            
            coefficient_file = fullfile(roothpath_CAPTURE, 'coefficients.mat');
            
            if do_extra_features
                % In case you want to do some extra features
                savefilename_extra = fullfile(roothpath_CAPTURE, 'myextratsnefeature', 'extraMLFeatures.mat');
                eigenposture_save_filder = fullfile(roothpath_CAPTURE, 'myextratsnefeature');
                mkdir(eigenposture_save_filder)
            end
            
            savefilename_features = fullfile(roothpath_CAPTURE);
            MLmatobjfile = fullfile(savefilename_features, 'myMLfeatures.mat');
            
            %% Create behavioral features
            % This determines the set of frames to use -- in general if the animal is
            % resting for too long it will cause errors
            mocapstruct.modular_cluster_properties.clipped_index{8} = 1:size(mocapstruct.aligned_mean_position, 1);
            % to control the wavelet parameters, you can change the properties in the
            % compute_wl_transform_features file
            
            if ~exist(MLmatobjfile, 'file') || overwrite_MLmatobjfile
                fprintf('    Creating behavioral features for %s\n', animal_identifier);
                MLmatobj = create_behavioral_features(mocapstruct, coefficient_file, overwrite_coefficient, linkname);
                save(MLmatobjfile, 'MLmatobj', '-v7.3')
            else
                fprintf('    Loading ML features for %s\n', animal_identifier);
                MLmatobj = load(MLmatobjfile, 'MLmatobj');
                MLmatobj = MLmatobj.MLmatobj;
            end
            
            %% Compute analysis structure
            fprintf('    Computing tsne features for %s\n', animal_identifier);
            % subselect a particular set of features
            analysisstruct = compute_tsne_features(MLmatobj, mocapstruct, analysisparams);
            
            %% Create extra features
            if do_extra_features
                fprintf('    Creating extra behavioral features for %s\n', animal_identifier);
                ML_extra_obj = create_extra_behavioral_features(mocapstruct, 'concate_mice', savefilename_extra, overwrite_coefficient, eigenposture_save_filder);
                jt_features_extra = load_extra_tsne_features(mocapstruct,ML_extra_obj,analysisparams);
                analysisstruct.extra_jt_features = jt_features_extra;
            else
                jt_features_extra = [];
            end
            % save analysisstruct in animal-condition folder
            save(fullfile(animal_condition_folder, 'analysisstruct.mat'), "analysisstruct");
        
            %% Store results for this animal-condition
            all_analysisstruct{end+1} = analysisstruct;
            all_jt_features_extra{end+1} = jt_features_extra;
            all_animal_identifiers{end+1} = animal_identifier;
            
            fprintf('    Completed analysis for %s\n', animal_identifier);
            
        end
        
    catch me
        fprintf('Error processing condition %s, animal %s: %s\n', condition_name, animal_id, me.message);
        continue
    end
end

%% Summary
fprintf('\n=== Analysis Summary ===\n');
fprintf('Total animal-condition combinations processed: %d\n', length(all_analysisstruct));
fprintf('Animal identifiers:\n');
for i = 1:length(all_animal_identifiers)
    fprintf('  %d: %s\n', i, all_animal_identifiers{i});
end

%% Save results
results_file = fullfile(rootpath, 'per_animal_analysis_results.mat');
save(results_file, 'all_analysisstruct', 'all_jt_features_extra', 'all_animal_identifiers', '-v7.3');
fprintf('Results saved to: %s\n', results_file);

fprintf('\n=== Analysis Complete ===\n');
fprintf('Each animal-condition has been processed individually.\n');
fprintf('Results stored in individual folders under: %s\n', rootpath);
fprintf('Next step: Use these results for concatenated analysis or individual animal studies.\n');

%% Functions
% First, let's create the condition_to_letter function to handle all conditions
function letter = condition_to_letter(condition_name)
    switch lower(condition_name)
        case 'saline'
            letter = 'S';
        case 'formalin'
            letter = 'F';
        case 'baseline'
            letter = 'B';
        case 'sham'
            letter = 'H';  % H for sHam to avoid confusion with S for Saline
        case 'sni'
            letter = 'N';  % N for sNi
        case 'car'
            letter = 'C';
        case 'gbp'
            letter = 'G';
        otherwise
            error('Unknown condition: %s', condition_name);
    end
end
