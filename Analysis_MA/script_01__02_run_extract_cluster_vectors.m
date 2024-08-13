%% Preamble
% This script will take all the analysed animals (in animal_list)
% and find the clusters id in the analysisstruct data.
% This script should run after the script_03__GOOD_analysis_concat_predictions.m.
% The output of this script is a .mat file containing a structure with cluster id for each animal, in a variable called cluster_vector.
% At the same time, the structure output should have the conditions.

%% Parameters

%{
 animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
 
%}

clc, clear

GC = general_configs(); % load general configurations

%% Load Data
% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');

% Load predictions
load(GC.filename_predictions, 'predictions', 'animal_condition_identifier');

% Load ratception structure
load(GC.filename_ratception, 'ratception_struct');



%% Prepare Data

% get the animal list and the condition
animals_and_conditions = unique(animal_condition_identifier, 'stable');

% Extract unique animal IDs and conditions
[animal_ids, ~, animal_indices] = unique(cellfun(@(x) x(1:end-2), animals_and_conditions, 'UniformOutput', false));
conditions = cellfun(@(x) x(end), animals_and_conditions, 'UniformOutput', false);

% Initialize cluster vector for each animal
save_folder = GC.preprocessing_rootpath;
if ~exist(save_folder, 'dir')
    mkdir(save_folder)
end

% check if file exists and if you want to run it again
clusters_struct_file = fullfile(save_folder, 'clusters_struct.mat');

%{
 if exist(clusters_struct_file, 'file')
    answer= input("The file already exists, you want to load it? Y/N [Y]:",'s');
    if strcmp(answer, 'Yes')   
        cluster_data = load(clusters_struct_file);
        disp('DONE')       
    else
        disp('DONE')
        return
    end
end 
%}


% TODO: check if we add more animals, to just append to the struct the new animals and condition


% set variables:
% Given frame rate
ori_frame_rate = GC.frame_rate; 
upsamplig_factor = GC.repfactor;
%upsamplig_factor = ceil(analysisstruct.mocapstruct_reduced_agg{1, 1}.fps/ori_frame_rate); 

long_animal_frames_identifier = repelem(animal_condition_identifier,upsamplig_factor);


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

% Loop through the animals
clusters_struct = struct();
for animal_idx = 1:length(animal_list)
    animal_ID = animal_list{animal_idx};
    clusters = get_clusters(animal_list, animal_ID, cond_inds, analysisstruct, upsamplig_factor, conditions);
    % store the clusters in a sturcture
    animal_ID_in_struct = ['JH_', animal_ID];
    clusters_struct.(animal_ID_in_struct) = clusters;

end
% Store the conditions in the same structure
clusters_struct.conditions = conditions;

% Save the structure
save(clusters_struct_file, 'clusters_struct')

disp('DONE')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions
function clusters = get_clusters(animal_list, animal_ID, cond_inds, analysisstruct, upsamplig_factor, conditions)
    % 1. Find index of animal '328'
    animal_idx = find(strcmp(animal_list, animal_ID));

    % 2. Get the condition indices for animal '328'
    specific_cond_inds = find(ismember(cond_inds,animal_idx));

    % 3. Extract specific frames and outcomes based on the condition indices
    % analysisstruct.frames_with_good_tracking is at least diff 50 frames
    frames = analysisstruct.frames_with_good_tracking{1,1}(specific_cond_inds);
    outcomes = analysisstruct.annot_reordered_matched{1, 1}(specific_cond_inds); % sampled every 50
    
    % get total number of frames: for this we need to search for the prediction.mat file in the animal folder
    % load the prediction file
    % load(fullfile("D:\test_CAPTURE",animal_ID, 'predictions.mat' ), "predictions");
    this_exp_cond = conditions{animal_idx};
    if strcmp(this_exp_cond, 'F')
        folder_exp_cond = 'PFA';
    elseif strcmp(this_exp_cond, 'S')
        folder_exp_cond = 'saline';
    end

    disp(['Running Animal ', animal_ID, ' - condition ', folder_exp_cond])
    % server_folder = fullfile("H:\DANNCE\6cam_behavior",folder_exp_cond, animal_ID,"DANNCE_ready\DANNCE\predict_results_net_8" );
    % use com to read the number of frames
    % TODO: check if this is correct, becuase for some reasong for an example animal com has 216001 frames,
    % but the predictions has 216000
    % com_filename = fullfile(server_folder, 'com3d_used.mat');
    % COM = load(com_filename, "com");
    n_frames_total = 30*60*100; % TODO, fix later


    % Initialize the clusters array
    clusters = zeros(1, n_frames_total * upsamplig_factor);

    for i = 1:length(frames)
        % frame_val = ceil(frames_this_idx(i) / upsamplig_factor);

        frame_val = frames(i);

        if i == 1
            relative_start = frame_val; % This will be adjusted below
            if animal_idx > 1 && frame_val > analysisstruct.tsnegranularity  % 50
                previous_frame = analysisstruct.frames_with_good_tracking{1,1}(specific_cond_inds(1) - 1); % diff of 50 frames 
                % relative_start = frame_val - ceil(previous_frame / 3);
                relative_start = frame_val - ceil(previous_frame);
                if relative_start == analysisstruct.tsnegranularity % means it's only one frame
                    relative_start = 1;
                elseif relative_start > analysisstruct.tsnegranularity % meaning at the beginning there are no clusters
                    % relative_start = previous_frame + 50; TODO, still I don't
                    % know how to go on with this
                    keyboard %#ok<KEYBOARDFUN>
                end
            end
            clusters(1:relative_start) = outcomes(i);
            % clusters(relative_start:frame_val) = outcomes(i);

            current_index = relative_start + 1;
        else

            % set the gap: cehck if the previous frame has a gap of 50
            gap = ceil((frames(i) - frames(i-1)));
            if gap > ceil(analysisstruct.tsnegranularity) 
                % keyboard % TODO
                clusters(current_index:current_index+gap-1) = 0; % fill gaps with zeros
                current_index = current_index + gap;
            else
                clusters(current_index:current_index+gap-1) = outcomes(i);
                % current_index = current_index + gap-1;
                current_index = current_index + gap;

            end
        end
    end
    % Downsample the clusters
    clusters = downsample(clusters, upsamplig_factor);
end 