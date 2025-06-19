%% Preamble
% This script will take all the analysed animals (in animal_list)
% and find the clusters id in the analysisstruct data.
% This script should run after the script_03__GOOD_analysis_concat_predictions.m.
% The output of this script is a .mat file containing a structure with cluster id for each animal, in a variable called cluster_vector.
% At the same time, the structure output should have the conditions.

%% Parameters
logger('Starting extract cluster vectors script', 'INFO');
%{
 animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
 
%}

clc, clear

GC = general_configs(); % load general configurations

%% Load Data
logger('Loading data', 'INFO');
% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');

% Load predictions
load(GC.filename_predictions, 'predictions', 'animal_condition_identifier');

% Load ratception structure
load(GC.filename_ratception, 'ratception_struct');



%% Prepare Data
logger('Preparing data', 'INFO');
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
clusters_struct_file = fullfile(save_folder, 'clusters_struct_high_density.mat');

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
logger('Looping through animals to extract clusters', 'INFO');
clusters_struct = struct();
for animal_idx = 1:length(animal_list)
    animal_ID = animal_list{animal_idx};
    clusters = get_clusters(animal_list, animal_ID, cond_inds, analysisstruct.highdensity_analysisstruct, upsamplig_factor, conditions);
    % store the clusters in a sturcture
    animal_ID_in_struct = ['ID_',animal_ID];
    clusters_struct.(animal_ID_in_struct) = clusters;

end
% Store the conditions in the same structure
clusters_struct.conditions = conditions;

% Save the structure
logger('Saving clusters structure', 'INFO');
save(clusters_struct_file, 'clusters_struct')

disp('DONE')
logger('Extract cluster vectors script completed', 'INFO');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions
