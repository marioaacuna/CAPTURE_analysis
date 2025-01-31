%% Parameters
logger('Starting extra features', 'INFO');
%{
 animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
 
%}

clc, clear

GC = general_configs(); % load general configurations

%% Load Data
logger('Loading data', 'INFO');
% % Load analysis structure
% load(GC.filename_analysis, 'analysisstruct');
% 
% % Load predictions
% load(GC.filename_predictions, 'predictions', 'animal_condition_identifier');

% Load ratception structure
load(GC.filename_ratception, 'ratception_struct');
mocapstruct = ratception_struct;
%% Run extra features
clear ratception_struct
mocapstruct.modular_cluster_properties.clipped_index{8} = 1:size(mocapstruct.aligned_mean_position,1 );

% Set inputs
temp_dir = 'D:\CAPTURE\_temp\250130_extra_features';
savefilename =fullfile(temp_dir,'myMLfeatures.mat');
directory_here = temp_dir;
overwrite_coefficient=0;

MLmatobj_extra = create_extra_behavioral_features(mocapstruct,'concate_mice',savefilename,overwrite_coefficient,directory_here);
jt_features_extra = load_extra_tsne_features(mocapstruct,MLmatobj_extra,analysisparams);

