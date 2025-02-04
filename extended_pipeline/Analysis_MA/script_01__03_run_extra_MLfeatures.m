%% Parameters
logger('Starting extra features', 'INFO');
%{
 animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
 
%}

clc, clear

GC = general_configs(); % load general configurations

%% Load Data
logger('Loading data', 'INFO');
% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');
% 
% % Load predictions
% load(GC.filename_predictions, 'predictions', 'animal_condition_identifier');

answer = questdlg('Create new ML extra features or Overwrite existing? - you need a lot of RAM Y/N [Y]', ...
    'Overwrite Confirmation', ...
    'Yes', 'No', 'No');  % 'No' is default button
if strcmp(answer, 'Yes')
    overwriteML = true;
else
    overwriteML = false;
end

% Load ratception structure
load(GC.filename_ratception, 'ratception_struct');
mocapstruct = ratception_struct;
clear ratception_struct
%% Run extra features

mocapstruct.modular_cluster_properties.clipped_index{8} = 1:size(mocapstruct.aligned_mean_position,1 );

% Set inputs
temp_dir = 'D:\CAPTURE\_temp\250130_extra_features';
savefilename =fullfile(temp_dir,'myMLfeatures.mat');
directory_here = temp_dir;
overwrite_coefficient=0;

if overwriteML
    MLmatobj_extra = create_extra_behavioral_features(mocapstruct,'concate_mice',savefilename,overwrite_coefficient,directory_here);
else
    MLmatobj_extra =matfile(savefilename);
end

%% load extra features
% analysisparams.tsnegranularity = 50; 
% MLmatobj_extra = load(savefilename);
jt_features_extra = load_extra_tsne_features(mocapstruct,MLmatobj_extra,analysisstruct);

% Save jt_features_extra 
extra_filename = fullfile(temp_dir, 'jt_features_extra.mat');
save(extra_filename, "jt_features_extra")
%% plot dim red
% look at tsne of these added features
zvals_extra = tsne(jt_features_extra);
% or the combination
%zvals_extra_combined = tsne(cat(2,analysisstruct.jt_features,jt_features_extra));
figure(2)
plot(zvals_extra(:,1),zvals_extra(:,2),'ob','MarkerFaceColor','b', 'MarkerSize',1)
analysisstruct.zValues_extra = zvals_extra;

zvals_extra_combined = tsne(cat(2,analysisstruct.jt_features,jt_features_extra));

figure(3)
plot(zvals_extra_combined(:,1),zvals_extra_combined(:,2),'ob','MarkerFaceColor','b', 'MarkerSize',1)
