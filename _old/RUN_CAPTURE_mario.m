%% Overall preamble
% the main folder is here:
% T:\Mario\DANNCE\predictions\CAPTURE
% WE have to run the preprocessing of the DANNCE data first. 
% then we can run CAPTURE.
% 
clc, clear, close all

do_animation_poses = 0;
plot_poses = 1;
overwrite_ratception = 0;
overwrite_coefficient= 0;
overwrite_MLmatobjfile = 0; % for some reason the tsne doe not work if we load it


% the default is 'predictions'
folfer_to_look_predictions = 'predict_results';

%% 1. You have to create a ratception struct. 
animal_ID = '322';
% date = '230508';
% animal_name = 'mario_mouse'; % this is the animal to be used for preprocessing and on. NOT the animal ID of interest
% animal_name = 'mario_mouse_14_pts'; % this is the animal to be used for preprocessing and on with only 14 pints.. tail out. NOT the animal ID of interest
animal_name = 'mario_mouse22';
% rootpath = fullfile('T:\Marta\test_Formalin\dannce\040822', animal_ID,'DANNCE\predict_results\3_predicted_with_normal_and_pain_frames');
rootpath = uigetdir('H:\DANNCE\6cam_behavior\', 'Folder prediction');
% rootpath = fullfile('H:\DANNCE\230508', animal_ID,'DANNCE_ready','DANNCE',folfer_to_look_predictions);
%rootpath = fullfile('T:\Marta\test_Formalin\dannce\040822', animal_ID,'DANNCE',folfer_to_look_predictions);
filename_predictions = fullfile(rootpath, "predictions.mat");
filename_ratception = fullfile(rootpath, "ratception_prediction.mat");
% inputs for features

%% Preprocess
% input parameters
% input_params: a struct containing experiment specific
% information: the markers to use as
% SpineF (SpineF_marker) and SpineM (SpineM_marker)
% to align the animal, the relative framerate to 300
% Hz (repfactor = 300/experiment framerate),
% conversion factor (to scale the outputs)
input_params = struct();
input_params.SpineF_marker = 'SpineF';
input_params.SpineM_marker = 'SpineM';
% input_params.repfactor = 300/30;
input_params.repfactor = round(300/120);

input_params.conversion_factor = 1;
% Run prepro if it doesn't exist
if ~exist(filename_ratception, "file") || overwrite_ratception
    disp('%% Running Pre-Pro %% ')
    ratception_struct = preprocess_dannce(filename_predictions,filename_ratception,animal_name,input_params);
else
    disp('Loading previously analysied prepro data')
   load(filename_ratception);
end
%% Preamble CAPTURE
% function [analysisstruct,hierarchystruct] =  CAPTURE_quickdemo(inputfile,ratnames,coefficientfilename,linkname)
% File to generate tsne features and run reembedding on a mouse
%      inputfile: a .mat file that contains a preprocessed dannce struct
%                 (see preprocess_dannce)
%      ratnames: a string containing the name of the experiment/rat to be
%                used in saving files
%      coefficientnames: preexisting names of coefficients or file to save
%                      tsne coefficients to. 
%      linkname: the name of the animal (ie kyle_mouse or 'rats' or 'bird'
%                 to be used in computing tsne features)
%      
% This repository contains contributions from the following FEX/Open source contributions which are included:
%Chronux
%Pca Randomized
%MTimesX
%othercolor
%Motionmapper
%Structvars
% ---------------------------
% (C) Jesse D Marshall 2020
%     Harvard University 
% 
% load mocap file
% if isempty(inputfile)
% datafile = ...
%     load('C:\Users\Jesse Marshall\Documents\GitHub\Movement_analysis\Cortex_analysis\DemoRepo\Data\nolj_Recording_day8_caff1_nolj_imputed.mat');
% other_mat = load('D:\nolj_Recording_day8_caff1_nolj_imputed.mat');
% else
%     datafile = load(inputfile);
%     if isstruct(datafile) && numel(fieldnames(datafile)) == 1
%         fname = fieldnames(datafile);
%         datafile = datafile.(fname{1});
%     end
% end
% if isempty(ratnames)
%     ratname = 'myrat';
% else
%     ratname = ratnames;
% end

%% Run Capture from here on
% inputs
roothpath_CAPTURE = fullfile(rootpath, 'CAPTURE');
if ~exist(roothpath_CAPTURE, 'dir'), mkdir(roothpath_CAPTURE); end
coefficient_file = fullfile(roothpath_CAPTURE,'coefficients.mat');
linkname = 'mario_mouse22';
ratname ='myrat';% 'test_mouse';

%% Load Mocapstruct
mocapstruct = ratception_struct;

% In case you want to do some extra features
savefilename_extra = fullfile(rootpath, 'CAPTURE', 'myextratsnefeature', 'extraMLFeatures.mat');
eigenposture_save_filder = fullfile(rootpath, 'CAPTURE', 'myextratsnefeature');
% savefilename ='myextratsnefeature';
% directory_here = pwd;
%feature filename and whether or not to overwrite
savefilename_features = fullfile(roothpath_CAPTURE);
MLmatobjfile = fullfile(savefilename_features, 'myMLfeatures.mat');

%visualize the mocap data
% animate_markers_nonaligned_fullmovie_demo(mocapstruct,1:10:88000, []);

%% Create behavioral features
%this determines the set of frames to use -- in general if the animal is
%resting for too long it will cause errors
mocapstruct.modular_cluster_properties.clipped_index{8} = 1:size(mocapstruct.aligned_mean_position,1 );

% to control the wavelet parameters, you can change the properties in the
% compute_wl_transform_features file
if ~exist(MLmatobjfile,'file') || overwrite_MLmatobjfile
    MLmatobj = create_behavioral_features(mocapstruct,coefficient_file,overwrite_coefficient,linkname);
    save(MLmatobjfile, 'MLmatobj')
else
    disp('Loading ML features')
    MLmatobj = load(MLmatobjfile, 'MLmatobj');
    MLmatobj = MLmatobj.MLmatobj;
end
keyboard
%% perform a tsne embedding subselecting every 50 frames
% analysisparams.tsnegranularity = 50;
% 60,50,30 and 15 look fine
% analysisparams.tsnegranularity = 30;
analysisparams.tsnegranularity = 50;% 120; 250(1.5, corr_thr) ; 100(1.5)- very detailed; 100(2.5)-very good
% Mario: it seems that high number of frames (between 100 - 120) yield better results

%subselect a particular set of features
analysisstruct = compute_tsne_features(MLmatobj,mocapstruct,analysisparams);

%run tsne
disp('%% Running TSNE %%')
zvals = tsne(analysisstruct.jt_features, "Perplexity",30, 'verbose',1); %perplexity 90 works well too (less nr of clusters), but maybe not recommended due to few nr of frames (see length(analysisstruct.jt_features))
% %%%% new
% mouse_keypoints = analysisstruct.jt_features';
% mouse_keypoints_stand = zscore(mouse_keypoints);
% [~,score_features] = pca(mouse_keypoints_stand', 'NumComponents',10);
% score_features = score_features';
% zvals = tsne(score_features', "Perplexity",30, 'Verbose',1, 'NumDimensions',2, 'NumPCAComponents',2);
% %%%%
% [~, zvals] = pca(analysisstruct.jt_features, Centered=false, NumComponents=2);
disp('Done')

figure(1)
plot(zvals(:,1),zvals(:,2),'ob','MarkerFaceColor','b')
title(num2str(analysisparams.tsnegranularity))

% % Plot the manually selected frames
% manual_filename = fullfile(['T:\Marta\test_Formalin\dannce\040822\',animal_ID, 'manual_labeling.xlsx']);
% T = readtable(manual_filename);
% events = [T.frameStart, T.frameStop];
% n_events = sum(~isnan(events(:,1)));
% hold on
% events = [T.frameStart, T.frameStop];
% for  iv = 1:n_events
%     these_events = events(iv,1):events(iv,2);
%     frs = floor(these_events*(1/(analysisparams.tsnegranularity/10)));
%     plot(zvals(frs,1),zvals(frs,2), 'or', 'MarkerFaceColor', 'r')
% 
% end
hold off
%%




analysisstruct.zValues = zvals;

%% clustering parameters
disp('%% INIT clustering %%')
analysisstruct.params.density_res = 1001; %resolution of the map
analysisstruct.params.density_width = 9;%2; %density kernel in tsne space 0.5 for pca; 2.5 tsne .  !! 5 so far works well
analysisstruct.params.expansion_factor = 1.1; %add a little room to the map after kernel smoothing
analysisstruct.params.density_threshold = 1*10^(-5); %remove regions in plots with low density
analysisstruct.matchedconds = {[unique(analysisstruct.condition_inds)]}; %if running over multiple conditions
analysisstruct.conditions_to_run = [unique(analysisstruct.condition_inds)];
analysisstruct.tsnegranularity = analysisparams.tsnegranularity;

params.reorder=1;
analysisstruct = compute_analysis_clusters_demo(analysisstruct,params);
disp('%% Done clustering %%')

%% behavior plots and movies
analysisstruct.conditionnames = ratname;
analysisstruct.ratnames = ratname;
analysisstruct.filesizes = {size(mocapstruct.aligned_mean_position,1 );};

%% plot a tsne map -- see plotting script for parameter definitions
h1=figure(609);
clf;
params.nameplot=1;
params.density_plot =0;
params.watershed = 1;
params.sorted = 1;
params.markersize = 1;
params.coarseboundary =0;
params.do_coarse = 0;
% plot tsne
plot_clustercolored_tsne(analysisstruct,1,params.watershed,h1,params)
set(gcf,'Position',([100 100 1100 1100]))

% bird specific axes
axisparams.zlim = ([200 300]);
axisparams.xlim = ([-400 400]);
axisparams.ylim = ([-400 400]);
%% Plot all clusters
% - Mario:
% Plot the clusters
[cls, c_idx, r] = unique(analysisstruct.annot_reordered{end}, 'stable');

% 
% animate_markers_nonaligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},...
%     find(analysisstruct.annot_reordered{end}==56),[],axisparams);

%
% animate_markers_nonaligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},...
%     find(analysisstruct.annot_reordered{end}==10),[],axisparams);
if do_animation_poses
    h= figure(370);
    clf;

    for ic = 1:numel(cls)
        this_cls = cls(ic);
        fprintf('ic = %i - ', this_cls)
        animate_markers_aligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},...
            find(analysisstruct.annot_reordered{end}==this_cls),h,['cl nr :  ', num2str(this_cls)]);
    end
end


if plot_poses
    % h= figure(370);
    % clf;

    for ic = 1:numel(cls)
        figure
        this_cls = cls(ic);
        fprintf('ic = %i - \n', this_cls)
        plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1},...
            find(analysisstruct.annot_reordered{end}==this_cls),['cl nr :  ', num2str(this_cls)]);
        title(this_cls)
    end
end



%% or use extnded set of 140 features
% if strcmp(ratname,'myrat')
%     % Mario: this might not work properly (check the main code create_extra_behavioral_features, and study it)
%     if ~exist(eigenposture_save_filder, 'dir'), mkdir(eigenposture_save_filder), end
%     MLmatobj_extra = create_extra_behavioral_features(mocapstruct,'myrat',savefilename_extra,overwrite_coefficient,eigenposture_save_filder);
%     % Mario: so far only extracted once saved (modify later)
% %     MLmatobj_extra = load(savefilename_extra);
% %     MLmatobj_extra = MLmatobj_extra.ML_features;
% %     % Mario_ maybe run:
% %     % check the coefficient_file
% %     %load coeff file:
% %     % Try with the appearance coeff
% % 
% %     MLmatobj_extra = compute_appendage_pc_demos(mocapstruct,MLmatobj_extra,coefficient_file,overwrite_coefficient);
% %     %compute the wavelet transform
% %     tic
% %     MLmatobj_extra = compute_wl_transform_features_demo(mocapstruct,MLmatobj_extra,coeff_file,overwrite_coefficient);
% %     toc
% 
%     jt_features_extra = load_extra_tsne_features(mocapstruct,MLmatobj_extra,analysisparams);
% 
%     % look at tsne of these added features
%     zvals_extra = tsne(jt_features_extra);
%     % or the combination
%     %zvals_extra = tsne(cat(2,analysisstruct.jt_features,jt_features_extra));
%     figure(2)
%     plot(zvals_extra(:,1),zvals_extra(:,2),'ob','MarkerFaceColor','b')
%     analysisstruct.zValues_extra = zvals;
% 
% end

%% run sequence and state analysis
params.do_show_pdistmatrix =1;
params.decimation_factor = 1; %downsample if needed to save on memory
% params.decimation_factor = 5; %downsample if needed to save on memory
params.doclustering = 1;

%clustering parameters
params.corr_threshold = 0.2;% 0.2
params.clustercutoff =0.1;%0.12; % 0.65
analysisstruct.plotdirectory = '';
%timescale to use, in seconds
params.timescales = [0.25];%[1./4 2]; 
% params.timescales = [0.1];%[1./4 2]; 
% params.timescales = [0.05];%[1./4 2]; 


analysisstruct.conditionnames = {'test'};
analysisstruct.ratname = {ratname};

hierarchystruct=   find_sequences_states_demo(analysisstruct,1,params);
% identofy the unique cls
[seq_cls, seq_c_idx, seq_r] = unique(hierarchystruct.clustered_behavior{1}, 'stable');
%% the default case

% [maxval] = max(mocapstruct.markers_preproc.SpineM,[],1);
% [minval] =   min(mocapstruct.markers_preproc.SpineM,[],1);
% buffactor_axis = 1.3; %1.1
% buffactor_arena = 1.02;
% 
% 
% th = 0:pi/100:2*pi;
% xcent = (maxval(1)+minval(1))./2;
% ycent = (maxval(1)+minval(1))./2;
% 
% 
% zlimvals = [-50 250];
% xlimvals = [-(304*buffactor_axis- xcent) (304*buffactor_axis+ xcent)];
% ylimvals = [-(304*buffactor_axis- ycent) (304*buffactor_axis+ ycent)];
% axisparams = struct();
% axisparams.xlim = xlimvals;
% axisparams.ylim = ylimvals;
% axisparams.zlim = zlimvals;
h=figure(370);
CL = cell(length(seq_c_idx),1);
for seq_ic = 1:numel(seq_cls)
    this_cls = seq_cls(seq_ic);    fprintf('ic = %i - ', this_cls)
     CL(seq_ic) =  {find(hierarchystruct.clustered_behavior{1}==this_cls)};
    if this_cls==0,  fprintf('\n'),continue, end
    animate_markers_nonaligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},...
        find(hierarchystruct.clustered_behavior{1}==this_cls), h, [], ['ic =  ',num2str(this_cls)]);
   
end
% end
%visualize
%% craete a figure with subplots depicting the clusters





%internal: check dependencies
%[fList,pList] = matlab.codetools.requiredFilesAndProducts('CAPTURE_quickdemo.m');

