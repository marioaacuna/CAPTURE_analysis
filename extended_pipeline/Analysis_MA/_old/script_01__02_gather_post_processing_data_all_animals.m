%% preamble
% This script will take each animal's predictions from dannce and run the mocap and
% MLmatobj analysis for the subsequent hierarchy analysis. Note that the
% clustering of posses is done in script_03, where all animals are grouped together.

%% INIT
close all
clc
clear
% List of animals

overwrite_ratception = 0; % will do or not mocap
overwrite_MLmatobjfile = 0; % Overwrite extarcted features
overwrite_coefficient = 0; % Overwrite dim red coeffs 

% so far for these test, we have the data on local drive D
rootpath = 'D:\test_CAPTURE';

animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
conditions = {'F', 'F', 'S', 'S', 'N', 'S', 'F', 'F', 'F'};

% Input variables
input_params = struct();
input_params.SpineF_marker = 'SpineF';
input_params.SpineM_marker = 'SpineM';
input_params.repfactor = round(300/120);
ratname ='myrat';% 'test_mouse';
animal_name = 'mario_mouse22';


%Input for ML features
linkname = 'mario_mouse22';
ratname ='myrat';% 'test_mouse';

% inputs for clustering
density_res = 1001; %resolution of the map
density_width = 3;%2; %density kernel in tsne space 0.5 for pca; 2.5 tsne .  !! 5 so far works well
expansion_factor = 1.1; %add a little room to the map after kernel smoothing
density_threshold = 1*10^(-5); %remove regions in plots with low density



%% Loop through animals
for iid = 1:length(animal_list)
    animal_ID = animal_list{iid};
    fprintf('%%%% ANALYSING %s %%%%', animal_ID)
    roothpath_CAPTURE = fullfile(rootpath, animal_ID, 'CAPTURE');
    filename_predictions = fullfile(rootpath, animal_ID, "predictions.mat");
    filename_ratception = fullfile(rootpath,animal_ID, "ratception_prediction.mat");
    filename_hierarchy = fullfile(roothpath_CAPTURE, "hierarchystruct.mat");
    filename_analysis = fullfile(roothpath_CAPTURE, "analysissctruct.mat");
    % Load/run preprocessing
    if ~exist(filename_ratception, 'file') || overwrite_ratception
        disp('%% Running Pre-Pro %% ')
        ratception_struct = preprocess_dannce(filename_predictions,filename_ratception,animal_name,input_params);
    else
        disp('Loading previously analysied prepro data')
        load(filename_ratception);
    end
    % Rename Mocapstruct
    mocapstruct = ratception_struct;
    clear ratception_struct
    % Load/run ML features
    if ~exist(roothpath_CAPTURE, 'dir'), mkdir(roothpath_CAPTURE); end
    coefficient_file = fullfile(roothpath_CAPTURE,'coefficients.mat');
    linkname = 'mario_mouse22';
    ratname ='myrat';% 'test_mouse';

    % In case you want to do some extra features
    savefilename_extra = fullfile(rootpath, animal_ID,'CAPTURE', 'myextratsnefeature', 'extraMLFeatures.mat');
    eigenposture_save_filder = fullfile(rootpath,animal_ID, 'CAPTURE', 'myextratsnefeature');
    if  ~exist(savefilename_extra,'file'), mkdir(savefilename_extra), end
    if  ~exist(eigenposture_save_filder,'file'), mkdir(eigenposture_save_filder), end

    %feature filename and whether or not to overwrite
    savefilename_features = fullfile(roothpath_CAPTURE);
    MLmatobjfile = fullfile(savefilename_features, 'myMLfeatures.mat');
    
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
    % % MLmatobj_extra = create_extra_behavioral_features(mocapstruct,'myrat',savefilename_extra,overwrite_coefficient,eigenposture_save_filder);
    % Perform tsne embedding subselecting every 50 frames
    analysisparams.tsnegranularity = 50;% 120; 250(1.5, corr_thr) ; 100(1.5)- very detailed; 100(2.5)-very good
    analysisstruct = compute_tsne_features(MLmatobj,mocapstruct,analysisparams);
    zvals = tsne(analysisstruct.jt_features, "Perplexity",30, 'verbose',0); %perplexity 90 works well too (less nr of clusters), but maybe not recommended due to few nr of frames (see length(analysisstruct.jt_features))
    analysisstruct.zValues = zvals;
    % Run clustering
    analysisstruct.params.density_res = density_res; %resolution of the map
    analysisstruct.params.density_width = density_width;%2; %density kernel in tsne space 0.5 for pca; 2.5 tsne .  !! 5 so far works well
    analysisstruct.params.expansion_factor = expansion_factor; %add a little room to the map after kernel smoothing
    analysisstruct.params.density_threshold = density_threshold; %remove regions in plots with low density
    analysisstruct.matchedconds = {[unique(analysisstruct.condition_inds)]}; %if running over multiple conditions
    analysisstruct.conditions_to_run = [unique(analysisstruct.condition_inds)];
    analysisstruct.tsnegranularity = analysisparams.tsnegranularity;

    params.reorder=1;
    analysisstruct = compute_analysis_clusters_demo(analysisstruct,params);
    
    analysisstruct.conditionnames = ratname;
    analysisstruct.ratnames = ratname;
    analysisstruct.filesizes = {size(mocapstruct.aligned_mean_position,1 );};

    % Run hierarchy structure
    params.do_show_pdistmatrix = 0;
    params.decimation_factor = 1; %downsample if needed to save on memory
    % params.decimation_factor = 5; %downsample if needed to save on memory
    params.doclustering = 1;
    %clustering parameters
    params.corr_threshold = 0.3;% 0.2, paper is 0.3
    params.clustercutoff =0.65;%0.12; % 0.65
    analysisstruct.plotdirectory = '';
    %timescale to use, in seconds
    params.timescales = [1./4 2];  % this is in minutes, check the find_sequences_state_demo code
    % paper 15 s and 12 s -> 1/4, 2
   
    analysisstruct.conditionnames = {'test'};
    analysisstruct.ratname = {ratname};
    % Save analysis structure
    save(filename_analysis, "analysisstruct")
    condition =1;
    hierarchystruct=   find_sequences_states_demo(analysisstruct,condition,params);
    save(filename_hierarchy, "hierarchystruct")
end
keyboard




