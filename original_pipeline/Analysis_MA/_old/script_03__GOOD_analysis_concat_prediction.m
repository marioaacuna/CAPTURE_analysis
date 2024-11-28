%% Add pramble


%% preINIT
close all
clc
clear
%% INIT
overwrite_pred_concat = 0; % in case you want to add more animals, set this to 1
overwrite_ratception = 0; % will do or not mocap
overwrite_MLmatobjfile = 0; % Overwrite extarcted features
overwrite_coefficient = 0; % Overwrite dim red coeffs 
overwrite_zvals = 0; % Overwrite tsne zvals
do_extra_features = 0;

init_frame_rate = 100; % effective frame rate of videos

if ispc
    % so far for these test, we have the data on local drive D
    rootpath = 'D:\CAPTURE';
else
    rootpath = '/Users/mario/Library/Mobile Documents/com~apple~CloudDocs/_todo/SpontaneousPainProject/240704_data'; % I have updated where the files are to gather them more quickly
end

if ~exist("rootpath", 'dir'), mkdir(rootpath); end
 
% animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
% conditions = {'F', 'F', 'S', 'S', 'N', 'S', 'F', 'F', 'F'};

% Define the animals and conditions
animals = {'AK_552', 'AK_553', 'AK_665', 'AK_667'}; % animal ids
conditions = {'0', '1'}; % 0 for condition 'S', and 1 for condition 'F'
dannce_path = 'D:\DANNCE';
% Initialize structures
agg_predictions = struct();
animal_condition_identifier = {};


analysis_filename = fullfile(rootpath,'raw_concat_analysis.mat' );
filename_predictions = fullfile(rootpath, 'agg_predictions.mat');

if ~any(exist(filename_predictions, 'file'))
    run_pred_concat = 1 ;
else
    run_pred_concat = 0;
end

if run_pred_concat || overwrite_pred_concat
    disp('%% Running concatenation')
    % Initialize the aggregate predictions structure
    agg_predictions = struct();
    animal_condition_identifier = {};
  
    % Current offset for tracking frames
    offset = 0;


    % Loop through each animal
    for i = 1:length(animals)
        animal_id = animals{i};

        % Loop through each condition for the current animal
        for j = 1:length(conditions)
            condition = conditions{j};

            % Construct the path to the predictions file
            load_path = fullfile(dannce_path, animal_id, condition,'DANNCE\predict_results', 'predictions.mat');

            % Check if the file exists
            if ~exist(load_path, 'file')
                warning('File not found: %s', load_path);
                continue;
            end

            % Load predictions for the current animal and condition
            load(load_path); % Assuming the loaded structure is called 'predictions'

            body_parts = fieldnames(predictions); % Get the list of body parts

            % Process each body part
            for k = 1:length(body_parts)
                part_name = body_parts{k};

                if isfield(agg_predictions, part_name)
                    if ~strcmp(part_name, 'sampleID')
                        % Concatenate if the field already exists in agg_predictions
                        agg_predictions.(part_name) = cat(1, agg_predictions.(part_name), predictions.(part_name));
                    else
                        % For sampleID, adjust the values before concatenating
                        last_sample_id = 0;
                        if ~isempty(agg_predictions.(part_name))
                            last_sample_id = agg_predictions.(part_name)(end);
                        end
                        agg_predictions.(part_name) = cat(2, agg_predictions.(part_name), predictions.(part_name) + last_sample_id);
                    end
                else
                    % Otherwise, initialize the field in agg_predictions
                    agg_predictions.(part_name) = predictions.(part_name);
                end
            end

            % Update the animal_condition_identifier
            num_frames_for_current_combo = size(predictions.(body_parts{1}), 1);
            condition_letter = condition_to_letter(condition);
            animal_condition_identifier = [animal_condition_identifier; repmat({sprintf('%s_%s', animal_id, condition_letter)}, num_frames_for_current_combo, 1)];
        end
    end
    % At this point, agg_predictions contains concatenated body part data for all animals
    % and animal_frames_identifier indicates the animal for each frame.
    
    % save the data
 
    predictions = agg_predictions;
    % predictions.predictions = agg_predictions;
    save(filename_predictions, 'predictions', 'animal_condition_identifier')
else
    disp('predictions previously concatenated, now loading them')
    load(filename_predictions)
end

%%
long_animal_frames_identifier = repelem(animal_condition_identifier,3);
%% INIT ratception procedure
filename_ratception = fullfile(rootpath, "ratception_prediction.mat");

animal_name = 'mario_mouse22';

input_params = struct();
input_params.SpineF_marker = 'SpineF';
input_params.SpineM_marker = 'SpineM';
% input_params.repfactor = 300/30;
input_params.repfactor = round(300/init_frame_rate);
input_params.conversion_factor = 1;

% Run prepro if it doesn't exist
if ~exist(filename_ratception, "file") || overwrite_ratception
    disp('%% Running Pre-Pro %% ')
    ratception_struct = preprocess_dannce(filename_predictions,filename_ratception,animal_name,input_params);
else
    disp('Loading previously analysied prepro data')
   load(filename_ratception);
end

%%
roothpath_CAPTURE = fullfile(rootpath, 'CAPTURE');

coefficient_file = fullfile(roothpath_CAPTURE,'coefficients.mat');
linkname = 'mario_mouse22';
ratname ='myrat';% 'test_mouse';

%% Load Mocapstruct
mocapstruct = ratception_struct;
if do_extra_features
    % In case you want to do some extra features
    savefilename_extra = fullfile(rootpath, 'CAPTURE', 'myextratsnefeature', 'extraMLFeatures.mat');
    eigenposture_save_filder = fullfile(rootpath, 'CAPTURE', 'myextratsnefeature');
end
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

if ~exist(fullfile(rootpath, 'CAPTURE'),'dir' )
    mkdir(fullfile(rootpath, 'CAPTURE'))
end

if ~exist(MLmatobjfile,'file') || overwrite_MLmatobjfile
    MLmatobj = create_behavioral_features(mocapstruct,coefficient_file,overwrite_coefficient,linkname);
    save(MLmatobjfile, 'MLmatobj', '-v7.3')
else
    disp('Loading ML features')
    MLmatobj = load(MLmatobjfile, 'MLmatobj');
    MLmatobj = MLmatobj.MLmatobj;
end

%%
analysisparams.tsnegranularity = 25;% 50:default ;120; 250(1.5, corr_thr) ; 100(1.5)- very detailed; 100(2.5)-very good
% Mario: it seems that high number of frames (between 100 - 120) yield better results

%subselect a particular set of features
analysisstruct = compute_tsne_features(MLmatobj,mocapstruct,analysisparams);

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

% analysisstruct.condition_inds = cond_inds;

% conditions ie. actua conditions, concatenating animals in same cond
condition_inds = cond_inds; % sorting per condition
for iff = 1:length(animal_list_used_after_analysis)
     animal_ID = animal_list_used_after_analysis{iff};
    condition_inds(iff) = double(endsWith(animal_ID, '_F') +1);
   
end
% analysisstruct.condition_inds = condition_inds;
% this does not work, see  compute_tsne_features.m line 144

%% Plot dim red
zvals_filename = fullfile(roothpath_CAPTURE, 'zvals.mat');
perplexity = 200; % 75;

if ~exist(zvals_filename, 'file') || overwrite_zvals

    %run tsne
    disp('%% Running TSNE %%')
    zvals = tsne(analysisstruct.jt_features, "Perplexity",perplexity, 'Exaggeration', 20,'verbose',1,'LearnRate', 1200); %perplexity 90 works well too (less nr of clusters), but maybe not recommended due to few nr of frames (see length(analysisstruct.jt_features))
    % save zvals to then read later if necessary
    save(zvals_filename, 'zvals','-mat')

    % Plot TSNE per animal
    figure(1)
    gscatter(zvals(:,1), zvals(:,2), cond_inds)
    % plot(zvals(:,1),zvals(:,2),'ob','MarkerFaceColor','b', 'MarkerSize',2)
    title({['Granu: ',num2str(analysisparams.tsnegranularity)], ['Perp: ', num2str(perplexity)]})
    set(gcf,'Position',([100 100 1100 1100]))
    set(gcf, 'color', 'w')
else
    disp(' Loading TSNE zvals')
    load(zvals_filename)
end

disp('Done TSNE')

%% clustering parameters

disp('%% INIT clustering %%')
analysisstruct.zValues = zvals;
analysisstruct.params.density_res = 1001; %resolution of the map
analysisstruct.params.density_width = 1.5;%default:3; 2; %density kernel in tsne space 0.5 for pca; 2.5 tsne .  !! 5 so far works well
analysisstruct.params.expansion_factor = 1.1; %add a little room to the map after kernel smoothing
analysisstruct.params.density_threshold = 1*10^(-5); %remove regions in plots with low density
analysisstruct.matchedconds = {[unique(analysisstruct.condition_inds)]}; %if running over multiple conditions
analysisstruct.conditions_to_run = [unique(analysisstruct.condition_inds)];
analysisstruct.tsnegranularity = analysisparams.tsnegranularity;

params.reorder=1;
analysisstruct = compute_analysis_clusters_demo(analysisstruct,params); % check line 248, cluster_tsne_map.m
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
set(h1,'Position',([100 100 1100 1100]))
  
% bird specific axes
axisparams.zlim = ([200 300]);
axisparams.xlim = ([-400 400]);
axisparams.ylim = ([-400 400]);

cluster_figure_filename = fullfile(roothpath_CAPTURE, 'Tsne_clusters.pdf');
export_fig(cluster_figure_filename, '-pdf', h1)

%% Plot cluster poses
[cls, c_idx, r] = unique(analysisstruct.annot_reordered{end}, 'stable');

plot_poses = 1;
if plot_poses
    % h= figure(370);
    % clf;

    fig_poses = figure('pos', [10,300,1500,1900]);
    nclus = numel(cls);
    n_rows = ceil(sqrt(nclus));
    n_cols = ceil(sqrt(nclus));
    for ic = 1:numel(cls)
        subplot(n_rows, n_cols, ic)
        this_cls = cls(ic);
        fprintf('ic = %i - \n', this_cls)
        plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1},...
            find(analysisstruct.annot_reordered{end}==this_cls),['cl nr :  ', num2str(this_cls)]);
        title(this_cls)
    end
end
% save cluster plot
cluster_poses_figure_filename = fullfile(roothpath_CAPTURE, 'Poses_clusters.pdf');
export_fig(cluster_poses_figure_filename, '-pdf', fig_poses)

%% run sequence and state analysis
params.do_show_pdistmatrix =1;
params.decimation_factor = 1; %downsample if needed to save on memory
% params.decimation_factor = 5; %downsample if needed to save on memory
params.doclustering = 1;

%clustering parameters
params.corr_threshold = 0.2;% 0.2
params.clustercutoff =0.65;%0.12; % 0.65
analysisstruct.plotdirectory = '';
%timescale to use, in seconds
params.timescales = [1./4 2];  % this is in minutes, check the find_sequences_state_demo code 

analysisstruct.conditionnames = {'test'};
analysisstruct.ratname = {ratname};
condition =1;

% Uncomment the following lines to to sequence analysis
% hierarchystruct=   find_sequences_states_demo(analysisstruct,condition,params);
% % identofy the unique cls
% [seq_cls, seq_c_idx, seq_r] = unique(hierarchystruct.clustered_behavior{1}, 'stable');

% h=figure(370);
% CL = cell(length(seq_c_idx),1);
% for seq_ic = 1:numel(seq_cls)
%     this_cls = seq_cls(seq_ic);    fprintf('ic = %i - ', this_cls)
%      CL(seq_ic) =  {find(hierarchystruct.clustered_behavior{1}==this_cls)};
%     if this_cls==0,  fprintf('\n'),continue, end
%     animate_markers_nonaligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},...
%         find(hierarchystruct.clustered_behavior{1}==this_cls), h, [], ['ic =  ',num2str(this_cls)]);
% 
% end



%% save analysis
save(analysis_filename, 'analysisstruct' , '-v7.3')


%% This is the end of the analysis script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Functions
function letter = condition_to_letter(condition)
    if strcmp(condition, '0')
        letter = 'S';
    elseif strcmp(condition, '1')
        letter = 'F';
    else
        error('Unknown condition: %s', condition);
    end
end