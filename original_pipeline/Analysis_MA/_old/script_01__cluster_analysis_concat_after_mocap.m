% This script will concatenate effectively the analysis struct generated
% after running compute_tsne_features. Then will evaluate the clusters by
% the compute_analysis_clusters_demo script.

clear
clc
close all

% Initialize
all_jt_features = [];
all_frames = [];
animal_frames_identifier = cell(0,0);
animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
conditions = {'F', 'F', 'S', 'S', 'N', 'S', 'F', 'F', 'F'};
conditions_idx = zeros(1,length(conditions));
conditions_idx(ismember(conditions, 'F')) = 1;
conditions_idx(ismember(conditions, 'S')) = 2;
conditions_idx(ismember(conditions, 'N')) = 3;

% Initialize aggregated structures
agg_analysisstruct.tsnefeat_name = {};
% agg_analysisstruct.frames_with_good_tracking = cell(1,length(unique(conditions_idx)));
agg_analysisstruct.frames_with_good_tracking = cell(1,length(animal_list));
% agg_analysisstruct.frames_with_good_tracking = cell(1,1);
agg_analysisstruct.frames_tracking_appendages = cell(1,length(animal_list));
agg_analysisstruct.subset_of_points_to_plot_tsne_capped = cell(1,length(animal_list));
agg_analysisstruct.subset_of_points_to_plot_tsne_move = cell(1,length(animal_list));
% agg_analysisstruct.subset_of_points_to_plot_tsne_capped = cell(1,length(unique(conditions_idx)));
% agg_analysisstruct.subset_of_points_to_plot_tsne_move = cell(1,length(unique(conditions_idx)));
agg_analysisstruct.condition_inds = [];
agg_analysisstruct.jt_features = [];
agg_analysisstruct.jt_features_raw = [];
agg_analysisstruct.jt_features_mean = [];  % Assuming you want to concatenate this and take an overall mean later
% agg_analysisstruct.jt_features_std = [];   % Same assumption as mean
agg_analysisstruct.file_sizes = cell(1,length(animal_list));
agg_analysisstruct.mocapstruct_reduced_agg = cell(1,length(animal_list));
% agg_analysisstruct.file_sizes = cell(1,length(unique(conditions_idx)));
% agg_analysisstruct.mocapstruct_reduced_agg = cell(1,length(unique(conditions_idx)));
agg_analysisstruct.coarse_annotation_mat = [];
analysisparams.tsnegranularity = 50;% 120; 250(1.5, corr_thr) ; 100(1.5)- very detailed; 100(2.5)-very good

% Current offset, starts at 0
offset = 0;
AF = []; % actual frames to tracK
for k = 1:length(animal_list)
    animal_name = animal_list{k};
    disp(animal_name)
    % Load MLmatobj and ratception_struct for the current animal
    load(['D:\test_CAPTURE\',animal_name,'\CAPTURE\myMLfeatures.mat']);
    load(['D:\test_CAPTURE\',animal_name,'\ratception_prediction.mat']);

    mocapstruct = ratception_struct;
    mocapstruct.modular_cluster_properties.clipped_index{8} = 1:size(mocapstruct.aligned_mean_position,1);

    % analysisparams.tsnegranularity = 50;
    analysisstruct = compute_tsne_features(MLmatobj, mocapstruct, analysisparams);

    % Aggregate each field of analysisstruct
    fields = fieldnames(agg_analysisstruct);
    for f = 1:length(fields)
        thisf = fields{f};

        % Handle special cases
        switch thisf
            case 'tsnefeat_name'
                agg_analysisstruct.tsnefeat_name = analysisstruct.tsnefeat_name;
                continue;  % Skip to next iteration

            % case 'mocapstruct_reduced_agg'
            %     % For now, simply copying the new data; we will handle concatenation next
            %     if k==1
            %         agg_analysisstruct.mocapstruct_reduced_agg = analysisstruct.mocapstruct_reduced_agg;
            %     end
            %     % Extracting body parts and concatenating
            %     body_parts_preproc = fieldnames(analysisstruct.mocapstruct_reduced_agg{1}.markers_aligned_preproc);
            %     body_parts_not_aligned = fieldnames(analysisstruct.mocapstruct_reduced_agg{1}.markers_preproc);
            %     for bp = 1:length(body_parts_preproc)
            %         if k==1
            %             agg_analysisstruct.mocapstruct_reduced_agg{1}.markers_aligned_preproc.(body_parts_preproc{bp}) =[];
            %         end
            %         agg_analysisstruct.mocapstruct_reduced_agg{1}.markers_aligned_preproc.(body_parts_preproc{bp}) = ...
            %             cat(1, agg_analysisstruct.mocapstruct_reduced_agg{1}.markers_aligned_preproc.(body_parts_preproc{bp}), ...
            %             analysisstruct.mocapstruct_reduced_agg{1}.markers_aligned_preproc.(body_parts_preproc{bp}));
            %     end
            %     for bp = 1:length(body_parts_not_aligned)
            %         if k==1
            %             agg_analysisstruct.mocapstruct_reduced_agg{1}.markers_preproc.(body_parts_not_aligned{bp}) = [];
            %         end
            %         agg_analysisstruct.mocapstruct_reduced_agg{1}.markers_preproc.(body_parts_not_aligned{bp}) = ...
            %             cat(1, agg_analysisstruct.mocapstruct_reduced_agg{1}.markers_preproc.(body_parts_not_aligned{bp}), ...
            %             analysisstruct.mocapstruct_reduced_agg{1}.markers_preproc.(body_parts_not_aligned{bp}));
            %     end
                % continue;  % Skip to next iteration
            case 'condition_inds'
                 % for using per conditions
                 % agg_analysisstruct.condition_inds =[agg_analysisstruct.condition_inds, analysisstruct.condition_inds * conditions_idx(k)];
                 % for using per animal
                 agg_analysisstruct.condition_inds =[agg_analysisstruct.condition_inds, analysisstruct.condition_inds * (k)];

                % if k == 1
                %     agg_analysisstruct.condition_inds = analysisstruct.condition_inds;
                % else
                %     agg_analysisstruct.condition_inds = cat(2,agg_analysisstruct.condition_inds, analysisstruct.condition_inds +k-1);
                % end
                continue
            case 'subset_of_points_to_plot_tsne_capped'
                % this_condition = conditions_idx(k);
                % if isempty(agg_analysisstruct.subset_of_points_to_plot_tsne_capped{1,this_condition})
                %     agg_analysisstruct.subset_of_points_to_plot_tsne_capped(1,this_condition) = analysisstruct.subset_of_points_to_plot_tsne_capped;
                % else
                %     agg_analysisstruct.subset_of_points_to_plot_tsne_capped(1,this_condition) =...
                %         {cat(2, agg_analysisstruct.subset_of_points_to_plot_tsne_capped{1,this_condition},...
                %         analysisstruct.subset_of_points_to_plot_tsne_capped{1})} ; % + max(agg_analysisstruct.subset_of_points_to_plot_tsne_capped{1,this_condition}))}
                % end
                % [agg_analysisstruct.subset_of_points_to_plot_tsne_capped{1,this_condition}, analysisstruct.subset_of_points_to_plot_tsne_capped{1}];
                agg_analysisstruct.subset_of_points_to_plot_tsne_capped(1,k) = analysisstruct.subset_of_points_to_plot_tsne_capped;
                % Works for individual animals
                continue
            case 'subset_of_points_to_plot_tsne_move'
                % this_condition = conditions_idx(k);
                % if isempty(agg_analysisstruct.subset_of_points_to_plot_tsne_move{1,this_condition})
                %     agg_analysisstruct.subset_of_points_to_plot_tsne_move(1,this_condition) = analysisstruct.subset_of_points_to_plot_tsne_move;
                % else
                %     agg_analysisstruct.subset_of_points_to_plot_tsne_move(1,this_condition) =...
                %         {cat(2, agg_analysisstruct.subset_of_points_to_plot_tsne_move{1,this_condition},...
                %         analysisstruct.subset_of_points_to_plot_tsne_move{1})} ; % + max(agg_analysisstruct.subset_of_points_to_plot_tsne_move{1,this_condition}))};
                % end
                % [agg_analysisstruct.subset_of_points_to_plot_tsne_capped{1,this_condition},
                % analysisstruct.subset_of_points_to_plot_tsne_capped{1}];
                % Does not work

                agg_analysisstruct.subset_of_points_to_plot_tsne_move(1,k) = analysisstruct.subset_of_points_to_plot_tsne_move;
                % Works for individual animals
                continue
            case 'mocapstruct_reduced_agg'
                % keyboard
                % this_condition = conditions_idx(k);
                % if isempty(agg_analysisstruct.mocapstruct_reduced_agg{1,this_condition})
                %     agg_analysisstruct.mocapstruct_reduced_agg(1,this_condition) = analysisstruct.mocapstruct_reduced_agg;
                % else
                %     agg_analysisstruct.mocapstruct_reduced_agg(1,this_condition) =...
                %         {cat(2, agg_analysisstruct.mocapstruct_reduced_agg{1,this_condition},...
                %         analysisstruct.mocapstruct_reduced_agg{1} + max(agg_analysisstruct.mocapstruct_reduced_agg{1,this_condition}))};
                % end

                % Second attempt: concatenating by condition
                % this_condition = conditions_idx(k);
                % if isempty(agg_analysisstruct.mocapstruct_reduced_agg{1,this_condition})
                %     agg_analysisstruct.mocapstruct_reduced_agg(1,this_condition) = analysisstruct.mocapstruct_reduced_agg;
                %     % end
                %     % Extracting body parts and concatenating
                % else
                %     body_parts_preproc = fieldnames(analysisstruct.mocapstruct_reduced_agg{1}.markers_aligned_preproc);
                %     body_parts_not_aligned = fieldnames(analysisstruct.mocapstruct_reduced_agg{1}.markers_preproc);
                %     for bp = 1:length(body_parts_preproc)
                %         % if k==1
                %         % agg_analysisstruct.mocapstruct_reduced_agg{1}.markers_aligned_preproc.(body_parts_preproc{bp}) =[];
                %         % end
                %         agg_analysisstruct.mocapstruct_reduced_agg{this_condition}.markers_aligned_preproc.(body_parts_preproc{bp}) = ...
                %             cat(1, agg_analysisstruct.mocapstruct_reduced_agg{this_condition}.markers_aligned_preproc.(body_parts_preproc{bp}), ...
                %             analysisstruct.mocapstruct_reduced_agg{1}.markers_aligned_preproc.(body_parts_preproc{bp}));
                %     end
                %     for bp = 1:length(body_parts_not_aligned)
                %         % if k==1
                %         % agg_analysisstruct.mocapstruct_reduced_agg{1}.markers_preproc.(body_parts_not_aligned{bp}) = [];
                %         % end
                %         agg_analysisstruct.mocapstruct_reduced_agg{this_condition}.markers_preproc.(body_parts_not_aligned{bp}) = ...
                %             cat(1, agg_analysisstruct.mocapstruct_reduced_agg{this_condition}.markers_preproc.(body_parts_not_aligned{bp}), ...
                %             analysisstruct.mocapstruct_reduced_agg{1}.markers_preproc.(body_parts_not_aligned{bp}));
                %     end
                % end
                %  % concatenating by animal
                agg_analysisstruct.mocapstruct_reduced_agg(1,k) = analysisstruct.mocapstruct_reduced_agg;
                continue

            case 'frames_with_good_tracking'
                % for conditions:
                %  this_condition = conditions_idx(k);
                % if isempty(agg_analysisstruct.frames_with_good_tracking{1,this_condition})
                %     agg_analysisstruct.frames_with_good_tracking(1,this_condition) = analysisstruct.frames_with_good_tracking;
                % else
                %     agg_analysisstruct.frames_with_good_tracking(1,this_condition) =...
                %         {cat(1, agg_analysisstruct.frames_with_good_tracking{1,this_condition},...
                %         analysisstruct.frames_with_good_tracking{1} + max(agg_analysisstruct.frames_with_good_tracking{1,this_condition}))};
                % end

                % For animals
                agg_analysisstruct.frames_with_good_tracking(1,k) = analysisstruct.frames_with_good_tracking;

                continue
        case'frames_tracking_appendages'
             agg_analysisstruct.frames_tracking_appendages(1,k) = {analysisstruct.frames_tracking_appendages};
             continue



        end

        % Determine concatenation dimension based on the size of the original field
        if iscell(analysisstruct.(fields{f}))
            is_cell = 1;
        else
            is_cell=0;
        end


        if size(analysisstruct.(fields{f}), 2) > size(analysisstruct.(fields{f}), 1)
            concat_dim = 2;  % V concatenation
        else
            if is_cell
                if size(cell2mat(analysisstruct.(fields{f})),1) >1

                    concat_dim = 1;  % H concatenation
                else
                    concat_dim = 2;
                end
            else
                concat_dim = 1;
            end


        end

        % if iscell(analysisstruct.(fields{f}))
        %     agg_analysisstruct.(fields{f}) = cat(concat_dim, agg_analysisstruct.(fields{f}), analysisstruct.(fields{f}));
        % else
        %     agg_analysisstruct.(fields{f}) = cat(concat_dim, agg_analysisstruct.(fields{f}), analysisstruct.(fields{f}));
        % end

        % If it's a cell but the content is numeric, handle the content
        if iscell(analysisstruct.(fields{f})) && isnumeric(analysisstruct.(fields{f}){1})
            agg_analysisstruct.(fields{f}){1} = cat(concat_dim, agg_analysisstruct.(fields{f}){1}, analysisstruct.(fields{f}){1});

            % If it's just a numeric array, concatenate as normal
        elseif ~iscell(analysisstruct.(fields{f}))
            agg_analysisstruct.(fields{f}) = cat(concat_dim, agg_analysisstruct.(fields{f}), analysisstruct.(fields{f}));

            % For other cell types, concatenate the cells themselves
        else
            agg_analysisstruct.(fields{f}) = cat(concat_dim, agg_analysisstruct.(fields{f}), analysisstruct.(fields{f}));
        end
    end


    % Scenario 2: Using frames_with_good_tracking
    adjusted_frames = cell2mat(analysisstruct.frames_with_good_tracking) + offset;
    all_frames = [all_frames; adjusted_frames];

    % Append the animal_name to the identifier list for the current batch of frames
    animal_frames_identifier = [animal_frames_identifier; repmat({animal_name}, length(adjusted_frames), 1)];

    % Update the offset for next iteration
    offset = adjusted_frames(end);

    actual_frames = analysisstruct.frames_tracking_appendages;
    AF = [AF; actual_frames];

    % clear analysisstruct
end

T =table(animal_frames_identifier,AF, 'VariableNames',{'animal_ID', 'Actual_frame'});
agg_analysisstruct.frame_data = T;

% At the end, all_jt_features contains the concatenated features and all_frames has the adjusted frame IDs.
disp('%% Running TSNE %%')
zvals = tsne(agg_analysisstruct.jt_features, "Perplexity",75, 'verbose',1, 'NumDimensions',3);

figure, gscatter(zvals(:,1), zvals(:,2), animal_frames_identifier)
if size(zvals,2) > 2
    zvals = zvals(:,1:2);
end
agg_analysisstruct.zValues = zvals;

%% clustering parameters
disp('%% INIT clustering %%')
agg_analysisstruct.params.density_res = 1001; %resolution of the map
agg_analysisstruct.params.density_width = 5;%3;9,2; %density kernel in tsne space 0.5 for pca; 2.5 tsne .  !! 5 so far works well
agg_analysisstruct.params.expansion_factor = 1.1; %add a little room to the map after kernel smoothing
agg_analysisstruct.params.density_threshold = 1*10^(-5); %remove regions in plots with low density
agg_analysisstruct.matchedconds = {[unique(agg_analysisstruct.condition_inds)]}; %if running over multiple conditions
agg_analysisstruct.conditions_to_run = [unique(agg_analysisstruct.condition_inds)];
agg_analysisstruct.tsnegranularity = analysisparams.tsnegranularity;

params.reorder=1;
agg_analysisstruct = compute_analysis_clusters_demo(agg_analysisstruct,params);
disp('%% Done clustering %%')
ratname ='myrat';
%% behavior plots and movies
agg_analysisstruct.conditionnames = ratname;
agg_analysisstruct.ratnames = ratname;
agg_analysisstruct.filesizes = {size(mocapstruct.aligned_mean_position,1 );};
%% save
save('H:\DANNCE\_CAPTURE_results_tests\agg_analysisstruct.mat', 'agg_analysisstruct');
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
plot_clustercolored_tsne(agg_analysisstruct,1,params.watershed,h1,params)
set(gcf,'Position',([100 100 1100 1100]))


%% in case we run extra features
% savefilename_extra = fullfile('D:\test_CAPTURE', 'myextratsnefeature', 'extraMLFeatures.mat');
% eigenposture_save_filder = fullfile('D:\test_CAPTURE', 'myextratsnefeature');
% if ~exist(eigenposture_save_filder, 'dir'), mkdir(eigenposture_save_filder), end
% overwrite_coefficient = 1;
% MLmatobj_extra = create_extra_behavioral_features(mocapstruct,'myrat',savefilename_extra,overwrite_coefficient,eigenposture_save_filder);
% 
% jt_features_extra = load_extra_tsne_features(mocapstruct,MLmatobj_extra,analysisparams);
% 
% % look at tsne of these added features
% zvals_extra = tsne(jt_features_extra);
% % or the combination
% %zvals_extra = tsne(cat(2,analysisstruct.jt_features,jt_features_extra));
% figure(2)
% plot(zvals_extra(:,1),zvals_extra(:,2),'ob','MarkerFaceColor','b')
% analysisstruct.zValues_extra = zvals;
%%
[cls, c_idx, r] = unique(agg_analysisstruct.annot_reordered{end}, 'stable');
% [cls, c_idx, r] = unique(agg_analysisstruct.annot_reordered{end});
% so far this doesnt work anymore, since the mocap struct is organized per
% animal
plot_poses = 1;
if plot_poses
    % h= figure(370);
    % clf;

    figure('pos', [10,300,1500,1900])
    nclus = numel(cls);
    n_rows = ceil(sqrt(nclus));
    n_cols = ceil(sqrt(nclus));
    for ic = 1:numel(cls)
        subplot(n_rows, n_cols, ic)
        this_cls = cls(ic);
        fprintf('ic = %i - \n', this_cls)
        plot_mean_cluster_aligned(agg_analysisstruct.mocapstruct_reduced_agg{1},...
            find(agg_analysisstruct.annot_reordered{end}==this_cls),['cl nr :  ', num2str(this_cls)]);
        title(this_cls)
    end
end


%% run sequence and state analysis
params.do_show_pdistmatrix =1;
params.decimation_factor = 3; %downsample if needed to save on memory
params.doclustering = 1;

%clustering parameters
params.corr_threshold = 0.3;
params.clustercutoff = 0.65;
agg_analysisstruct.plotdirectory = '';
%timescale to use, in seconds (in minutes, effectively)
params.timescales = [1./4 2];

agg_analysisstruct.conditionnames = animal_list;
% agg_analysisstruct.conditionnames = {'F', 'S', 'N'};
agg_analysisstruct.ratname = ratname;

%%
hierarchystruct=   find_sequences_states_demo(agg_analysisstruct,([1:9]),params);
%% Loop through animals and save the hierachy structure
for iid = 1:numel(animal_list) % so far this is done per animal. But at some point it will be done per condition
    filename_hierarchy = fullfile('D:\test_CAPTURE',animal_name,'CAPTURE', 'hierarchystruct.mat');
    if ~exist(filename_hierarchy, 'file') || overwrite_hierarchy
        hierarchystruct=   find_sequences_states_demo(agg_analysisstruct,iid,params);
        save(filename_hierarchy,'filename_hierarchy')
        % else
        % load(filename_hierarchy)
    end

end
%%
% Assuming you have a list 'cluster_sequence' that denotes the cluster (state) at each frame
% For example: cluster_sequence = [1, 2, 2, 1, 3, 2, ...];
% choose the timescale (1 or 2)
cluster_sequence = hierarchystruct.clustered_behavior{2};
cluster_sequence = cluster_sequence +1; % nomenclature starts from 0

unique_clusters = unique(cluster_sequence);
n_unique_clusters = length(unique_clusters);
cluster_to_index = containers.Map('KeyType', 'double', 'ValueType', 'double');

for i = 1:n_unique_clusters
    cluster_to_index(unique_clusters(i)) = i;
end

transition_matrix = zeros(n_unique_clusters, n_unique_clusters);

for i = 1:length(cluster_sequence) - 1
    current_cluster = cluster_sequence(i);
    next_cluster = cluster_sequence(i+1);
    
    current_index = cluster_to_index(current_cluster);
    next_index = cluster_to_index(next_cluster);
    
    transition_matrix(current_index, next_index) = transition_matrix(current_index, next_index) + 1;
end

row_sums = sum(transition_matrix, 2);
transition_matrix = bsxfun(@rdivide, transition_matrix, row_sums);

% Normalize each row to sum to 1
% row_sums = sum(transition_matrix, 2);
% transition_matrix = bsxfun(@rdivide, transition_matrix, row_sums);

% Display the transition matrix
disp('Transition Matrix:');
% disp(transition_matrix);
figure
% If you want to visualize this matrix
imagesc(transition_matrix);
colorbar;
xlabel('Next Cluster');
ylabel('Current Cluster');
title('Markovian Transition Matrix');

script_04_test_transition_matrix_analysis