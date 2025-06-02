% load the analysisstructure
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
mocapstruct = ratception_struct;
analysisstruct_to_use = analysisstruct;

% do clustering:

analysisstruct_to_use.params.density_res = 1001; %, 1001; GC.density_res; %resolution of the map
analysisstruct_to_use.params.density_width = 3;%GC.density_width;% 1 default
analysisstruct_to_use.params.expansion_factor = GC.expansion_factor; %add a little room to the map after kernel smoothing
analysisstruct_to_use.params.density_threshold = GC.density_threshold; %remove regions in plots with low density
analysisstruct_to_use.matchedconds = {[unique(analysisstruct.condition_inds)]}; %if running over multiple conditions
analysisstruct_to_use.conditions_to_run = [unique(analysisstruct.condition_inds)];


params.reorder=1;
analysisstruct_to_use = compute_analysis_clusters_demo(analysisstruct_to_use,params); % check line 248, cluster_tsne_map.m
disp('%% Done clustering %%')

%%
%% behavior plots and movies
analysisstrict_to_use.conditionnames = ratname;
analysisstruct_to_use.ratnames = ratname;
analysisstruct_to_use.filesizes = {size(mocapstruct.aligned_mean_position,1 );};

%% plot a tsne map -- see plotting script for parameter definitions
h1=figure(609);
clf;
params.nameplot=1;
params.density_plot =0;
params.watershed = 1;
params.sorted = 1;
params.markersize = 2;
params.coarseboundary =0;
params.do_coarse = 0;
% plot tsne
plot_clustercolored_tsne(analysisstruct_to_use,1,params.watershed,h1,params)
set(h1,'Position',([100 100 1100 1100]))
export_folder = fullfile(GC.temp_root);

fig_name = ['High-dens_clustering'];
exportgraphics(gcf, fullfile(export_folder, [fig_name,'.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');

% bird specific axes
axisparams.zlim = ([200 300]);
axisparams.xlim = ([-400 400]);
axisparams.ylim = ([-400 400]);

%% make video to see each cluster

h=figure(370);
CL = cell(length(seq_c_idx),1);
for seq_ic = 1:numel(seq_cls)
    this_cls = seq_cls(seq_ic);    fprintf('ic = %i - ', this_cls)
     CL(seq_ic) =  {find(hierarchystruct.clustered_behavior{1}==this_cls)};
    if this_cls==0,  fprintf('\n'),continue, end
    animate_markers_nonaligned_fullmovie_demo(analysisstruct_to_use.mocapstruct_reduced_agg{1},...
        find(hierarchystruct.clustered_behavior{1}==this_cls), h, [], ['ic =  ',num2str(this_cls)]);

end

[cls, c_idx, r] = unique(analysisstruct_to_use.annot_reordered{end}, 'stable');

fig_poses = figure('pos', [10,300,1500,1900]);
nclus = numel(cls);

for ic = 1:numel(cls)
    % subplot(n_rows, n_cols, ic)
    this_cls = cls(ic);
    frames_to_plot = find(analysisstruct_to_use.annot_reordered{end}==this_cls);
    frames_to_plot = frames_to_plot(1:min(frames_to_plot(end), 1000)); % limit to 1000 frames
    fprintf('ic = %i - \n', this_cls)
    animate_markers_nonaligned_fullmovie_demo(analysisstruct_to_use.mocapstruct_reduced_agg{1},...
        frames_to_plot,fig_poses, [],['cl nr :  ', num2str(this_cls)]);
    title(this_cls)
end
