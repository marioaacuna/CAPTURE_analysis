% Script to run clusterting baesd on concatenated analysis struct
clear
close all
clc

global GC
GC = general_configs();
% load analysisstruct
load('/home/mario/Documents/Projects/CAPTURE/data/0_preprocessing/concatenated_analysis_results.mat') % TODO: change it to GC
%% Load zvals
disp('%% INIT clustering %%')
% for now load zvalues alraedy saved. (so far not saved in analysis struct)
load('/home/mario/Documents/Projects/CAPTURE/data/0_preprocessing/zvals.mat')
%% Run clustering
% rename to keep consistency
analysisstruct = concatenated_analysisstruct;
clear concatenated_analysisstruct
%%

% for inter animal analysis:
animal_list = unique(frame_to_animal_condition, 'stable');
cond_inds = zeros(1,length(frame_to_animal_condition)); % sorting per animal
for iid = 1:length(animal_list)
    animal_ID = animal_list{iid};
    idx = ismember(frame_to_animal_condition, animal_ID);
    cond_inds(idx) = iid;
end
%%

analysisstruct.zValues = zvals_extra_all;
analysisstruct.params.density_res = GC.density_res; %resolution of the map
analysisstruct.params.density_width = 0.75;%GC.density_width;% 1 default
analysisstruct.params.expansion_factor = GC.expansion_factor; %add a little room to the map after kernel smoothing
analysisstruct.params.density_threshold = GC.density_threshold; %remove regions in plots with low density

analysisstruct.condition_inds = cond_inds;
analysisstruct.matchedconds = {unique(condition_indices)}; %if running over multiple conditions
analysisstruct.conditions_to_run = unique(analysisstruct.condition_inds );

params.reorder=1;
%% CLustering
analysisstruct = compute_analysis_clusters_demo(analysisstruct,params); % check line 248, cluster_tsne_map.m
disp('%% Done clustering %%')

%%
%% behavior plots and movies
ratname ='myrat';% 'test_mouse';

analysisstruct.conditionnames = ratname;
analysisstruct.ratnames = ratname;
% analysisstruct.filesizes = {size(mocapstruct.aligned_mean_position,1 );};

%% Save analysisstruct
analysis_filename = GC.filename_analysis;
save(analysis_filename, 'analysisstruct' , '-v7.3')

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

cluster_figure_filename = fullfile(GC.figure_folder, 'Tsne_clusters.pdf');
% export_fig(cluster_figure_filename, '-pdf', h1)
%% save
exportgraphics(h1, '/mnt/VMs/share/new_clusters_only_extras.pdf')
%%
%% Plot cluster poses
% condition to plot
idxcond = find(endsWith(analysisstruct.concatenation_info.animal_identifiers, 'S'));
idxcond = idxcond(1); % take the Nth animal
[cls, c_idx, r] = unique(analysisstruct.annot_reordered{idxcond}, 'stable');

plot_poses = 1;
if plot_poses
    % h= figure(370);
    % clf;

    figure('pos', [10,100,1500,1200]);
    nclus = numel(cls);
    n_rows = ceil(sqrt(nclus));
    n_cols = ceil(sqrt(nclus));
    for ic = 1:numel(cls)
        subplot(n_rows, n_cols, ic)
        this_cls = cls(ic);
        fprintf('ic = %i - \n', this_cls)
        plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1},...
            find(analysisstruct.annot_reordered{idxcond}==this_cls),['cl nr :  ', num2str(this_cls)]);
        title(this_cls)
    end
end


zvals_all = tsne([concatenated_analysisstruct.jt_features(1:1:end,:), concatenated_analysisstruct.extra_jt_features(1:1:end,end-5:end)], ...
'Perplexity', 500, ...
'Exaggeration', 4, ...
'Verbose', 1);