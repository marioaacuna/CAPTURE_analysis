%%
clusters_to_plot = [191,143, 100, 101,102,111,106,117,196,104,105,198,198,111,109,117,119,192,194,193,197,115,120,146,148,168,200];
cls = unique(clusters_to_plot);


    fig_poses = figure('pos', [10,300,1500,1900]);
    nclus = numel(cls);
    n_rows = ceil(sqrt(nclus));
    n_cols = ceil(sqrt(nclus));
    for ic = 1:nclus
        subplot(n_rows, n_cols, ic)
        this_cls = cls(ic);
        fprintf('ic = %i - \n', this_cls)
        plot_mean_cluster_aligned(analysisstruct.mocapstruct_reduced_agg{1},...
            find(analysisstruct.annot_reordered{end}==this_cls),['cl nr :  ', num2str(this_cls)]);
        title(this_cls)
    end


%%
 h=figure(370,'pos', [10,300,1500,1900]);
 seq_cls = 143;
% CL = cell(length(seq_c_idx),1);
for seq_ic = 1:numel(seq_cls)
    this_cls = seq_cls(seq_ic);    fprintf('ic = %i - ', this_cls)
     % CL(seq_ic) =  {find(hierarchystruct.clustered_behavior{1}==this_cls)};
    if this_cls==0,  fprintf('\n'),continue, end
    animate_markers_nonaligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},...
        find(analysisstruct.annot_reordered{end}==this_cls), h, [], ['ic =  ',num2str(this_cls)]);

end

%%
% in case some things are not loaded

filename_predictions = GC.filename_predictions;
    load(filename_predictions)
long_animal_frames_identifier = repelem(animal_condition_identifier,3);
animal_list_used_after_analysis =  long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});

%%
clc
 try
close(figure_cl_to_take_amim, figure_traces)

 end
 seq_cls = 289; % 171, 170, 174, 195
% CL = cell(length(seq_c_idx),1);
trace = analysisstruct.annot_reordered{end}==seq_cls;
figure_traces= figure('pos', [100,50, 1000,500]); hold on
plot(trace); title(num2str(seq_cls))
plot(endsWith(animal_list_used_after_analysis, 'F'), 'r')

hold off
 figure_cl_to_take_amim=figure('pos', [10,300,1000,850]);

for seq_ic = 1:numel(seq_cls)
    this_cls = seq_cls(seq_ic);    fprintf('ic = %i - ', this_cls)
     % CL(seq_ic) =  {find(hierarchystruct.clustered_behavior{1}==this_cls)};
    if this_cls==0,  fprintf('\n'),continue, end
    animate_markers_aligned_fullmovie_demo(analysisstruct.mocapstruct_reduced_agg{1},...
        find(analysisstruct.annot_reordered{end}==this_cls), figure_cl_to_take_amim, ['ic =  ',num2str(this_cls)]);

end
