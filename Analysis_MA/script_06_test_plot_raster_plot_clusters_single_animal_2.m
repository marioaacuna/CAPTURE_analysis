clc, clear

% load analysisstruct
analysisstruct_filename = 'D:\test_CAPTURE\CAPTURE\analysisstruct_clusters.mat';
load(analysisstruct_filename)
agg_predictions_filename = 'D:\test_CAPTURE\agg_predictions.mat';
load(agg_predictions_filename)
% set variables:
animal_ID = '328';
% Given frame rate
ori_frame_rate = 120; % fps
upsamplig_factor = ceil(analysisstruct.mocapstruct_reduced_agg{1, 1}.fps/ori_frame_rate); 
% check main code for that ceil(analysisstruct.mocapstruct_reduced_agg{1, 1}.fps/120) : 120 : real FPS
% Constants
n_frames_total = 216000; % total number of frames recorded (120hz*30min)

long_animal_frames_identifier = repelem(animal_frames_identifier,upsamplig_factor);
animal_list_used_after_analysis =  long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});
% analysisstruct.condition_inds = [];
% for inter animal analysis:

animal_list = unique(animal_frames_identifier, 'stable');

cond_inds = zeros(1,length(analysisstruct.condition_inds));
for iid = 1:length(animal_list)
    animal_id= animal_list{iid};
    
    idx = ismember(animal_list_used_after_analysis, animal_id);
    cond_inds(idx) = iid;

end


%% Run
% Given data
% animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};

% 1. Find index of animal '328'
animal_idx = find(strcmp(animal_list, animal_ID));

% 2. Get the condition indices for animal '328'
specific_cond_inds = find(ismember(cond_inds,animal_idx));

% 3. Extract specific frames and outcomes based on the condition indices
% analysisstruct.frames_with_good_tracking is at least diff 50 frames

frames = analysisstruct.frames_with_good_tracking{1,1}(specific_cond_inds);
% frames_this_idx = floor(frames(1)/3) : 1: ceil(frames(end)/3) ;
frames_this_idx = floor(frames(1))  : 1: ceil(frames(end)) ;

outcomes = analysisstruct.annot_reordered_matched{1, 1}(specific_cond_inds); % sampled every 50
outcomes_rep = repelem(outcomes, 50);
% Initialize the clusters array
clusters = zeros(1, n_frames_total * upsamplig_factor);

for i = 1:length(frames)
    % frame_val = ceil(frames_this_idx(i) / upsamplig_factor);

    frame_val = frames(i);

    if i == 1
        relative_start = frame_val; % This will be adjusted below
        if animal_idx > 1 && frame_val > analysisstruct.tsnegranularity  % 50
            previous_frame = analysisstruct.frames_with_good_tracking{1,1}(specific_cond_inds(1) - 1); % diff of 50 frames 
            % relative_start = frame_val - ceil(previous_frame / 3);
            relative_start = frame_val - ceil(previous_frame);
            if relative_start == 50 % means it's only one frame
                relative_start = 1;
            elseif relative_start > 50 % meaning at the beginning there are no clusters
                % relative_start = previous_frame + 50; TODO, still I don't
                % know how to go on with this
                keyboard %#ok<KEYBOARDFUN>
            end
        end
        clusters(1:relative_start) = outcomes(i);
        % clusters(relative_start:frame_val) = outcomes(i);

        current_index = relative_start + 1;
    else

        % set the gap: cehck if the previous frame has a gap of 50
        gap = ceil((frames(i) - frames(i-1)));




        if gap > ceil(analysisstruct.tsnegranularity) 
             % keyboard % TODO
            clusters(current_index:current_index+gap-1) = 0; % fill gaps with zeros
            current_index = current_index + gap;
        else
            clusters(current_index:current_index+gap-1) = outcomes(i);
            % current_index = current_index + gap-1;
            current_index = current_index + gap;

        end
    end
end
% Downsample the clusters
clusters = downsample(clusters, upsamplig_factor);

save(fullfile("D:\test_CAPTURE",animal_ID, 'cluster_vector.mat' ), "clusters");
%% Plot
% Create the raster plot
figure_raster=figure('pos', [312,564,1514,420], 'color', 'k');
hold on;
for i = 1:length(clusters)
    if clusters(i) ~= 0 % We can skip 0s or plot them, based on your preference
        plot(i, clusters(i), 'color','w', 'Marker','square', 'MarkerSize', 2); % change 'k.' to another color or marker if desired
    end
end
hold off;

% Set title and axis labels
title('Raster Plot of Clusters');
xlabel('Frames');
ylabel('Cluster ID');

% You might want to set y-axis limits or ticks for clarity
ylim([0 max(clusters)+1]); % Adjust as necessary
yticks(0:10:max(clusters)); % Adjust the interval as necessary
% You might want to set y-axis limits or ticks for clarity
% xlim([0 length(clusters)+1]); % Adjust as necessary
% yticks(0:10:max(clusters)); % Adjust the interval as necessary

% % You might also want to adjust the x-axis ticks for clarity
% xticks(0:10000:length(clusters)); % example interval, adjust as necessary
box off;
set(gca, 'TickDir', 'out');
ax = gca;
ax.Color = [0 0 0 ];
ax.XColor = 'white';
ax.YColor = 'white';
% Compute total number of frames for 30 minutes
total_frames = length(clusters);


% Set x-axis ticks for every 5 minutes
xticks(0:5*60*ori_frame_rate:total_frames);
xticklabels(0:5:30);
xlabel('Time (minutes)');
%% save figure
% export_fig('C:\Users\acuna\iCloudDrive\Grant_proposals\SNF_project_grant_2023\Figures\raster_cluster_328_2.pdf', figure_raster)
export_fig(['raster_cluster_', animal_ID,' .pdf'], figure_raster)

