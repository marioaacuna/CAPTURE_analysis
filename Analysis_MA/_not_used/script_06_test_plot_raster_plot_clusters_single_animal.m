% Given data
animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
animal_ID = '328';
% Load analysis structure
analysisstruct_filename = fullfile("D:\test_CAPTURE\CAPTURE\analysisstruct_clusters.mat");
load(analysisstruct_filename)

%load predictions
rootpath = 'D:\test_CAPTURE';
filename_predictions = fullfile(rootpath, 'agg_predictions.mat');
load(filename_predictions, 'predictions', 'animal_frames_identifier')
% Set conds ids
long_animal_frames_identifier = repelem(animal_frames_identifier,3);
animal_list_used_after_analysis =  long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});
cond_inds = zeros(1,length(analysisstruct.condition_inds));
for iid = 1:length(animal_list)
    animal_id = animal_list{iid};
    
    idx = ismember(animal_list_used_after_analysis, animal_id);
    cond_inds(idx) = iid;

end


% Constants
n_frames_total = 216000;

% 1. Find index of animal '328'
animal_idx = find(strcmp(animal_list, animal_ID));

% 2. Get the condition indices for animal '328'
specific_cond_inds = find(ismember(cond_inds,animal_idx));

% 3. Extract specific frames and outcomes based on the condition indices
frames = analysisstruct.frames_with_good_tracking{1,1}(specific_cond_inds);
outcomes = analysisstruct.annot_reordered_matched{1, 1}(specific_cond_inds);
% outcomes = analysisstruct.annotation_vec{1, 2}(specific_cond_inds);


% Initialize the clusters array
clusters = zeros(1, n_frames_total);

for i = 1:length(frames)
    frame_val = ceil(frames(i) / 3);
    
    if i == 1
        relative_start = frame_val; % This will be adjusted below
        if animal_idx > 1 && frame_val > 50
            previous_frame = analysisstruct.frames_with_good_tracking{1,1}(specific_cond_inds(1) - 1);
            relative_start = frame_val - ceil(previous_frame / 3);
        end
        % clusters(relative_start) = outcomes(i);
        clusters(relative_start:frame_val) = outcomes(i);

        current_index = relative_start + 1;
    else
        gap = ceil((frames(i) - frames(i-1)) / 3);
        if gap > 17 % = ceil(50/3)
            clusters(current_index:current_index+gap-1) = 0; % fill gaps with zeros
            current_index = current_index + gap;
        else
            clusters(current_index:current_index+gap-1) = outcomes(i);
            current_index = current_index + gap-1;
        end
    end
end

% save clusters
save("D:\test_CAPTURE\328\cluster_vector_reordered.mat", "clusters")

%% Plot
% Given clusters

% Create the raster plot
figure_raster=figure;
hold on;
for i = 1:length(clusters)
    if clusters(i) ~= 0 % We can skip 0s or plot them, based on your preference
        plot(i, clusters(i), 'color','k', 'Marker','square', 'MarkerSize', 2); % change 'k.' to another color or marker if desired
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
% Given frame rate
frame_rate = 120; % fps

% Compute total number of frames for 30 minutes
total_frames = length(clusters);


% Set x-axis ticks for every 5 minutes
xticks(0:5*60*frame_rate:total_frames);
xticklabels(0:5:30);
xlabel('Time (minutes)');
%% save figure
export_fig('C:\Users\acuna\iCloudDrive\Grant_proposals\SNF_project_grant_2023\Figures\raster_cluster_328.pdf', figure_raster)

