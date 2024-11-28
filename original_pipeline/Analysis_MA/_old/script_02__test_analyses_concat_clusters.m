% load agg_analysisstrct
%
%%
cluster_assignments = agg_analysisstruct.annot_reordered{end};
animal_identifiers = animal_frames_identifier; % Assuming it's in the format described

% animal_to_condition = containers.Map(animal_list, conditions);
unique_clusters = unique(cluster_assignments);
n = length(unique_clusters);
count_formalin = zeros(n, length(animal_list));
count_saline = zeros(n, length(animal_list));


for i = 1:length(cluster_assignments)
    cluster_index = find(unique_clusters == cluster_assignments(i));
    animal_index = find(strcmp(animal_list, animal_identifiers(i)));

    current_condition = getCondition(animal_identifiers{i}, animal_list, conditions);
    
    if strcmp(current_condition, 'F')
        count_formalin(cluster_index, animal_index) = count_formalin(cluster_index, animal_index) + 1;
    elseif strcmp(current_condition, 'S')
        count_saline(cluster_index, animal_index) = count_saline(cluster_index, animal_index) + 1;
    end
end



% 
% 
% for i = 1:length(cluster_assignments)
%     cluster_index = find(unique_clusters == cluster_assignments(i));
%     animal_index = find(strcmp(animal_list, animal_identifiers(i)));
% 
%     if strcmp(animal_to_condition(animal_identifiers(i)), 'F')
%         count_formalin(cluster_index, animal_index) = count_formalin(cluster_index, animal_index) + 1;
%     elseif strcmp(animal_to_condition(animal_identifiers(i)), 'S')
%         count_saline(cluster_index, animal_index) = count_saline(cluster_index, animal_index) + 1;
%     end
% end
% 

sum_formalin = sum(count_formalin, 2);
sum_saline = sum(count_saline, 2);


num_formalin = sum(strcmp(conditions, 'F'));
num_saline = sum(strcmp(conditions, 'S'));

sum_formalin = sum_formalin / num_formalin;
sum_saline = sum_saline / num_saline;


difference = sum_formalin - sum_saline;

figure
bar(difference)
title('Difference in Cluster Appearances: Formalin vs. Saline')
xlabel('Cluster')
ylabel('Difference in Appearances')


%% do stats per animal
% Assuming the following data structures are available:
% - agg_analysisstruct.annot_reordered{end}: index of each cluster for all points
% - animal_list: list of animal IDs
% - conditions: list of conditions corresponding to animal_list
% - animal_frames_identifier: cell array of size (nframes x 1) indicating the animal ID for each frame

cluster_indices = agg_analysisstruct.annot_reordered{end};
unique_clusters = unique(cluster_indices);

% Initialize storage for cluster counts
cluster_counts = zeros(length(unique_clusters), length(animal_list));

% 1. Organize the data
% Count the appearances of each cluster for each animal
for i = 1:length(animal_list)
    animal_id = animal_list{i};
    frames_for_this_animal = strcmp(animal_frames_identifier, animal_id);
    clusters_for_this_animal = cluster_indices(frames_for_this_animal);
    for c = 1:length(unique_clusters)
        cluster_counts(c, i) = sum(clusters_for_this_animal == unique_clusters(c));
    end
end

% 2. Group by condition
saline_indices = strcmp(conditions, 'S');
formalin_indices = strcmp(conditions, 'F');

saline_counts = cluster_counts(:, saline_indices);
formalin_counts = cluster_counts(:, formalin_indices);

mean_saline = mean(saline_counts, 2);
mean_formalin = mean(formalin_counts, 2);
std_saline = std(saline_counts, 0, 2);
std_formalin = std(formalin_counts, 0, 2);

% 3. Statistical testing
p_values = zeros(length(unique_clusters), 1);
for c = 1:length(unique_clusters)
    [~, p_values(c)] = ttest2(saline_counts(c, :), formalin_counts(c, :));

     % p_values(c) = ranksum(saline_counts(c, :), formalin_counts(c, :));
end

% You can now interpret the p_values to see which clusters have significantly different appearances between conditions.
significant_clusters = find(p_values < 0.05);  % Assuming alpha = 0.05








function condition = getCondition(animal_id, animal_list, conditions)
    idx = find(strcmp(animal_list, animal_id));
    if isempty(idx)
        condition = 'Unknown'; % Or any default value
    else
        condition = conditions{idx};
    end
end
