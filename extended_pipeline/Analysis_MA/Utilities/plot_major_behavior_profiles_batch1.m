%% Behavioral Analysis Pipeline for MoCap Data
% Extract and visualize temporal behavioral patterns
clc, clear, close all
global GC
GC = general_configs();

%% Set publication-quality figure defaults
set(groot, 'defaultFigureColor', 'white');
set(groot, 'defaultAxesColor', 'white');
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultTextColor', 'black');
set(groot, 'defaultAxesXColor', 'black');
set(groot, 'defaultAxesYColor', 'black');
set(groot, 'defaultAxesZColor', 'black');

%% 1. Load and prepare data
logger('Loading data', 'INFO');
project_path = fullfile(GC.project_path, "data/0_preprocessing_BSHFN_300hz/");
% load(GC.filename_analysis, 'analysisstruct');
load(fullfile(project_path, 'agg_predictions.mat'), 'predictions', 'animal_condition_identifier');
load(fullfile(project_path, 'ratception_prediction.mat'), 'ratception_struct');

% % Downsample from 300Hz to 100Hz if needed
% if size(ratception_struct.markers_aligned_preproc.SpineM, 1) ~= length(animal_condition_identifier)
%     % Downsample markers to match condition identifier
%     marker_fields = fieldnames(ratception_struct.markers_aligned_preproc);
%     for i = 1:length(marker_fields)
%         ratception_struct.markers_aligned_preproc.(marker_fields{i}) = ...
%             ratception_struct.markers_aligned_preproc.(marker_fields{i})(1:3:end, :);
% 
%         ratception_struct.markers_preproc.(marker_fields{i}) = ...
%             ratception_struct.markers_preproc.(marker_fields{i})(1:3:end, :);
%     end
% end

%% 3. Process all animals
unique_conditions = {'B', 'S', 'F', 'H', 'N'};
unique_animals = unique(cellfun(@(x) x(1:end-2), animal_condition_identifier, 'UniformOutput', false));

% behavior_names = {'rearing', 'grooming', 'left_paw_licking', 'walking', 'quiet'};
behavior_names = {'rearing', 'grooming', 'walking', 'quiet'};


all_behaviors = struct();
condition_labels = {};

% Get unique animal/condition combinations
unique_animal_conditions = unique(animal_condition_identifier, 'stable');

for i = 1:length(unique_animal_conditions)
    current_animal_condition = unique_animal_conditions{i};
    [behaviors, ~] = extract_behavioral_features_from_mocap(ratception_struct, current_animal_condition, animal_condition_identifier);
    
    % Store behaviors for each animal-condition
    all_behaviors.(['ID_',current_animal_condition]) = behaviors;
    condition_labels{i} = current_animal_condition(end);
end

% Define color mapping for behaviors (reordered to match behavior_names)
behavior_colors = [
    1 0.5 0;    % rearing (orange)
    1 0 1;      % grooming (magenta)
    1 0 0;      % left_paw_attention (red)
    0 1 1;      % walking (cyan)
    0.5 0 1;    % quiet (purple)
];
behavior_colors = behavior_colors(1:length(behavior_names),:);


%% 4. Bottom Heatmap - Behavioral Timeline by Condition and Animal

% Get unique conditions and animal-condition combinations
%unique_conditions = {'B', 'S', 'F', 'H', 'N'};
%unique_animal_conditions = unique(animal_condition_identifier);


% Time binning for 30 min at 100Hz (180,000 total frames)
expected_frames_30min = 30 * 60 * 100; % 180,000 frames
bin_size = 2;
n_time_bins = expected_frames_30min / bin_size; %

% Get animal IDs without 'ID_' prefix for matching
animal_ids = fieldnames(all_behaviors);
clean_animal_ids = cellfun(@(x) x(4:end), animal_ids, 'UniformOutput', false); % Remove 'ID_' prefix

% Organize data by condition
condition_data = struct();
for i = 1:length(unique_conditions)
    cond = unique_conditions{i};
    cond_animals = unique_animal_conditions(contains(unique_animal_conditions, ['_' cond]));
    condition_data.(cond) = cond_animals;
end

% Create heatmap matrix
total_animals = length(unique_animal_conditions);
heatmap_matrix = zeros(total_animals, n_time_bins);

% Fill matrix with dominant behaviors
row_idx = 1;
condition_boundaries = [];

fig_heatmap = figure('Position', [100, 100, 1200, 600], 'Color', 'white');
set(gca, 'Color', 'white');

for i = 1:length(unique_conditions)
    cond = unique_conditions{i};
    cond_animals = condition_data.(cond);
    
    for j = 1:length(cond_animals)
        animal_id = cond_animals{j};
        
        % Find corresponding fieldname in all_behaviors
        matching_field = animal_ids{ismember(clean_animal_ids, animal_id)};
        
        if ~isempty(matching_field)
            % Get behavioral data for this animal
            max_animal_frames = length(all_behaviors.(matching_field).rearing);
            behavior_matrix = zeros(length(behavior_names), max_animal_frames);
            
            for k = 1:length(behavior_names)
                behavior_data = all_behaviors.(matching_field).(behavior_names{k});
                behavior_matrix(k, 1:length(behavior_data)) = behavior_data;
            end
            
            % Bin the data and find dominant behavior per bin
            for bin = 1:min(n_time_bins, ceil(max_animal_frames / bin_size))
                start_frame = (bin-1) * bin_size + 1;
                end_frame = min(bin * bin_size, max_animal_frames);
                
                if start_frame <= max_animal_frames
                    % Sum behaviors in this bin
                    bin_sums = sum(behavior_matrix(:, start_frame:end_frame), 2);
                    [max_val, dominant_behavior] = max(bin_sums);
                    
                    % Only assign if there's actual behavior (not all zeros)
                    if max_val > 0
                        heatmap_matrix(row_idx, bin) = dominant_behavior;
                    end
                end
            end
        end
        row_idx = row_idx + 1;
    end
    
    % Store boundary for condition separation
    condition_boundaries = [condition_boundaries, row_idx - 1];
end

% Create the heatmap
imagesc(heatmap_matrix);
colormap(behavior_colors);

% Add condition labels and boundaries
ytick_positions = [];
ytick_labels = {};
for i = 1:length(unique_conditions)
    if i == 1
        start_row = 1;
    else
        start_row = condition_boundaries(i-1) + 1;
    end
    end_row = condition_boundaries(i);
    
    % Add horizontal lines to separate conditions
    if i > 1
        hold on;
        plot([0.5, n_time_bins + 0.5], [start_row - 0.5, start_row - 0.5], 'k-', 'LineWidth', 2);
    end
    
    % Set y-axis labels
    mid_position = (start_row + end_row) / 2;
    ytick_positions = [ytick_positions, mid_position];
    ytick_labels{end+1} = unique_conditions{i};
end

% Format axes
set(gca, 'YTick', ytick_positions, 'YTickLabel', ytick_labels, 'FontSize', 12, 'FontName', 'Arial');
set(gca, 'XTick', []); % Remove x-axis ticks
xlabel('Time (s)', 'FontSize', 14, 'FontName', 'Arial', 'Color', 'black');
ylabel('Conditions', 'FontSize', 14, 'FontName', 'Arial', 'Color', 'black');
title('Behavioral Timeline by Condition (30 min @ 100Hz)', 'FontSize', 16, 'FontName', 'Arial', 'Color', 'black');

% Add time scale bar (200s as in reference image)
time_scale_length = 200; % seconds
scale_bins = time_scale_length / (bin_size / 100); % convert to bins (20 bins for 200s)
x_pos = n_time_bins - scale_bins - 50;
y_pos = total_animals + 2;

hold on;
plot([x_pos, x_pos + scale_bins], [y_pos, y_pos], 'k-', 'LineWidth', 6);
text(x_pos + scale_bins/2, y_pos + 1, '200 s', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial', 'Color', 'black');

% Add colorbar legend
cb = colorbar;
cb.Ticks = 1:length(behavior_names);
cb.TickLabels = behavior_names;
cb.Label.String = 'Behaviors';
cb.Label.FontSize = 14;
cb.Label.FontName = 'Arial';
cb.Label.Color = 'black';
cb.FontSize = 12;
cb.FontName = 'Arial';

% Adjust layout
xlim([0.5, n_time_bins + 0.5]);
ylim([0.5, total_animals + 0.5]);

% export figure
fig_filename = fullfile(GC.figure_folder,'heatmap_mocap.pdf');
exportgraphics(fig_heatmap, fig_filename)

%% 5. Behavioral Proportion Analysis
fig_proportions = figure('Position', [100, 700, 1200, 400], 'Color', 'white');
set(gca, 'Color', 'white');

% Calculate proportions per condition using proper animal/condition matching
condition_proportions = zeros(length(unique_conditions), length(behavior_names));
condition_errors = zeros(length(unique_conditions), length(behavior_names));

% Get animal IDs from all_behaviors fieldnames
animal_ids = fieldnames(all_behaviors);
clean_animal_ids = cellfun(@(x) x(4:end), animal_ids, 'UniformOutput', false); % Remove 'ID_' prefix

for i = 1:length(unique_conditions)
    cond = unique_conditions{i};
    % Find animals for this condition from unique_animal_conditions
    cond_animals = unique_animal_conditions(contains(unique_animal_conditions, ['_' cond]));
    
    if ~isempty(cond_animals)
        cond_props = [];
        for j = 1:length(cond_animals)
            animal_condition = cond_animals{j};
            
            % Find corresponding fieldname in all_behaviors
            matching_idx = ismember(clean_animal_ids, animal_condition);
            if any(matching_idx)
                matching_field = animal_ids{matching_idx};
                
                animal_props = zeros(1, length(behavior_names));
                for k = 1:length(behavior_names)
                    behavior_data = all_behaviors.(matching_field).(behavior_names{k});
                    animal_props(k) = mean(behavior_data);  % Use mean instead of sum/length
                end
                
                % % Debug: Print animal proportions for first few animals
                % if j <= 3
                %     fprintf('Animal %s (%s): [%.3f %.3f %.3f %.3f %.3f]\n', ...
                %         animal_condition, cond, animal_props);
                % end
                cond_props = [cond_props; animal_props];
            end
        end
        
        if ~isempty(cond_props)
            condition_proportions(i, :) = mean(cond_props, 1);
            condition_errors(i, :) = std(cond_props, [], 1) / sqrt(size(cond_props, 1));
        end
    end
end

% Plot proportions
bar_positions = 1:length(behavior_names);
bar_width = 0.15;

colors = lines(length(unique_conditions));
hold on;

for i = 1:length(unique_conditions)
    pos_offset = (i - (length(unique_conditions)+1)/2) * bar_width;
    errorbar(bar_positions + pos_offset, condition_proportions(i, :), ...
             condition_errors(i, :), 'o', 'LineWidth', 2, 'MarkerSize', 8, ...
             'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
end

set(gca, 'XTick', bar_positions, 'XTickLabel', behavior_names, 'FontSize', 12, 'FontName', 'Arial');
set(gca, 'FontSize', 12, 'FontName', 'Arial');
ylabel('Proportion of Time', 'FontSize', 14, 'FontName', 'Arial', 'Color', 'black');
xlabel('Behaviors', 'FontSize', 14, 'FontName', 'Arial', 'Color', 'black');
legend(unique_conditions, 'Location', 'best', 'FontSize', 12, 'FontName', 'Arial', 'TextColor', 'black');
title('Behavioral Proportions by Condition', 'FontSize', 16, 'FontName', 'Arial', 'Color', 'black');
grid on;
set(gca, 'GridColor', 'black', 'GridAlpha', 0.3);

% export figure
fig_filename = fullfile(GC.figure_folder,'proportions_mocap.pdf');
exportgraphics(fig_proportions, fig_filename)


%% 6. 3D Scatter Plot (One Animal Example per Condition)
downsample_factor = 10; 
fig_examples = figure('Position', [200, 200, 1500, 1000], 'Color', 'white');

% Get animal IDs and clean them for proper matching
animal_ids = fieldnames(all_behaviors);
clean_animal_ids = cellfun(@(x) x(4:end), animal_ids, 'UniformOutput', false); % Remove 'ID_' prefix

% Select first animal for each condition
% example_animals = {'ID_1641_B','ID_1389_S', 'ID_1639_F', 'ID_1383_H', 'ID_1636_N'};
example_animals = {};

for i = 1:length(unique_conditions)
    cond = unique_conditions{i};
    % Find animals for this condition from unique_animal_conditions
    cond_animals = unique_animal_conditions(contains(unique_animal_conditions, ['_' cond]));

    if ~isempty(cond_animals)
        % Take the first animal for this condition
        animal_to_take = cond_animals{5};
        % Find corresponding fieldname in all_behaviors
        matching_idx = ismember(clean_animal_ids, animal_to_take);
        if any(matching_idx)
            example_animals{i} = animal_ids{matching_idx};
        end
    end
end

% Create subplots for each condition
n_conditions = length(unique_conditions);
n_cols = ceil(sqrt(n_conditions));
n_rows = ceil(n_conditions / n_cols);

for i = 1:length(example_animals)
    if ~isempty(example_animals{i})
        subplot(n_rows, n_cols, i);
        
        example_animal = example_animals{i};
        clean_example_id = example_animal(4:end); % Remove 'ID_' prefix
        
        % Find frames for this specific animal/condition
        example_frame_indices = find(strcmp(animal_condition_identifier, clean_example_id));
        spine_pos = ratception_struct.markers_preproc.SpineF(example_frame_indices, :); % Use NON-ALIGNED data for trajectory
        
        % Create behavior color coding using the new behavior_colors from heatmap
        behavior_colors_3d = zeros(size(spine_pos, 1), 3);
        for j = 1:length(behavior_names)
            behavior_mask = all_behaviors.(example_animal).(behavior_names{j});
            if length(behavior_mask) == size(spine_pos, 1)
                behavior_colors_3d(behavior_mask, :) = repmat(behavior_colors(j, :), sum(behavior_mask), 1);
            end
        end
        
        % Handle frames with no dominant behavior (gray)
        no_behavior_mask = sum(behavior_colors_3d, 2) == 0;
        behavior_colors_3d(no_behavior_mask, :) = 0.3;
        
        scatter3(spine_pos(1:downsample_factor:end,1), spine_pos(1:downsample_factor:end,2), spine_pos(1:downsample_factor:end,3), 2, behavior_colors_3d(1:downsample_factor:end,:), 'filled', 'MarkerEdgeColor', 'none');
        xlabel('X Position (mm)', 'FontSize', 12, 'FontName', 'Arial', 'Color', 'black');
        ylabel('Y Position (mm)', 'FontSize', 12, 'FontName', 'Arial', 'Color', 'black');
        zlabel('Z Position (mm)', 'FontSize', 12, 'FontName', 'Arial', 'Color', 'black');
        title(sprintf('Condition %s - Animal %s', unique_conditions{i}, clean_example_id), 'FontSize', 14, 'FontName', 'Arial', 'Color', 'black');
        grid on;
        set(gca, 'GridColor', 'black', 'GridAlpha', 0.3);
        set(gca, 'FontSize', 10, 'FontName', 'Arial');
        set(gca, 'Color', 'white');
        view(45, 45);
        
        % % Set consistent axis limits for all subplots
        % if i == 1
        %     % Store limits from first plot
        %     x_lims = xlim;
        %     y_lims = ylim;
        %     z_lims = zlim;
        % else
        %     % Apply consistent limits
        %     xlim(x_lims);
        %     ylim(y_lims);
        %     zlim(z_lims);
        % end
    end
end

% Add overall title
sgtitle('Behavioral Trajectories - One Example per Condition', 'FontSize', 18, 'FontName', 'Arial', 'Color', 'black');

% Add legend for behaviors (outside the subplots)
legend_fig = figure('Position', [1400, 200, 200, 400], 'Color', 'white');
set(gca, 'Color', 'white');
legend_handles = [];
for i = 1:length(behavior_names)
    legend_handles(i) = scatter(NaN, NaN, 50, behavior_colors(i, :), 'filled');
    hold on;
end
legend(legend_handles, behavior_names, 'Location', 'best', 'FontSize', 12, 'FontName', 'Arial', 'TextColor', 'black');
axis off;
title('Behavior Legend', 'FontSize', 14, 'FontName', 'Arial', 'Color', 'black');


% export figure
fig_filename = fullfile(GC.figure_folder,'example_mocap.pdf');
exportgraphics(fig_examples, fig_filename)

%% 7. Statistical Analysis
fprintf('\n=== BEHAVIORAL ANALYSIS SUMMARY ===\n');

% Get animal IDs and clean them for proper matching
animal_ids = fieldnames(all_behaviors);
clean_animal_ids = cellfun(@(x) x(4:end), animal_ids, 'UniformOutput', false); % Remove 'ID_' prefix

% ANOVA for each behavior across conditions
for i = 1:length(behavior_names)
    behavior_name = behavior_names{i};
    
    % Prepare data for ANOVA
    group_data = [];
    group_labels = [];
    
    for j = 1:length(unique_conditions)
        cond = unique_conditions{j};
        % Find animals for this condition from unique_animal_conditions
        cond_animals = unique_animal_conditions(contains(unique_animal_conditions, ['_' cond]));
        
        for k = 1:length(cond_animals)
            animal_condition = cond_animals{k};
            
            % Find corresponding fieldname in all_behaviors
            matching_idx = ismember(clean_animal_ids, animal_condition);
            if any(matching_idx)
                matching_field = animal_ids{matching_idx};
                animal_prop = mean(all_behaviors.(matching_field).(behavior_name));
                group_data = [group_data; animal_prop];
                group_labels = [group_labels; {cond}];
            end
        end
    end
    
    if length(unique(group_labels)) > 1 && length(group_data) > 1
        [p_value, ~, stats] = anova1(group_data, group_labels, 'off');
        fprintf('%s:  p = %.3f \n', behavior_name, ...
                 p_value);
    end
end

fprintf('=== ANALYSIS COMPLETE ===\n');