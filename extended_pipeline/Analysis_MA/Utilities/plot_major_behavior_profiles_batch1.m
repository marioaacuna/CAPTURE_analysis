%% Behavioral Analysis Pipeline for MoCap Data
% This script takes concatenations of BSHFN data from all batches and
% performs detection of several coarse behaviors
% Additionally, it analysis genetal meobility features (distance travelled, 
% mean velocity, etc)

%% INIT
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
% Upsammple animal_condition_identifier
rep_factor = size(ratception_struct.aligned_mean_position,1) / length(animal_condition_identifier);
animal_condition_identifier = repelem(animal_condition_identifier, rep_factor);


unique_conditions = {'B', 'S', 'F', 'H', 'N'};
unique_animals = unique(cellfun(@(x) x(1:end-2), animal_condition_identifier, 'UniformOutput', false));

% Get unique animal/condition combinations
unique_animal_conditions = unique(animal_condition_identifier, 'stable');

%% 3. Process coarse behaviors all animals

% behavior_names = {'rearing', 'grooming', 'left_paw_licking', 'walking', 'quiet'};
behavior_names = {'rearing', 'grooming', 'walking', 'quiet'};


all_behaviors = struct();
condition_labels = {};


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

fprintf('=== ANALYSIS COARSE BEHAVIORS COMPLETE ===\n');


%% INIT movement/displacement analysis %%

fprintf('\n === Displacement analysis ===\n');

% Sampling frequency and conversion factors
sampling_freq = 300; % Hz
arena_diameter = 24; % cm: not used
center_diameter = 10; % cm
center_radius = center_diameter / 2; % cm
velocity_threshold = 1; % cm/s - threshold to consider mouse as "moving"

% Initialize storage for movement metrics
movement_metrics = struct();

% Get animal IDs and clean them for proper matching
animal_ids = fieldnames(all_behaviors);
clean_animal_ids = cellfun(@(x) x(4:end), animal_ids, 'UniformOutput', false); % Remove 'ID_' prefix

% 1. Calculate per-animal/condition movement metrics
fprintf('Calculating movement metrics for each animal...\n');

for i = 1:length(unique_animal_conditions)
    current_animal_condition = unique_animal_conditions{i};
    if strcmp(current_animal_condition, '1636_B'), continue, end
    % Find corresponding fieldname in all_behaviors
    matching_idx = ismember(clean_animal_ids, current_animal_condition);
    if any(matching_idx)
        matching_field = animal_ids{matching_idx};
        
        % Get frame indices for this animal
        frame_indices = find(strcmp(animal_condition_identifier, current_animal_condition));
        
        if ~isempty(frame_indices)
            % Extract SpineF position data (x, y, z coordinates in mm)
            spine_pos = ratception_struct.markers_preproc.SpineF(frame_indices, :);
            
            % Convert from mm to cm
            spine_pos_cm = spine_pos / 10;
            
            % 1. Calculate total distance moved (2D movement in x-y plane)
            if size(spine_pos_cm, 1) > 1
                displacement_2d = diff(spine_pos_cm(:, 1:2));
                frame_distances = sqrt(sum(displacement_2d.^2, 2));
                total_distance = sum(frame_distances); % in cm
                
                % 2. Calculate velocities
                dt = 1/sampling_freq; % time between frames in seconds
                velocities = frame_distances / dt; % cm/s
                
                % Mean velocity when moving (above threshold)
                moving_mask = velocities > velocity_threshold;
                if any(moving_mask)
                    mean_velocity_moving = mean(velocities(moving_mask));
                    time_spent_moving = sum(moving_mask) / sampling_freq; % seconds
                    percent_time_moving = (sum(moving_mask) / length(velocities)) * 100;
                else
                    mean_velocity_moving = 0;
                    time_spent_moving = 0;
                    percent_time_moving = 0;
                end
                
                % Overall mean velocity
                mean_velocity_overall = mean(velocities);
                
                % 3. Calculate time spent in center (draft implementation)
                % Find center of arena for this animal (mean position as approximation)
                % animal_center_x = mean(spine_pos_cm(:, 1));
                % animal_center_y = mean(spine_pos_cm(:, 2));

                % Find the boundaries of movement
                x_min = min(spine_pos_cm(:,1));
                x_max = max(spine_pos_cm(:,1));
                y_min = min(spine_pos_cm(:,2));
                y_max = max(spine_pos_cm(:,2));

                % Calculate the center of the arena
                animal_center_x = (x_min + x_max) / 2;
                animal_center_y = (y_min + y_max) / 2;

                % Calculate distances from animal's estimated arena center
                distances_from_center = sqrt((spine_pos_cm(:, 1) - animal_center_x).^2 + ...
                                           (spine_pos_cm(:, 2) - animal_center_y).^2);
                
                % Frames spent in center (within center_radius)
                in_center_mask = distances_from_center <= center_radius;
                time_in_center = sum(in_center_mask) / sampling_freq; % seconds
                percent_time_in_center = (sum(in_center_mask) / length(distances_from_center)) * 100;
                
                % 4. Time spent rearing (using existing behavioral data)
                rearing_data = all_behaviors.(matching_field).rearing;
                time_rearing = sum(rearing_data) / sampling_freq; % seconds
                percent_time_rearing = (sum(rearing_data) / length(rearing_data)) * 100;
                
                % 5. Movement complexity metrics
                % 5.1 Calculate path tortuosity (ratio of actual path to straight-line distance)
                if size(spine_pos_cm, 1) > 2
                    % Total path length (sum of all segments)
                    total_path_length = total_distance; % already calculated above
                    
                    % Straight-line distance from start to end
                    start_pos = spine_pos_cm(1, 1:2);
                    end_pos = spine_pos_cm(end, 1:2);
                    straight_line_distance = sqrt(sum((end_pos - start_pos).^2));
                    
                    % Tortuosity (1 = straight line, >1 = more tortuous)
                    if straight_line_distance > 0.1 % Avoid division by very small numbers
                        path_tortuosity = total_path_length / straight_line_distance;
                    else
                        path_tortuosity = 1; % If mouse didn't move much, consider it straight
                    end
                else
                    path_tortuosity = 1;
                end
                
                % 5.2 Compute angular velocity distribution
                if size(spine_pos_cm, 1) > 2
                    % Calculate heading angles between consecutive movements
                    movement_vectors = diff(spine_pos_cm(:, 1:2));
                    heading_angles = atan2(movement_vectors(:, 2), movement_vectors(:, 1));
                    
                    % Calculate angular changes (unwrap to handle 2π discontinuities)
                    unwrapped_angles = unwrap(heading_angles);
                    angular_changes = abs(diff(unwrapped_angles));
                    
                    % Angular velocity (rad/s)
                    angular_velocities = angular_changes / dt;
                    
                    % Statistics
                    mean_angular_velocity = mean(angular_velocities);
                    std_angular_velocity = std(angular_velocities);
                    max_angular_velocity = max(angular_velocities);
                    
                    % Count sharp turns (>90 degrees per frame)
                    sharp_turns = sum(angular_changes > pi/2);
                    turn_rate = sharp_turns / (length(angular_changes) / sampling_freq); % turns per second
                else
                    mean_angular_velocity = 0;
                    std_angular_velocity = 0;
                    max_angular_velocity = 0;
                    turn_rate = 0;
                end
                
                % 5.3 Analyze movement entropy/predictability
                if size(spine_pos_cm, 1) > 10
                    % Discretize movement directions into bins for entropy calculation
                    n_direction_bins = 8; % 8 cardinal/ordinal directions (N, NE, E, SE, S, SW, W, NW)
                    
                    % Calculate movement directions
                    movement_vectors = diff(spine_pos_cm(:, 1:2));
                    movement_angles = atan2(movement_vectors(:, 2), movement_vectors(:, 1));
                    
                    % Convert to degrees and normalize to 0-360
                    movement_angles_deg = mod(rad2deg(movement_angles), 360);
                    
                    % Bin the directions
                    bin_edges = linspace(0, 360, n_direction_bins + 1);
                    [direction_counts, ~] = histcounts(movement_angles_deg, bin_edges);
                    
                    % Calculate probabilities
                    direction_probs = direction_counts / sum(direction_counts);
                    direction_probs = direction_probs(direction_probs > 0); % Remove zero probabilities
                    
                    % Shannon entropy (higher = more unpredictable movement)
                    if ~isempty(direction_probs)
                        movement_entropy = -sum(direction_probs .* log2(direction_probs));
                        % Normalize by maximum possible entropy
                        max_entropy = log2(n_direction_bins);
                        normalized_entropy = movement_entropy / max_entropy;
                    else
                        movement_entropy = 0;
                        normalized_entropy = 0;
                    end
                    
                    % Movement predictability (inverse of normalized entropy)
                    movement_predictability = 1 - normalized_entropy;
                else
                    movement_entropy = 0;
                    normalized_entropy = 0;
                    movement_predictability = 1; % Very predictable if no movement
                end
                
            else
                % Handle case with insufficient data
                total_distance = 0;
                mean_velocity_moving = 0;
                mean_velocity_overall = 0;
                time_spent_moving = 0;
                percent_time_moving = 0;
                time_in_center = 0;
                percent_time_in_center = 0;
                time_rearing = 0;
                percent_time_rearing = 0;
                path_tortuosity = 1;
                mean_angular_velocity = 0;
                std_angular_velocity = 0;
                max_angular_velocity = 0;
                turn_rate = 0;
                movement_entropy = 0;
                normalized_entropy = 0;
                movement_predictability = 1;
            end
            
            % Store metrics
            movement_metrics.(['ID_', current_animal_condition]) = struct(...
                'total_distance_cm', total_distance, ...
                'mean_velocity_moving_cm_s', mean_velocity_moving, ...
                'mean_velocity_overall_cm_s', mean_velocity_overall, ...
                'time_spent_moving_s', time_spent_moving, ...
                'percent_time_moving', percent_time_moving, ...
                'time_in_center_s', time_in_center, ...
                'percent_time_in_center', percent_time_in_center, ...
                'time_rearing_s', time_rearing, ...
                'percent_time_rearing', percent_time_rearing, ...
                'path_tortuosity', path_tortuosity, ...
                'mean_angular_velocity_rad_s', mean_angular_velocity, ...
                'std_angular_velocity_rad_s', std_angular_velocity, ...
                'max_angular_velocity_rad_s', max_angular_velocity, ...
                'turn_rate_per_s', turn_rate, ...
                'movement_entropy', movement_entropy, ...
                'normalized_entropy', normalized_entropy, ...
                'movement_predictability', movement_predictability, ...
                'condition', current_animal_condition(end) ...
            );
            
            fprintf('Animal %s: Distance=%.1f cm, Vel=%.1f cm/s, Tortuous=%.2f, Entropy=%.2f, Turns=%.1f/s\n', ...
                current_animal_condition, total_distance, mean_velocity_moving, ...
                path_tortuosity, normalized_entropy, turn_rate);
        end
    end
end

%% 5. Movement Metrics Summary and Visualization

% Organize data by condition for statistical analysis
condition_movement_data = struct();
movement_metric_names = {'total_distance_cm', 'mean_velocity_moving_cm_s', ...
                         'percent_time_in_center', 'percent_time_rearing', ...
                         'path_tortuosity', ...
                         'movement_predictability'};


metric_titles = {'Total Distance (cm)', 'Mean Velocity When Moving (cm/s)', ...
                'Time in Center (%)', 'Time Rearing (%)', ...
                'Path Tortuosity',  'Movement Predictability'}; % 'Turn Rate (turns/s)', 'Movement Entropy (normalized)'



for i = 1:length(unique_conditions)
    cond = unique_conditions{i};
    condition_movement_data.(cond) = struct();
    
    % Initialize arrays for each metric
    for j = 1:length(movement_metric_names)
        condition_movement_data.(cond).(movement_metric_names{j}) = [];
    end
end

% Fill condition data
movement_animal_ids = fieldnames(movement_metrics);
for i = 1:length(movement_animal_ids)
    animal_data = movement_metrics.(movement_animal_ids{i});
    cond = animal_data.condition;
    
    for j = 1:length(movement_metric_names)
        metric_name = movement_metric_names{j};
        condition_movement_data.(cond).(metric_name) = [condition_movement_data.(cond).(metric_name), animal_data.(metric_name)];
    end
end

% Create summary figure
fig_movement = figure('Position', [100, 100, 1800, 1200], 'Color', 'white');
set(gca, 'Color', 'white');


colors = lines(length(unique_conditions));

for metric_idx = 1:length(movement_metric_names)
    subplot(3, 3, metric_idx);
    
    metric_name = movement_metric_names{metric_idx};
    
    % Prepare data for bar plot
    means = zeros(1, length(unique_conditions));
    errors = zeros(1, length(unique_conditions));
    
    for cond_idx = 1:length(unique_conditions)
        cond = unique_conditions{cond_idx};
        data = condition_movement_data.(cond).(metric_name);
        
        if ~isempty(data)
            means(cond_idx) = mean(data);
            errors(cond_idx) = std(data) / sqrt(length(data)); % SEM
        end
    end
    
    % Create bar plot with error bars
    bar_handles = bar(1:length(unique_conditions), means, 'FaceColor', 'flat');
    hold on;
    
    % Color bars by condition
    for cond_idx = 1:length(unique_conditions)
        bar_handles.CData(cond_idx, :) = colors(cond_idx, :);
    end
    
    errorbar(1:length(unique_conditions), means, errors, 'k.', 'LineWidth', 1.5);
    
    set(gca, 'XTick', 1:length(unique_conditions), 'XTickLabel', unique_conditions);
    ylabel(metric_titles{metric_idx}, 'FontSize', 12, 'FontName', 'Arial', 'Color', 'black');
    xlabel('Condition', 'FontSize', 12, 'FontName', 'Arial', 'Color', 'black');
    title(metric_titles{metric_idx}, 'FontSize', 14, 'FontName', 'Arial', 'Color', 'black');
    grid off;
    box off
    set(gca, 'GridColor', 'black', 'GridAlpha', 0.3);
    set(gca, 'FontSize', 10, 'FontName', 'Arial');
    set(gca, 'Color', 'white');
    set(gca, 'TickDir', 'out');
end

sgtitle('Movement & Complexity Metrics by Condition', 'FontSize', 18, 'FontName', 'Arial', 'Color', 'black');

% Export figure
fig_filename = fullfile(GC.figure_folder,'movement_metrics_mocap.pdf');
exportgraphics(fig_movement, fig_filename);

%% 6. Statistical Analysis for Movement Metrics
fprintf('\n=== MOVEMENT METRICS STATISTICAL ANALYSIS ===\n');

for metric_idx = 1:length(movement_metric_names)
    metric_name = movement_metric_names{metric_idx};
    
    % Prepare data for ANOVA
    group_data = [];
    group_labels = [];
    
    for cond_idx = 1:length(unique_conditions)
        cond = unique_conditions{cond_idx};
        data = condition_movement_data.(cond).(metric_name);
        
        if ~isempty(data)
            group_data = [group_data; data(:)];
            group_labels = [group_labels; repmat({cond}, length(data), 1)];
        end
    end
    
    if length(unique(group_labels)) > 1 && length(group_data) > 1
        [p_value, ~, stats] = anova1(group_data, group_labels, 'off');
        fprintf('%s: p = %.3f\n', metric_titles{metric_idx}, p_value);
    end
end

fprintf('\n=== MOVEMENT ANALYSIS COMPLETE ===\n');
fprintf('\nNOTE: Center calculation is based on mean position per animal.\n');
fprintf('This is a draft implementation that may need refinement based on your specific arena setup.\n');

%% 7. Summary of Movement Complexity Metrics
fprintf('\n=== MOVEMENT COMPLEXITY METRICS SUMMARY ===\n');
fprintf('Path Tortuosity: 1.0 = straight line, >1.0 = more tortuous/winding path\n');
fprintf('Angular Velocity: Measures rotational movement speed (rad/s)\n');
fprintf('Turn Rate: Number of sharp turns (>90°) per second\n');
fprintf('Movement Entropy: 0-1 scale, higher = more unpredictable movement patterns\n');
fprintf('Movement Predictability: 1-entropy, higher = more stereotypic behavior\n\n');

% Calculate and display overall condition summaries
fprintf('CONDITION SUMMARIES:\n');
for cond_idx = 1:length(unique_conditions)
    cond = unique_conditions{cond_idx};
    fprintf('\nCondition %s:\n', cond);
    
    % Calculate means for this condition
    total_dist = condition_movement_data.(cond).total_distance_cm;
    tortuosity = condition_movement_data.(cond).path_tortuosity;
    turn_rate = condition_movement_data.(cond).turn_rate_per_s;
    entropy = condition_movement_data.(cond).normalized_entropy;
    
    if ~isempty(total_dist)
        fprintf('  Distance: %.1f ± %.1f cm\n', mean(total_dist), std(total_dist));
        fprintf('  Tortuosity: %.2f ± %.2f\n', mean(tortuosity), std(tortuosity));
        fprintf('  Turn Rate: %.2f ± %.2f turns/s\n', mean(turn_rate), std(turn_rate));
        fprintf('  Movement Entropy: %.3f ± %.3f\n', mean(entropy), std(entropy));
        fprintf('  N animals: %d\n', length(total_dist));
    end
end

fprintf('\n=== COMPLETE BEHAVIORAL AND MOVEMENT ANALYSIS FINISHED ===\n');

