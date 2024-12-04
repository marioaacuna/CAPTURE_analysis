% Initialization
clear;
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Load Data
load(GC.filename_analysis, 'analysisstruct');

% Extract conditions
load(GC.filename_predictions, 'animal_condition_identifier');
input_params.repfactor = GC.repfactor;
upsampled_identifiers = repelem(animal_condition_identifier, input_params.repfactor);
good_frames = analysisstruct.frames_with_good_tracking{1, 1};
frame_identifiers = upsampled_identifiers(good_frames);
conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);
unique_conditions = unique(conditions);

%% Analyse angles.
% Get marker data
markers = analysisstruct.mocapstruct_reduced_agg{1, 1}.markers_aligned_preproc;

% Joint positions (centered to SpineM)
knee = markers.KneeL;
ankle = markers.AnkleL;
paw = markers.HindpawL;

% 1. Calculate angles between markers in 3D space
leg_angles = calculate_leg_angles(knee, ankle, paw);

% 2. Calculate velocities of markers in 3D space
%leg_velocities = calculate_velocity({knee, ankle, paw});

% 3. Joint velocities
fps = 8.3; % Adjust this later (granularity / expansion_factor))
knee_vel = calculate_velocity(knee, fps);
ankle_vel = calculate_velocity(ankle, fps);
paw_vel = calculate_velocity(paw, fps);


% TODO: continue with the separation of conditions and statistical analysis

% Separate data by conditions
fields = {'leg_angles', 'knee_vel', 'ankle_vel', 'paw_vel'};
data_struct = struct('leg_angles', leg_angles, ...
                     'knee_vel', knee_vel, 'ankle_vel', ankle_vel, 'paw_vel', paw_vel);

for f = 1:length(fields)
    fieldname = fields{f};
    data = data_struct.(fieldname);

    if strcmp(fieldname, 'leg_angles')
        % Analyze knee and ankle angles separately
        for angle_idx = 1:2
            angle_name = {'knee_angle', 'ankle_angle'};
            angle_data = data(:, angle_idx);

            baseline_data = angle_data(ismember(conditions, 'B'), :);
            % Initialize cell array to store data for each condition
            condition_data = cell(length(unique_conditions), 1);

            for c = 1:length(unique_conditions)
                condition = unique_conditions{c};
                condition_mask = strcmp(conditions, condition);
                condition_data{c} = angle_data(condition_mask, :) - mean(baseline_data);
            end

            % Perform statistical analysis
            all_data = [];
            group_labels = [];
            for c = 1:length(unique_conditions)
                all_data = [all_data; condition_data{c}];
                group_labels = [group_labels; repmat(unique_conditions(c), size(condition_data{c}, 1), 1)];
            end

            % Perform Kruskal-Wallis test
            [p, tbl, stats] = kruskalwallis(all_data, group_labels, 'off');

            % Display results
            fprintf('Kruskal-Wallis test for %s, p-value: %.4f\n', angle_name{angle_idx}, p);

            % Create figure for comparison
            figure('Position', [100 100 800 600], 'Color', 'w');

            % Calculate means and SEMs for each condition
            means = zeros(1, length(unique_conditions));
            sems = zeros(1, length(unique_conditions));
            for c = 1:length(unique_conditions)
                means(c) = median(condition_data{c});
                sems(c) = std(condition_data{c}) / sqrt(size(condition_data{c}, 1));
            end

            % Create bar plot with error bars
            b = bar(means, 'FaceColor', 'flat');
            hold on;
            errorbar(1:length(unique_conditions), means, sems, 'k', 'LineStyle', 'none', 'CapSize', 10);

            % Customize plot
            colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
                      [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
            for c = 1:length(unique_conditions)
                b.CData(c,:) = colors{c};
            end

            xlabel('Condition');
            ylabel(angle_name{angle_idx});
            title(sprintf('%s Across Conditions', angle_name{angle_idx}), 'Interpreter', 'none');
            set(gca, 'XTick', 1:length(unique_conditions), 'XTickLabel', unique_conditions);
            box off;
            set(gca, 'TickDir', 'out');

            % Add significance marker if test is significant
            if p < 0.05
                plot(1:length(unique_conditions), max(means + sems) * 1.1 * ones(1, length(unique_conditions)), 'k-');
                text(mean(1:length(unique_conditions)), max(means + sems) * 1.15, sprintf('p = %.3f', p), ...
                    'HorizontalAlignment', 'center');
            end
        end
    else
        baseline_data = data(ismember(conditions, 'B'), :);
        % Initialize cell array to store data for each condition
        condition_data = cell(length(unique_conditions), 1);

        for c = 1:length(unique_conditions)
            condition = unique_conditions{c};
            condition_mask = strcmp(conditions, condition);
            condition_data{c} = data(condition_mask, :) - mean(baseline_data);
        end

        % Perform statistical analysis
        all_data = [];
        group_labels = [];
        for c = 1:length(unique_conditions)
            all_data = [all_data; condition_data{c}];
            group_labels = [group_labels; repmat(unique_conditions(c), size(condition_data{c}, 1), 1)];
        end

        % Perform Kruskal-Wallis test
        [p, tbl, stats] = kruskalwallis(all_data(:, 1), group_labels, 'off'); % Example for first column

        % Display results
        fprintf('Kruskal-Wallis test for %s, p-value: %.4f\n', fieldname, p);

        % Create figure for comparison
        figure('Position', [100 100 800 600], 'Color', 'w');

        % Calculate means and SEMs for each condition
        means = zeros(1, length(unique_conditions));
        sems = zeros(1, length(unique_conditions));
        for c = 1:length(unique_conditions)
            means(c) = median(condition_data{c}(:, 1));  % Example for first column
            sems(c) = std(condition_data{c}(:, 1)) / sqrt(size(condition_data{c}, 1));
        end

        % Create bar plot with error bars
        b = bar(means, 'FaceColor', 'flat');
        hold on;
        errorbar(1:length(unique_conditions), means, sems, 'k', 'LineStyle', 'none', 'CapSize', 10);

        % Customize plot
        colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
                  [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
        for c = 1:length(unique_conditions)
            b.CData(c,:) = colors{c};
        end

        xlabel('Condition');
        ylabel(fieldname);
        title(sprintf('%s Across Conditions', fieldname), 'Interpreter', 'none');
        set(gca, 'XTick', 1:length(unique_conditions), 'XTickLabel', unique_conditions);
        box off;
        set(gca, 'TickDir', 'out');

        % Add significance marker if test is significant
        if p < 0.05
            plot(1:length(unique_conditions), max(means + sems) * 1.1 * ones(1, length(unique_conditions)), 'k-');
            text(mean(1:length(unique_conditions)), max(means + sems) * 1.15, sprintf('p = %.3f', p), ...
                'HorizontalAlignment', 'center');
        end
    end
end

disp('done')

function angles = calculate_leg_angles(knee, ankle, paw)
    % this function calculates the angles between the segments of the leg
    % knee, ankle and paw are the 3D coordinates of the markers
    % output is a matrix with the angles between the segments
    % where the first column is the knee angle and the second column is the ankle angle

    num_frames = size(knee, 1);
    angles = zeros(num_frames, 2); % [knee_angle, ankle_angle]
    
    for i = 1:num_frames
        % Vectors for segments
        thigh_vec = knee(i,:);  % From origin (SpineM) to knee
        shank_vec = ankle(i,:) - knee(i,:);
        foot_vec = paw(i,:) - ankle(i,:);
        
        % For knee angle
        % Get the cross product to determine rotation direction
        cross_knee = cross(thigh_vec, shank_vec);
        % Use sign of y-component (assuming sagittal plane primary motion)
        knee_direction = sign(cross_knee(2));
        % Calculate magnitude
        knee_mag = acosd(dot(-thigh_vec, shank_vec) / (norm(thigh_vec) * norm(shank_vec)));
        % Apply direction
        angles(i,1) = knee_direction * knee_mag;
        
        % For ankle angle
        cross_ankle = cross(shank_vec, foot_vec);
        ankle_direction = sign(cross_ankle(2));
        ankle_mag = acosd(dot(-shank_vec, foot_vec) / (norm(shank_vec) * norm(foot_vec)));
        angles(i,2) = ankle_direction * ankle_mag;
    end
end

function vel = calculate_velocity(position, fps)
    % Calculate velocity using central difference
    vel = diff(position) * fps;
    % Add duplicate of last velocity to match original length
    vel = [vel; vel(end,:)];
end