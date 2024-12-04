% Initialization
clear;
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Load Data
load(GC.filename_analysis, 'analysisstruct')
load(GC.filename_ratception, 'ratception_struct');
load(GC.filename_predictions, 'animal_condition_identifier');
input_params.repfactor = GC.repfactor;

% Preprocess data
markers_aligned_ds = load_aligned_markers(ratception_struct.markers_aligned_preproc, input_params.repfactor, 15);

% Extract conditions
frame_identifiers = animal_condition_identifier;
conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);
unique_conditions = unique(conditions);

%% Analyse angles.
% Get marker data
markers = markers_aligned_ds;

% Joint positions (centered to SpineM)
knee = markers.KneeL;
ankle = markers.AnkleL;
paw = markers.HindpawL;

% 1. Calculate angles between markers in 3D space
leg_angles = calculate_leg_angles(knee, ankle, paw);

% 2. Calculate velocities of markers in 3D space
fps = 100; % Adjust this later (granularity / expansion_factor))
knee_vel = calculate_velocity(knee, fps);
ankle_vel = calculate_velocity(ankle, fps);
paw_vel = calculate_velocity(paw, fps);

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
            angle_name = {'Knee Flexion/Extension', 'Ankle Dorsiflexion/Plantarflexion'};
            angle_data = data(:, angle_idx);

            baseline_data = angle_data(ismember(conditions, 'B'), :);
            % Initialize cell array to store data for each condition
            condition_data = cell(length(unique_conditions), 1);

            for c = 1:length(unique_conditions)
                condition = unique_conditions{c};
                condition_mask = strcmp(conditions, condition);
                condition_data{c} = angle_data(condition_mask, :);% - mean(baseline_data);
            end

            % Perform statistical analysis
            all_data = [];
            group_labels = [];
            for c = 1:length(unique_conditions)
                all_data = [all_data; condition_data{c}];
                group_labels = [group_labels; repmat(unique_conditions(c), size(condition_data{c}, 1), 1)];
            end

            % Perform Kruskal-Wallis test
            [p_kw, tbl_kw, stats_kw] = kruskalwallis(all_data, group_labels, 'off');
            
            % Perform ANOVA test
            [p_anova, tbl_anova, stats_anova] = anova1(all_data, group_labels, 'off');

            % Display results
            fprintf('Kruskal-Wallis test for %s, p-value: %.4f\n', angle_name{angle_idx}, p_kw);
            fprintf('ANOVA test for %s, p-value: %.4f\n', angle_name{angle_idx}, p_anova);

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
            if p_kw < 0.05
                plot(1:length(unique_conditions), max(means + sems) * 1.1 * ones(1, length(unique_conditions)), 'k-');
                text(mean(1:length(unique_conditions)), max(means + sems) * 1.15, sprintf('KW p = %.3f', p_kw), ...
                    'HorizontalAlignment', 'center');
            end
            if p_anova < 0.05
                plot(1:length(unique_conditions), max(means + sems) * 1.2 * ones(1, length(unique_conditions)), 'k--');
                text(mean(1:length(unique_conditions)), max(means + sems) * 1.25, sprintf('ANOVA p = %.3f', p_anova), ...
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
            condition_data{c} = data(condition_mask, :); %- mean(baseline_data);
        end

        % Perform statistical analysis
        all_data = [];
        group_labels = [];
        for c = 1:length(unique_conditions)
            all_data = [all_data; condition_data{c}];
            group_labels = [group_labels; repmat(unique_conditions(c), size(condition_data{c}, 1), 1)];
        end

        % Perform Kruskal-Wallis test
        [p_kw, tbl_kw, stats_kw] = kruskalwallis(all_data(:, 1), group_labels, 'off'); % Example for first column
        
        % Perform ANOVA test
        [p_anova, tbl_anova, stats_anova] = anova1(all_data(:, 1), group_labels, 'off');

        % Display results
        fprintf('Kruskal-Wallis test for %s, p-value: %.4f\n', fieldname, p_kw);
        fprintf('ANOVA test for %s, p-value: %.4f\n', fieldname, p_anova);

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
        if p_kw < 0.05
            plot(1:length(unique_conditions), max(means + sems) * 1.1 * ones(1, length(unique_conditions)), 'k-');
            text(mean(1:length(unique_conditions)), max(means + sems) * 1.15, sprintf('KW p = %.3f', p_kw), ...
                'HorizontalAlignment', 'center');
        end
        if p_anova < 0.05
            plot(1:length(unique_conditions), max(means + sems) * 1.2 * ones(1, length(unique_conditions)), 'k--');
            text(mean(1:length(unique_conditions)), max(means + sems) * 1.25, sprintf('ANOVA p = %.3f', p_anova), ...
                'HorizontalAlignment', 'center');
        end
    end
end

disp('done')

