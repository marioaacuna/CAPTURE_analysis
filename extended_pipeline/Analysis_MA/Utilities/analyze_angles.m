
function analyze_angles(markers, walking_bouts, conditions, unique_conditions, frame_identifiers, animal_condition_identifier, do_export_figure)
    global GC
    % Joint positions (centered to SpineM)
    spineM = markers.SpineM;
    spineF = markers.SpineF;

    % Left leg
    kneeL = markers.KneeL;
    ankleL = markers.AnkleL;
    pawL = markers.HindpawL;

    % Right leg
    kneeR = markers.KneeR;
    ankleR = markers.AnkleR;
    pawR = markers.HindpawR;

    % Calculate angles between markers in 3D space
    l_leg_angles = calculate_leg_angles(spineM, spineF, kneeL, ankleL, pawL); % knee, ankle, hip
    r_leg_angles = calculate_leg_angles(spineM, spineF, kneeR, ankleR, pawR); % knee, ankle, hip

    % Sample the walking bouts
    l_leg_angles_at_walking = l_leg_angles(walking_bouts,:);
    r_leg_angles_at_walking = r_leg_angles(walking_bouts,:);

    % Collect data for all angles and legs
    angle_names = {'Knee', 'Ankle', 'Hip'};
    data_struct = struct('left_leg', l_leg_angles_at_walking, ...
                         'right_leg', r_leg_angles_at_walking);

    fields = fieldnames(data_struct);
    n_angles = size(l_leg_angles,2);

    all_data = struct();
    for angle_idx = 1:n_angles
        angle_name = angle_names{angle_idx};
        all_data.(angle_name) = struct();
        for f = 1:length(fields)
            fieldname = fields{f};
            data = data_struct.(fieldname);
            angle_data = data(:, angle_idx);

            % Initialize cell array to store data for each condition
            condition_data = cell(length(unique_conditions), 1);
            for c = 1:length(unique_conditions)
                condition = unique_conditions{c};
                condition_mask = strcmp(conditions(walking_bouts), condition);
                condition_data{c} = angle_data(condition_mask, :);
            end
            all_data.(angle_name).(fieldname) = condition_data;
        end
    end

    % Plotting
    for angle_idx = 1:n_angles
        angle_name = angle_names{angle_idx};
        condition_data_left = all_data.(angle_name).left_leg;
        condition_data_right = all_data.(angle_name).right_leg;

        % Calculate means and SEMs for each condition
        means_left = zeros(1, length(unique_conditions));
        sems_left = zeros(1, length(unique_conditions));
        means_right = zeros(1, length(unique_conditions));
        sems_right = zeros(1, length(unique_conditions));
        for c = 1:length(unique_conditions)
            means_left(c) = median(condition_data_left{c});
            sems_left(c) = std(condition_data_left{c}) / sqrt(size(condition_data_left{c}, 1));
            means_right(c) = median(condition_data_right{c});
            sems_right(c) = std(condition_data_right{c}) / sqrt(size(condition_data_right{c}, 1));
        end

        % Create figure for comparison
        figure('Position', [100 100 800 600], 'Color', 'w');

        % Create bar plot with error bars
        b = bar([means_left; means_right]', 'grouped');
        hold on;
        errorbar((1:length(unique_conditions)) - 0.15, means_left, sems_left, 'k', 'LineStyle', 'none', 'CapSize', 10);
        errorbar((1:length(unique_conditions)) + 0.15, means_right, sems_right, 'k', 'LineStyle', 'none', 'CapSize', 10);

        % Customize plot
        colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
                  [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
        for k = 1:2
            b(k).FaceColor = 'flat';
            for c = 1:length(unique_conditions)
                b(k).CData(c,:) = colors{c};
            end
        end

        xlabel('Condition');
        ylabel(angle_name);
        title(sprintf('%s Across Conditions', angle_name), 'Interpreter', 'none');
        set(gca, 'XTick', 1:length(unique_conditions), 'XTickLabel', unique_conditions);
        box off;
        set(gca, 'TickDir', 'out');

        % Set y-axis limits
        if strcmp(angle_name, 'Knee')
            ylim([80 100]);
        elseif strcmp(angle_name, 'Ankle')
            ylim([90 130]);
        elseif strcmp(angle_name, 'Hip')
            ylim([50 70]);
        end
    end

    % Perform analysis per animal per condition
    animal_ids = cellfun(@(x) x(1:end-1), frame_identifiers, 'UniformOutput', false);
    unique_animals = unique(animal_ids);

    % Collect data for all angles and legs per animal per condition
    all_data_animal_condition = struct();
    for angle_idx = 1:n_angles
        angle_name = angle_names{angle_idx};
        all_data_animal_condition.(angle_name) = struct();
        for f = 1:length(fields)
            fieldname = fields{f};
            data = data_struct.(fieldname);
            angle_data = data(:, angle_idx);

            % Initialize cell array to store data for each animal per condition
            animal_condition_data = cell(length(unique_conditions), length(unique_animals));
            for c = 1:length(unique_conditions)
                condition = unique_conditions{c};
                condition_mask = strcmp(conditions(walking_bouts), condition);
                for a = 1:length(unique_animals)
                    animal = unique_animals{a};
                    animal_mask = strcmp(animal_ids(walking_bouts), animal);
                    combined_mask = condition_mask & animal_mask;
                    animal_condition_data{c, a} = angle_data(combined_mask, :);
                end
            end
            all_data_animal_condition.(angle_name).(fieldname) = animal_condition_data;
        end
    end    % Plotting per condition, averaging per animal (LEFT LEG ONLY)
    for angle_idx = 1:n_angles
        angle_name = angle_names{angle_idx};
        animal_condition_data_left = all_data_animal_condition.(angle_name).left_leg;

        % Calculate means and SEMs for each condition, averaging per animal (left leg only)
        means_left_condition = zeros(1, length(unique_conditions));
        sems_left_condition = zeros(1, length(unique_conditions));
        animal_means_matrix = zeros(length(unique_conditions), length(unique_animals)); % For ANOVA
        
        for c = 1:length(unique_conditions)
            animal_means_left = cellfun(@(x) nanmean(x), animal_condition_data_left(c, :));
            means_left_condition(c) = nanmean(animal_means_left);
            sems_left_condition(c) = nanstd(animal_means_left) / sqrt(sum(~isnan(animal_means_left)));
            animal_means_matrix(c, :) = animal_means_left;
        end

        % Prepare data for ANOVA
        anova_data = [];
        group_labels = [];
        for c = 1:length(unique_conditions)
            valid_data = animal_means_matrix(c, ~isnan(animal_means_matrix(c, :)));
            anova_data = [anova_data, valid_data];
            group_labels = [group_labels, repmat(c, 1, length(valid_data))];
        end        % Perform one-way ANOVA
        [p_anova, ~, stats] = anova1(anova_data, group_labels, 'off');
        
        % Perform Bonferroni post-hoc comparisons if ANOVA is significant
        p_bonferroni = [];
        if p_anova < 0.05
            [c_bonf, ~, ~, ~] = multcompare(stats, 'CType', 'tukey-kramer', 'Display', 'off');
            p_bonferroni = c_bonf(:, 6); % p-values from Bonferroni comparison
        end

        % Create figure for comparison
        figure('Position', [100 100 800 600], 'Color', 'w');

        % Create bar plot with error bars (left leg only)
        b = bar(means_left_condition, 'FaceColor', 'flat');
        hold on;
        errorbar(1:length(unique_conditions), means_left_condition, sems_left_condition, 'k', 'LineStyle', 'none', 'CapSize', 10);

        % Customize plot colors
        colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
                  [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
        for c = 1:length(unique_conditions)
            if c <= length(colors)
                b.CData(c,:) = colors{c};
            else
                b.CData(c,:) = [0.5 0.5 0.5]; % Default gray for extra conditions
            end
        end

        % Add significance stars for post-hoc comparisons
        if ~isempty(p_bonferroni)            % Get y-axis limits to position stars appropriately
            ylims = ylim;
            y_max = ylims(2);
            y_range = ylims(2) - ylims(1);
            star_height = y_max + 0.05 * y_range;
            
            % Counter for vertical positioning of multiple comparisons
            comparison_level = 0;
            
            for i = 1:size(c_bonf, 1)
                group1 = c_bonf(i, 1);
                group2 = c_bonf(i, 2);
                p_val = p_bonferroni(i);
                
                if p_val < 0.05
                    % Determine significance level
                    if p_val < 0.001
                        sig_text = '***';
                    elseif p_val < 0.01
                        sig_text = '**';
                    else
                        sig_text = '*';
                    end
                    
                    % Position for this comparison
                    current_star_height = star_height + comparison_level * 0.08 * y_range;
                    current_line_height = current_star_height + 0.02 * y_range;
                    
                    % Draw line connecting the groups
                    plot([group1, group2], [current_line_height, current_line_height], 'k-', 'LineWidth', 1);
                    plot([group1, group1], [current_line_height, current_line_height - 0.01 * y_range], 'k-', 'LineWidth', 1);
                    plot([group2, group2], [current_line_height, current_line_height - 0.01 * y_range], 'k-', 'LineWidth', 1);
                    
                    % Add significance text
                    text((group1 + group2) / 2, current_star_height, sig_text, ...
                         'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
                    
                    comparison_level = comparison_level + 1;
                end
            end
            
            % Adjust y-axis limits to accommodate significance markers
            if comparison_level > 0
                new_ylim_max = star_height + (comparison_level - 1) * 0.08 * y_range + 0.05 * y_range;
                ylim([ylims(1), new_ylim_max]);
            end
        end

        xlabel('Condition');
        ylabel([angle_name ' (Left Leg)']);
        title(sprintf('%s Across Conditions - Left Leg (p_{ANOVA} = %.4f)', angle_name, p_anova), 'Interpreter', 'none');
        set(gca, 'XTick', 1:length(unique_conditions), 'XTickLabel', unique_conditions);
        box off;
        set(gca, 'TickDir', 'out');

        % Set y-axis limits (modified to allow for significance markers)
        if strcmp(angle_name, 'Knee')
            ylim([80 105]); % Increased upper limit for significance markers
        elseif strcmp(angle_name, 'Ankle')
            ylim([90 140]); % Increased upper limit for significance markers
        elseif strcmp(angle_name, 'Hip')
            ylim([50 100]); % Increased upper limit for significance markers
        end

        if do_export_figure
            % Export folder
            export_folder = fullfile(GC.temp_root, 'figs_presentation_painAI');
            if ~exist(export_folder, 'dir')
                mkdir(export_folder);
            end
            fig_name = sprintf('angles_%s_leftleg', angle_name);
            exportgraphics(gcf, fullfile(export_folder, [fig_name,'.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
        end
            

    end
    
    disp('done')
end


%% This code is to check only.
% % Create a video reader object
% videoObj = VideoReader('K:\Mario\BioMed_students_2023\Anna\exp_6cam_miniscope_data\baseline\6cam_data\ID_1378\20240923\videos\Camera1\0.mp4');
% 
% % Get video properties
% numFrames = videoObj.NumFrames;
% 
% % Assuming walking_bouts is a logical array with same length as number of frames
% % If not already created, you'll need to generate this based on your analysis
% 
% % Initialize array to store extracted frames
% extractedFrames = {};
% frameIndices = [];
% 
% % Read and process frames
% for frameNum = 1:10000
%     % Read current frame
%     currentFrame = read(videoObj, frameNum);
% 
%     % Check if this frame has walking activity
% 
%     % Store the frame
%     extractedFrames{end+1} = currentFrame;
%     frameIndices = [frameIndices frameNum];
% 
% end
%figure
%for iv = 1:length(extractedFrames)
%    if walking_bouts(iv)
%        imshow(extractedFrames{iv});
%    end
%end
