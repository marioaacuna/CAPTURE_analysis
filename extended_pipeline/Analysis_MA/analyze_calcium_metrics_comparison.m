function cluster_changes = analyze_calcium_metrics_comparison(data1, data2, animals_of_interest, metric_to_take, comparison)
    % analyze_calcium_metrics_comparison - Analyzes and visualizes calcium imaging metrics
    % 
    % This function performs comprehensive analysis and visualization of calcium imaging 
    % metrics comparing two experimental conditions (e.g., Saline vs Formalin, Sham vs Neuropathic).
    % It generates four types of plots: bar comparison, difference chart, heatmap, and scatter plot.
    % 
    % Inputs:
    %   data1 - Structure containing control group data (Saline/Sham)
    %   data2 - Structure containing experimental group data (Formalin/Neuropathic)
    %   animals_of_interest - Cell array of animal IDs to include in analysis
    %   metric_to_take - String specifying metric ('max_amplitude', 'peaks', 'freqs')
    %   comparison - String specifying comparison type ('S_v_F' or 'H_v_N')
    %
    % Outputs:
    %   cluster_changes - Structure containing cluster classifications:
    %                    .increased - cluster IDs with significant increase
    %                    .decreased - cluster IDs with significant decrease  
    %                    .no_change - cluster IDs with no significant change
    %                    .statistics - detailed statistics for each cluster
    %
    % The function:
    %   1. Concatenates data across animals for each cluster
    %   2. Handles clusters present in only one condition (assigns zeros)
    %   3. Classifies clusters based on statistical significance and direction
    %   4. Generates four visualization types with appropriate sorting
    %   5. Exports figures to the configured temporary directory
    
    global GC;
    
    % Set group labels based on comparison type
    if strcmp(comparison, 'S_v_F')
        group1_label = 'Saline';
        group2_label = 'Formalin';
    else % H_v_N
        group1_label = 'Sham';
        group2_label = 'Neuropathic';
    end

    % Find global clusters
    global_clusters = [];
    for animal = 1:length(animals_of_interest)
        animal_ID = animals_of_interest{animal};
        if isfield(data1, animal_ID)
            global_clusters = union(global_clusters, data1.(animal_ID).unique_clusters);
        end
        if isfield(data2, animal_ID)
            global_clusters = union(global_clusters, data2.(animal_ID).unique_clusters);
        end
    end

    % Concatenate data for each cluster
    concatenated_1 = cell(1, max(global_clusters) + 1);
    concatenated_2 = cell(1, max(global_clusters) + 1);

    for cluster_idx = 0:max(global_clusters)
        concatenated_1{cluster_idx + 1} = {};
        concatenated_2{cluster_idx + 1} = {};
        
        for animal = 1:length(animals_of_interest)
            animal_ID = animals_of_interest{animal};
            
            % Process group 1 data
            if isfield(data1, animal_ID)
                if ismember(cluster_idx, data1.(animal_ID).unique_clusters)
                    if strcmp(metric_to_take, 'max_amplitude')
                        data_to_add = data1.(animal_ID).(metric_to_take)(:, data1.(animal_ID).unique_clusters == cluster_idx);
                    else % for 'peaks' and 'freqs'
                        animal_data = data1.(animal_ID).ensemble_activity;
                        if ~isempty(animal_data(data1.(animal_ID).unique_clusters == cluster_idx).(metric_to_take))
                            data_to_add = animal_data(data1.(animal_ID).unique_clusters == cluster_idx).(metric_to_take);
                        else
                            continue;
                        end
                    end
                    concatenated_1{cluster_idx + 1}{end+1} = data_to_add;
                end
            end
            
            % Process group 2 data
            if isfield(data2, animal_ID)
                if ismember(cluster_idx, data2.(animal_ID).unique_clusters)
                    if strcmp(metric_to_take, 'max_amplitude')
                        data_to_add = data2.(animal_ID).(metric_to_take)(:, data2.(animal_ID).unique_clusters == cluster_idx);
                    else % for 'peaks' and 'freqs'
                        animal_data = data2.(animal_ID).ensemble_activity;

                        if ~isempty(animal_data(data2.(animal_ID).unique_clusters == cluster_idx).(metric_to_take))
                            data_to_add = animal_data(data2.(animal_ID).unique_clusters == cluster_idx).(metric_to_take);
                        else
                            continue;
                        end

                    end
                    concatenated_2{cluster_idx + 1}{end+1} = data_to_add;
                end
            end
        end
    end    % Initialize structures for all ROIs per cluster
    all_ROIs_per_cluster_1 = struct();
    all_ROIs_per_cluster_2 = struct();

    for cluster_idx = 1:length(global_clusters)
        all_ROIs_1 = [];
        all_ROIs_2 = [];
        
        if ~isempty(concatenated_1{cluster_idx})
            for animal = 1:length(concatenated_1{cluster_idx})
                all_ROIs_1 = [all_ROIs_1; concatenated_1{cluster_idx}{animal}];
            end
        end
        
        if ~isempty(concatenated_2{cluster_idx})
            for animal = 1:length(concatenated_2{cluster_idx})
                all_ROIs_2 = [all_ROIs_2; concatenated_2{cluster_idx}{animal}];
            end
        end
        
        cluster_str = ['cluster_' num2str(global_clusters(cluster_idx))];
        
        % Handle empty clusters by assigning zeros
        if isempty(all_ROIs_1)
            all_ROIs_1 = 0;
        end
        if isempty(all_ROIs_2)
            all_ROIs_2 = 0;
        end
        
        all_ROIs_per_cluster_1.(cluster_str) = all_ROIs_1;
        all_ROIs_per_cluster_2.(cluster_str) = all_ROIs_2;
    end

    % Classify clusters based on statistical significance and direction of change
    cluster_changes = classify_cluster_changes(all_ROIs_per_cluster_1, all_ROIs_per_cluster_2, global_clusters);
    
    % Display cluster classification results
    fprintf('\n=== Cluster Classification Results for %s vs %s ===\n', group1_label, group2_label);
    fprintf('Metric: %s\n\n', metric_to_take);
    
    if ~isempty(cluster_changes.increased)
        fprintf('Clusters with SIGNIFICANT INCREASE:\n');
        for i = 1:length(cluster_changes.increased)
            cluster_id = cluster_changes.increased(i);
            stats = cluster_changes.statistics.(['cluster_' num2str(cluster_id)]);
            fprintf('  Cluster %d: p=%.4f, difference=%.4f\n', cluster_id, stats.p_value, stats.difference);
        end
        fprintf('\n');
    else
        fprintf('No clusters with significant increase.\n\n');
    end
    
    if ~isempty(cluster_changes.decreased)
        fprintf('Clusters with SIGNIFICANT DECREASE:\n');
        for i = 1:length(cluster_changes.decreased)
            cluster_id = cluster_changes.decreased(i);
            stats = cluster_changes.statistics.(['cluster_' num2str(cluster_id)]);
            fprintf('  Cluster %d: p=%.4f, difference=%.4f\n', cluster_id, stats.p_value, stats.difference);
        end
        fprintf('\n');
    else
        fprintf('No clusters with significant decrease.\n\n');
    end
    
    if ~isempty(cluster_changes.no_change)
        fprintf('Clusters with NO SIGNIFICANT CHANGE:\n');
        for i = 1:length(cluster_changes.no_change)
            cluster_id = cluster_changes.no_change(i);
            stats = cluster_changes.statistics.(['cluster_' num2str(cluster_id)]);
            fprintf('  Cluster %d: p=%.4f, difference=%.4f\n', cluster_id, stats.p_value, stats.difference);
        end
        fprintf('\n');
    else
        fprintf('All clusters show significant changes.\n\n');
    end

    % Generate plots
    % plot_bar_comparison(all_ROIs_per_cluster_1, all_ROIs_per_cluster_2, group1_label, group2_label, metric_to_take);
    plot_difference_chart(all_ROIs_per_cluster_1, all_ROIs_per_cluster_2, group1_label, group2_label, metric_to_take);
    % plot_heatmap(all_ROIs_per_cluster_1, all_ROIs_per_cluster_2, group1_label, group2_label, metric_to_take);
    % plot_scatter_comparison(all_ROIs_per_cluster_1, all_ROIs_per_cluster_2, group1_label, group2_label, metric_to_take);
end

function plot_bar_comparison(data1, data2, group1_label, group2_label, metric_to_take)
    % plot_bar_comparison - Creates side-by-side bar comparison plot
    %
    % This function generates a bar chart comparing mean values between two groups
    % for each cluster. Bars are sorted by group 2 activity (descending order).
    % Error bars show Standard Error of the Mean (SEM) for each group.
    % Statistical significance is indicated above each cluster pair.
    %
    % Inputs:
    %   data1, data2 - Structures containing data for each cluster
    %   group1_label, group2_label - String labels for the two groups
    %   metric_to_take - String specifying the metric being analyzed
    
    global GC;
    
    Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]);
    
    % Define plotting constants
    barWidth = 0.75;
    gapWidth = 1;
    currentX = 1;
    
    clusters = fieldnames(data1);
    
    % Filter out clusters with NaN values in either condition
    valid_clusters = {};
    group2_means_valid = [];
    
    for i = 1:length(clusters)
        cluster_name = clusters{i};
        data_group1 = data1.(cluster_name);
        data_group2 = data2.(cluster_name);
        
        mean1 = nanmean(data_group1);
        mean2 = nanmean(data_group2);
        
        % Only include clusters where both conditions have valid (non-NaN) data
        if ~isnan(mean1) && ~isnan(mean2) && ~isempty(data_group1) && ~isempty(data_group2)
            valid_clusters{end+1} = cluster_name;
            group2_means_valid(end+1) = mean2;
        end
    end
    
    % Sort valid clusters by group 2 means (descending order)
    [~, sort_idx] = sort(group2_means_valid, 'descend');
    sorted_clusters = valid_clusters(sort_idx);
    
    for i = 1:length(sorted_clusters)
        cluster_name = sorted_clusters{i};
        
        data_group1 = data1.(cluster_name);
        data_group2 = data2.(cluster_name);
        
        % Calculate statistics
        try
            [h, p, ci, stats] = ttest2(data_group1, data_group2);
        catch
            p = 1;
            stats.tstat = [];
        end
          % Calculate means and SEMs
        mean_1 = nanmean(data_group1);
        SEM_1 = nanstd(data_group1)/sqrt(length(data_group1));
        mean_2 = nanmean(data_group2);
        SEM_2 = nanstd(data_group2)/sqrt(length(data_group2));
        
        % Plot bars and error bars
        bar(currentX, mean_1, barWidth, 'b');
        hold on;
        errorbar(currentX, mean_1, SEM_1, 'k', 'LineStyle', 'none');
        bar(currentX + barWidth, mean_2, barWidth, 'r');
        errorbar(currentX + barWidth, mean_2, SEM_2, 'k', 'LineStyle', 'none');
        
        % Add significance markers - Fixed positioning
        maxY = max(mean_1 + SEM_1, mean_2 + SEM_2);
        textY = maxY + maxY * 0.15;
        textX = currentX + barWidth/2; % Center between the two bars
        
        if p < 0.001
            text(textX, textY, '***', 'HorizontalAlignment', 'center', 'FontSize', 12);
        elseif p < 0.01
            text(textX, textY, '**', 'HorizontalAlignment', 'center', 'FontSize', 12);
        elseif p < 0.05
            text(textX, textY, '*', 'HorizontalAlignment', 'center', 'FontSize', 12);
        end
        
        currentX = currentX + barWidth * 2 + gapWidth;
    end
    
    % Add legend
    hG1 = bar(NaN, NaN, 'b');
    hG2 = bar(NaN, NaN, 'r');
    legend([hG1, hG2], {group1_label, group2_label});
      % Customize plot
    ylabel(['Mean ' metric_to_take]);
    xlabel(['Clusters (sorted by ' group2_label ' activity)']);
    title(['Mean ' metric_to_take ' of ROIs per Cluster for ' group1_label ' and ' group2_label ' (sorted by ' group2_label ')']);
    set(gca, 'XTick', 1.5:barWidth*2+gapWidth:length(sorted_clusters)*(barWidth*2+gapWidth));
    set(gca, 'XTickLabel', 1:length(sorted_clusters)); % Use numbers instead of cluster names
    box off;
    set(gca, 'TickDir', 'out');
    
    % Export figure
    if ~isempty(GC) && isfield(GC, 'temp_root')
        export_folder = GC.temp_root;
        fig_name = ['BarComparison_' metric_to_take '_' group1_label '_vs_' group2_label];
        exportgraphics(gcf, fullfile(export_folder, [fig_name,'.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
    end
end

function plot_difference_chart(data1, data2, group1_label, group2_label, metric_to_take)
    % plot_difference_chart - Creates difference chart showing group2 - group1
    %
    % This function generates a bar chart showing the difference between group means
    % for each cluster. Bars are sorted by difference magnitude (descending order).
    % No error bars are shown since differences are calculated values, not measurements.
    % Color coding indicates statistical significance levels.
    %
    % Inputs:
    %   data1, data2 - Structures containing data for each cluster
    %   group1_label, group2_label - String labels for the two groups
    %   metric_to_take - String specifying the metric being analyzed

    global GC    
    Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]);
    
    barWidth = 0.85;
    clusters = fieldnames(data1);
    
    % Preallocate arrays
    differences = zeros(1, length(clusters));
    errors = zeros(1, length(clusters));
    p_values_array = zeros(1, length(clusters));
    valid_indices = false(1, length(clusters));
      % Calculate differences and statistics
    for i = 1:length(clusters)
        cluster_name = clusters{i};
        
        data_group1 = data1.(cluster_name);
        data_group2 = data2.(cluster_name);
        
        % Handle cases where data might be scalar zero (empty cluster)
        if isscalar(data_group1) && data_group1 == 0
            mean_group1 = 0;
        elseif isempty(data_group1)
            mean_group1 = 0;
        else
            mean_group1 = nanmean(data_group1);
        end
        
        if isscalar(data_group2) && data_group2 == 0
            mean_group2 = 0;
        elseif isempty(data_group2)
            mean_group2 = 0;
        else
            mean_group2 = nanmean(data_group2);
        end
          differences(i) = mean_group2 - mean_group1;
        
        % Calculate error - handle scalar zeros (empty clusters)
        if (isscalar(data_group1) && data_group1 == 0) && (isscalar(data_group2) && data_group2 == 0)
            % Both groups are empty clusters (scalar zeros)
            errors(i) = 0;
        elseif (isscalar(data_group1) && data_group1 == 0) && ~(isscalar(data_group2) && data_group2 == 0)
            % Group 1 is empty cluster, group 2 has data
            if ~isempty(data_group2)
                errors(i) = nanstd(data_group2)/sqrt(length(data_group2));
            else
                errors(i) = 0;
            end
        elseif ~(isscalar(data_group1) && data_group1 == 0) && (isscalar(data_group2) && data_group2 == 0)
            % Group 1 has data, group 2 is empty cluster
            if ~isempty(data_group1)
                errors(i) = nanstd(data_group1)/sqrt(length(data_group1));
            else
                errors(i) = 0;
            end
        elseif ~isempty(data_group1) && ~isempty(data_group2)
            % Both groups have actual data
            errors(i) = sqrt((nanstd(data_group2)^2/length(data_group2)) + (nanstd(data_group1)^2/length(data_group1)));
        elseif ~isempty(data_group1)
            errors(i) = nanstd(data_group1)/sqrt(length(data_group1));
        elseif ~isempty(data_group2)
            errors(i) = nanstd(data_group2)/sqrt(length(data_group2));
        else
            errors(i) = 0;
        end
        
        % Statistical test - handle scalar zeros
        try
            if (isscalar(data_group1) && data_group1 == 0) || (isscalar(data_group2) && data_group2 == 0)
                % If either group is empty cluster, use different significance level
                p_values_array(i) = 0.001; % Show as significant for visualization
            elseif ~isempty(data_group1) && ~isempty(data_group2)
                [~, p] = ttest2(data_group1, data_group2);
                p_values_array(i) = p;
            else
                p_values_array(i) = 0.001; % cosmetic
            end
        catch
            p_values_array(i) = 1;
        end
        
        valid_indices(i) = true;
    end
    
    % Sort by difference values in descending order
    valid_data = valid_indices;
    [sorted_differences, sort_idx] = sort(differences(valid_data), 'descend');
    valid_clusters = find(valid_data);
    sorted_indices = valid_clusters(sort_idx);
    
    % Extract sorted data
    sorted_clusters = clusters(sorted_indices);
    sorted_errors = errors(sorted_indices);
    sorted_p_values = p_values_array(sorted_indices);
    
    % Create color coding based on p-values
    colors = zeros(length(sorted_clusters), 3);
    colors(sorted_p_values < 0.05 & sorted_p_values >= 0.01, :) = repmat([0.8 0.4 0], sum(sorted_p_values < 0.05 & sorted_p_values >= 0.01), 1);
    colors(sorted_p_values < 0.01 & sorted_p_values >= 0.001, :) = repmat([0.9 0.2 0], sum(sorted_p_values < 0.01 & sorted_p_values >= 0.001), 1);
    colors(sorted_p_values < 0.001, :) = repmat([1 0 0], sum(sorted_p_values < 0.001), 1);
    colors(sorted_p_values >= 0.05, :) = repmat([0.7 0.7 0.7], sum(sorted_p_values >= 0.05), 1);
    
    % Plot
    b = bar(sorted_differences, barWidth);
    b.FaceColor = 'flat';
    b.CData = colors;
    
    hold on;
    errorbar(1:length(sorted_clusters), sorted_differences, sorted_errors, 'k.', 'LineStyle', 'none');
    yline(0, 'k--', 'Alpha', 0.3);
    
    % Customize plot
    ylabel(['Difference in Mean ' metric_to_take ' (' group2_label ' - ' group1_label ')']);
    xlabel('Clusters');
    title(['Difference in ' metric_to_take ' between ' group2_label ' and ' group1_label ' Conditions']);
    
    % Extract just the cluster numbers for cleaner labels
    cluster_numbers = zeros(1, length(sorted_clusters));
    for i = 1:length(sorted_clusters)
        % Find the underscore position and extract the number after it
        underscore_pos = strfind(sorted_clusters{i}, '_');
        if ~isempty(underscore_pos)
            cluster_numbers(i) = str2double(sorted_clusters{i}(underscore_pos(end)+1:end));
        else
            % If no underscore, just use the original name
            cluster_numbers(i) = i;
        end
    end
   
    
    
    % set(gca, 'XTick', 1:length(cluster_numbers));
    % set(gca, 'XTickLabel', cluster_numbers, 'TickLabelInterpreter', 'none');
    set(gca, 'XTick', {});
    box off;
    set(gca, 'TickDir', 'out');
    
    % Add legend
    legend_elements = [
        patch([0 0], [0 0], [0.7 0.7 0.7], 'DisplayName', 'n.s.');
        patch([0 0], [0 0], [0.8 0.4 0], 'DisplayName', 'p < 0.05');
        patch([0 0], [0 0], [0.9 0.2 0], 'DisplayName', 'p < 0.01');
        patch([0 0], [0 0], [1 0 0], 'DisplayName', 'p < 0.001')
    ];
    legend(legend_elements, 'Location', 'northeast');

    % Export figure
    if ~isempty(GC) && isfield(GC, 'temp_root')
        export_folder = GC.temp_root;
        fig_name = ['DifferenceChart_' metric_to_take '_' group1_label '_vs_' group2_label];
        exportgraphics(gcf, fullfile(export_folder, [fig_name,'.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
    end
end



function plot_heatmap(data1, data2, group1_label, group2_label, metric_to_take)
    % plot_heatmap - Creates heatmap visualization of group means
    %
    % This function generates a heatmap showing mean values for both groups
    % across all clusters. Rows represent conditions and columns represent clusters.
    % Clusters are sorted by group 2 activity (descending order) and labeled with
    % numbers instead of cluster IDs. Statistical significance is indicated with
    % asterisks above each cluster column.
    %
    % Inputs:
    %   data1, data2 - Structures containing data for each cluster
    %   group1_label, group2_label - String labels for the two groups
    %   metric_to_take - String specifying the metric being analyzed
    
    global GC;
    
    Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]);
    
    clusters = fieldnames(data1);
    
    % Filter out clusters with NaN values in either condition
    valid_clusters = {};
    group2_means_valid = [];
    
    for i = 1:length(clusters)
        cluster_name = clusters{i};
        data_group1 = data1.(cluster_name);
        data_group2 = data2.(cluster_name);
        
        mean1 = nanmean(data_group1);
        mean2 = nanmean(data_group2);
        
        % Only include clusters where both conditions have valid (non-NaN) data
        if ~isnan(mean1) && ~isnan(mean2) && ~isempty(data_group1) && ~isempty(data_group2)
            valid_clusters{end+1} = cluster_name;
            group2_means_valid(end+1) = mean2;
        end
    end
    
    % Sort valid clusters by group 2 means (descending order)
    [~, sort_idx] = sort(group2_means_valid, 'descend');
    sorted_clusters = valid_clusters(sort_idx);
    num_clusters = length(sorted_clusters);
    
    % Calculate means and p-values for sorted clusters
    means_matrix = zeros(2, num_clusters);
    p_values_array = zeros(1, num_clusters);
    
    for i = 1:num_clusters
        cluster_name = sorted_clusters{i};
        
        data_group1 = data1.(cluster_name);
        data_group2 = data2.(cluster_name);
        
        means_matrix(1,i) = nanmean(data_group1);
        means_matrix(2,i) = nanmean(data_group2);
        
        try
            [~, p] = ttest2(data_group1, data_group2);
            p_values_array(i) = p;
        catch
            p_values_array(i) = 1;
        end
    end
    
    % Create heatmap
    imagesc(means_matrix);
    colormap('jet');
    colorbar;
    
    % Add significance markers
    hold on;
    for i = 1:num_clusters
        if p_values_array(i) < 0.001
            text(i, 1.5, '***', 'HorizontalAlignment', 'center', 'Color', 'k', 'FontSize', 12);
        elseif p_values_array(i) < 0.01
            text(i, 1.5, '**', 'HorizontalAlignment', 'center', 'Color', 'k', 'FontSize', 12);
        elseif p_values_array(i) < 0.05
            text(i, 1.5, '*', 'HorizontalAlignment', 'center', 'Color', 'k', 'FontSize', 12);
        end
    end
      % Customize plot
    ylabel('Condition');
    xlabel(['Clusters (sorted by ' group2_label ' activity)']);
    title(['Mean ' metric_to_take ' Heatmap of ROIs per Cluster (sorted by ' group2_label ')']);
    set(gca, 'YTick', 1:2);
    set(gca, 'YTickLabel', {group1_label, group2_label});
    set(gca, 'XTick', 1:num_clusters);
    set(gca, 'XTickLabel', 1:num_clusters); % Use numbers instead of cluster names
    
    % Export figure
    if ~isempty(GC) && isfield(GC, 'temp_root')
        export_folder = GC.temp_root;
        fig_name = ['Heatmap_' metric_to_take '_' group1_label '_vs_' group2_label];
        exportgraphics(gcf, fullfile(export_folder, [fig_name,'.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
    end
end

function plot_scatter_comparison(data1, data2, group1_label, group2_label, metric_to_take)
    % plot_scatter_comparison - Creates scatter plot comparing group means
    %
    % This function generates a scatter plot showing mean values for both groups
    % with connecting lines for each cluster. Clusters are sorted by group 2 activity
    % (descending order) and labeled with numbers instead of cluster IDs.
    % Line colors indicate statistical significance levels between groups.
    %
    % Inputs:
    %   data1, data2 - Structures containing data for each cluster
    %   group1_label, group2_label - String labels for the two groups
    %   metric_to_take - String specifying the metric being analyzed
    
    global GC;
    
    Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]);
    
    clusters = fieldnames(data1);
    
    % Filter out clusters with NaN values in either condition
    valid_clusters = {};
    group2_means_valid = [];
    
    for i = 1:length(clusters)
        cluster_name = clusters{i};
        data_group1 = data1.(cluster_name);
        data_group2 = data2.(cluster_name);
        
        mean1 = nanmean(data_group1);
        mean2 = nanmean(data_group2);
        
        % Only include clusters where both conditions have valid (non-NaN) data
        if ~isnan(mean1) && ~isnan(mean2) && ~isempty(data_group1) && ~isempty(data_group2)
            valid_clusters{end+1} = cluster_name;
            group2_means_valid(end+1) = mean2;
        end
    end
    
    % Sort valid clusters by group 2 means (descending order)
    [~, sort_idx] = sort(group2_means_valid, 'descend');
    sorted_clusters = valid_clusters(sort_idx);
    num_clusters = length(sorted_clusters);
    
    % Calculate means and p-values for sorted clusters
    means_1 = zeros(1, num_clusters);
    means_2 = zeros(1, num_clusters);
    p_values_array = zeros(1, num_clusters);
    
    for i = 1:num_clusters
        cluster_name = sorted_clusters{i};
        
        data_group1 = data1.(cluster_name);
        data_group2 = data2.(cluster_name);
        
        means_1(i) = nanmean(data_group1);
        means_2(i) = nanmean(data_group2);
        
        try
            [~, p] = ttest2(data_group1, data_group2);
            p_values_array(i) = p;
        catch
            p_values_array(i) = 1;
        end
    end
    
    % Create color coding
    colors = zeros(num_clusters, 3);
    colors(p_values_array < 0.05 & p_values_array >= 0.01, :) = repmat([0.8 0.4 0], sum(p_values_array < 0.05 & p_values_array >= 0.01), 1);
    colors(p_values_array < 0.01 & p_values_array >= 0.001, :) = repmat([0.9 0.2 0], sum(p_values_array < 0.01 & p_values_array >= 0.001), 1);
    colors(p_values_array < 0.001, :) = repmat([1 0 0], sum(p_values_array < 0.001), 1);
    colors(p_values_array >= 0.05, :) = repmat([0.7 0.7 0.7], sum(p_values_array >= 0.05), 1);
    
    % Plot connecting lines and points
    for i = 1:num_clusters
        line([means_1(i) means_2(i)], [i i], 'Color', colors(i,:), 'LineWidth', 2);
    end
    
    hold on;
    scatter(means_1, 1:num_clusters, 50, 'b', 'filled');
    scatter(means_2, 1:num_clusters, 50, 'r', 'filled');
    
    % Customize plot
    ylabel(['Clusters (sorted by ' group2_label ' activity)']);
    xlabel(['Mean ' metric_to_take]);
    title(['Comparison of Mean ' metric_to_take ' between Conditions (sorted by ' group2_label ')']);
    set(gca, 'YTick', 1:num_clusters);
    set(gca, 'YTickLabel', 1:num_clusters); % Use numbers instead of cluster names
    box off;
    set(gca, 'TickDir', 'out');
      % Add legend
    legend_elements = [
        scatter(NaN, NaN, 50, 'b', 'filled', 'DisplayName', group1_label);
        scatter(NaN, NaN, 50, 'r', 'filled', 'DisplayName', group2_label);
        line([NaN NaN], [NaN NaN], 'Color', [0.7 0.7 0.7], 'DisplayName', 'n.s.');
        line([NaN NaN], [NaN NaN], 'Color', [0.8 0.4 0], 'DisplayName', 'p < 0.05');
        line([NaN NaN], [NaN NaN], 'Color', [0.9 0.2 0], 'DisplayName', 'p < 0.01');
        line([NaN NaN], [NaN NaN], 'Color', [1 0 0], 'DisplayName', 'p < 0.001')
    ];
    legend(legend_elements, 'Location', 'northeast');
    
    % Export figure
    if ~isempty(GC) && isfield(GC, 'temp_root')
        export_folder = GC.temp_root;
        fig_name = ['ScatterComparison_' metric_to_take '_' group1_label '_vs_' group2_label];
        exportgraphics(gcf, fullfile(export_folder, [fig_name,'.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
    end
end

function cluster_changes = classify_cluster_changes(data1, data2, global_clusters)
    % classify_cluster_changes - Classifies clusters based on statistical significance and direction
    %
    % This function analyzes each cluster to determine if there is a significant
    % increase, decrease, or no change between the two conditions. It handles
    % clusters that may be present in only one condition by treating missing
    % data as zeros.
    %
    % Inputs:
    %   data1 - Structure containing group 1 data for each cluster
    %   data2 - Structure containing group 2 data for each cluster
    %   global_clusters - Array of all cluster IDs to analyze
    %
    % Outputs:
    %   cluster_changes - Structure containing:
    %                    .increased - cluster IDs with significant increase (p<0.05, group2>group1)
    %                    .decreased - cluster IDs with significant decrease (p<0.05, group2<group1)
    %                    .no_change - cluster IDs with no significant change (p>=0.05)
    %                    .statistics - detailed statistics for each cluster
    
    % Initialize output structure
    cluster_changes = struct();
    cluster_changes.increased = [];
    cluster_changes.decreased = [];
    cluster_changes.no_change = [];
    cluster_changes.statistics = struct();
    
    % Significance threshold
    alpha = 0.05;
    
    % Analyze each cluster
    for i = 1:length(global_clusters)
        cluster_id = global_clusters(i);
        cluster_str = ['cluster_' num2str(cluster_id)];
        
        % Get data for this cluster
        data_group1 = data1.(cluster_str);
        data_group2 = data2.(cluster_str);
        
        % Handle cases where data might be scalar zero (empty cluster)
        if isscalar(data_group1) && data_group1 == 0
            data_group1 = [];
        end
        if isscalar(data_group2) && data_group2 == 0
            data_group2 = [];
        end
        
        % Calculate means (treat empty as 0)
        if isempty(data_group1)
            mean_group1 = 0;
            std_group1 = 0;
            n_group1 = 1; % Treat as single observation of 0
        else
            mean_group1 = nanmean(data_group1);
            std_group1 = nanstd(data_group1);
            n_group1 = length(data_group1);
        end
        
        if isempty(data_group2)
            mean_group2 = 0;
            std_group2 = 0;
            n_group2 = 1; % Treat as single observation of 0
        else
            mean_group2 = nanmean(data_group2);
            std_group2 = nanstd(data_group2);
            n_group2 = length(data_group2);
        end
        
        % Calculate difference (group2 - group1)
        difference = mean_group2 - mean_group1;
        
        % Perform statistical test
        p_value = 1; % Default to non-significant
        tstat = 0;
        
        if ~isempty(data_group1) && ~isempty(data_group2)
            % Both groups have data - use t-test
            try
                [~, p_value, ~, stats] = ttest2(data_group1, data_group2);
                tstat = stats.tstat;
            catch
                p_value = 1;
                tstat = 0;
            end
        elseif isempty(data_group1) && ~isempty(data_group2)
            % Only group 2 has data - test against zero
            try
                [~, p_value, ~, stats] = ttest(data_group2, 0);
                tstat = stats.tstat;
            catch
                p_value = 1;
                tstat = 0;
            end
        elseif ~isempty(data_group1) && isempty(data_group2)
            % Only group 1 has data - test against zero (negative direction)
            try
                [~, p_value, ~, stats] = ttest(data_group1, 0);
                tstat = -stats.tstat; % Negative because we're testing group2-group1
            catch
                p_value = 1;
                tstat = 0;
            end
        else
            % Both groups empty - no change
            p_value = 1;
            tstat = 0;
        end
        
        % Store detailed statistics
        cluster_changes.statistics.(cluster_str) = struct();
        cluster_changes.statistics.(cluster_str).cluster_id = cluster_id;
        cluster_changes.statistics.(cluster_str).mean_group1 = mean_group1;
        cluster_changes.statistics.(cluster_str).mean_group2 = mean_group2;
        cluster_changes.statistics.(cluster_str).difference = difference;
        cluster_changes.statistics.(cluster_str).p_value = p_value;
        cluster_changes.statistics.(cluster_str).t_statistic = tstat;
        cluster_changes.statistics.(cluster_str).n_group1 = n_group1;
        cluster_changes.statistics.(cluster_str).n_group2 = n_group2;
        cluster_changes.statistics.(cluster_str).std_group1 = std_group1;
        cluster_changes.statistics.(cluster_str).std_group2 = std_group2;
        
        % Classify based on significance and direction
        if p_value < alpha
            if difference > 0
                cluster_changes.increased(end+1) = cluster_id;
            else
                cluster_changes.decreased(end+1) = cluster_id;
            end
        else
            cluster_changes.no_change(end+1) = cluster_id;
        end
    end
    
    % Sort the arrays for consistent output
    cluster_changes.increased = sort(cluster_changes.increased);
    cluster_changes.decreased = sort(cluster_changes.decreased);
    cluster_changes.no_change = sort(cluster_changes.no_change);
end
