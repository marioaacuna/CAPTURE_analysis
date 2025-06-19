function [cluster_changes, plot_order] = analyze_calcium_metrics_comparison(data1, data2, animals_of_interest, metric_to_take, comparison, use_common_only)
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
    %   use_common_only - Boolean: true = only common clusters, false = all clusters (assigns zeros to missing)
    %
    % Outputs:
    %   cluster_changes - Structure containing cluster classifications:
    %                    .increased - cluster IDs with significant increase
    %                    .decreased - cluster IDs with significant decrease  
    %                    .no_change - cluster IDs with no significant change
    %                    .statistics - detailed statistics for each cluster
    %   plot_order - Vector of cluster IDs in the order they appear in the difference plot
    %               (sorted from highest increase to highest decrease)
    %
    % The function:
    %   1. Concatenates data across animals for each cluster
    %   2. Either uses only common clusters or handles missing clusters with zeros
    %   3. Classifies clusters based on statistical significance and direction
    %   4. Generates visualization with appropriate sorting
    %   5. Returns cluster order as shown in the plot
      global GC;
    
    % Set default for use_common_only if not provided
    if nargin < 6
        use_common_only = true; % Default to common clusters only
    end
    
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
    common_clusters = [];  % Track clusters to analyze

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
        
        if use_common_only
            % Only include clusters that have data in BOTH conditions
            if ~isempty(all_ROIs_1) && ~isempty(all_ROIs_2)
                all_ROIs_per_cluster_1.(cluster_str) = all_ROIs_1;
                all_ROIs_per_cluster_2.(cluster_str) = all_ROIs_2;
                common_clusters(end+1) = global_clusters(cluster_idx);
            end
        else
            % Include all clusters, assign zeros to empty ones
            if isempty(all_ROIs_1)
                all_ROIs_1 = 0;
            end
            if isempty(all_ROIs_2)
                all_ROIs_2 = 0;
            end
            all_ROIs_per_cluster_1.(cluster_str) = all_ROIs_1;
            all_ROIs_per_cluster_2.(cluster_str) = all_ROIs_2;
            common_clusters(end+1) = global_clusters(cluster_idx);
        end
    end
    
    % Update global_clusters to only include analyzed clusters
    global_clusters = common_clusters;
    
    if use_common_only
        fprintf('Analysis will include %d clusters present in both conditions.\n', length(global_clusters));
    else
        fprintf('Analysis will include %d clusters (empty clusters assigned zeros).\n', length(global_clusters));
    end
    
    if ~isempty(global_clusters)
        fprintf('Analyzed cluster IDs: %s\n', mat2str(global_clusters));
    else
        fprintf('Warning: No clusters found for analysis!\n');
        cluster_changes = struct('increased', [], 'decreased', [], 'no_change', [], 'statistics', struct());
        plot_order = [];
        return;
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
    end    % Generate plots
    % plot_bar_comparison(all_ROIs_per_cluster_1, all_ROIs_per_cluster_2, group1_label, group2_label, metric_to_take);
    plot_order = plot_difference_chart(all_ROIs_per_cluster_1, all_ROIs_per_cluster_2, group1_label, group2_label, metric_to_take);
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
    
    % All clusters should be valid since we only analyze common clusters
    group2_means_valid = [];
    
    for i = 1:length(clusters)
        cluster_name = clusters{i};
        data_group2 = data2.(cluster_name);
        mean2 = nanmean(data_group2);
        group2_means_valid(end+1) = mean2;
    end
    
    % Sort clusters by group 2 means (descending order)
    [~, sort_idx] = sort(group2_means_valid, 'descend');
    sorted_clusters = clusters(sort_idx);
    
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

function plot_order = plot_difference_chart(data1, data2, group1_label, group2_label, metric_to_take)
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
    %
    % Outputs:
    %   plot_order - Vector of cluster IDs in the order they appear in the plot
    %               (sorted from highest increase to highest decrease)

    global GC    
    Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]);
    
    barWidth = 0.85;
    clusters = fieldnames(data1);
    
    % Preallocate arrays
    differences = zeros(1, length(clusters));
    errors = zeros(1, length(clusters));
    p_values_array = zeros(1, length(clusters));
    valid_indices = false(1, length(clusters));    % Calculate differences and statistics
    for i = 1:length(clusters)
        cluster_name = clusters{i};
        
        data_group1 = data1.(cluster_name);
        data_group2 = data2.(cluster_name);
        
        % Calculate means (both groups should have data since we only analyze common clusters)
        mean_group1 = nanmean(data_group1);
        mean_group2 = nanmean(data_group2);
        
        differences(i) = mean_group2 - mean_group1;
        
        % Calculate error for both groups (both should have data)
        errors(i) = sqrt((nanstd(data_group2)^2/length(data_group2)) + (nanstd(data_group1)^2/length(data_group1)));
        
        % Statistical test
        try
            [~, p] = ttest2(data_group1, data_group2);
            p_values_array(i) = p;
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
    legend(legend_elements, 'Location', 'northeast');    % Export figure
    if ~isempty(GC) && isfield(GC, 'temp_root')
        export_folder = GC.temp_root;
        fig_name = ['DifferenceChart_' metric_to_take '_' group1_label '_vs_' group2_label];
        exportgraphics(gcf, fullfile(export_folder, [fig_name,'.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'none');
    end
    
    % Return cluster IDs in plot order (sorted from highest increase to highest decrease)
    plot_order = cluster_numbers;
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
    % increase, decrease, or no change between the two conditions. Only analyzes
    % clusters that are present in both conditions.
    %
    % Inputs:
    %   data1 - Structure containing group 1 data for each cluster
    %   data2 - Structure containing group 2 data for each cluster
    %   global_clusters - Array of cluster IDs to analyze (common clusters only)
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
        
        % Get data for this cluster (both should have data since we only analyze common clusters)
        data_group1 = data1.(cluster_str);
        data_group2 = data2.(cluster_str);
        
        % Calculate means and statistics
        mean_group1 = nanmean(data_group1);
        std_group1 = nanstd(data_group1);
        n_group1 = length(data_group1);
        
        mean_group2 = nanmean(data_group2);
        std_group2 = nanstd(data_group2);
        n_group2 = length(data_group2);
        
        % Calculate difference (group2 - group1)
        difference = mean_group2 - mean_group1;
        
        % Perform statistical test (t-test since both groups have data)
        p_value = 1; % Default to non-significant
        tstat = 0;
        
        try
            [~, p_value, ~, stats] = ttest2(data_group1, data_group2);
            tstat = stats.tstat;
        catch
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
