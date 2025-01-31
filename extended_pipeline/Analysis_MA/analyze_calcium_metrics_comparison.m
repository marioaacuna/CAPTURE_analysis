function analyze_calcium_metrics_comparison(data1, data2, animals_of_interest, metric_to_take, comparison)
    % analyze_calcium_metrics_comparison - Analyzes and visualizes calcium imaging metrics
    % 
    % Inputs:
    %   data1 - Structure containing control group data (Saline/Sham)
    %   data2 - Structure containing experimental group data (Formalin/Neuropathic)
    %   animals_of_interest - Cell array of animal IDs
    %   metric_to_take - String specifying metric ('max_amplitude', 'peaks', 'freqs')
    %   comparison - String specifying comparison type ('S_v_F' or 'H_v_N')
    
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
    end

    % Initialize structures for all ROIs per cluster
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
        all_ROIs_per_cluster_1.(cluster_str) = all_ROIs_1;
        all_ROIs_per_cluster_2.(cluster_str) = all_ROIs_2;
    end

    % Generate plots
    plot_bar_comparison(all_ROIs_per_cluster_1, all_ROIs_per_cluster_2, group1_label, group2_label, metric_to_take);
    plot_difference_chart(all_ROIs_per_cluster_1, all_ROIs_per_cluster_2, group1_label, group2_label, metric_to_take);
    plot_heatmap(all_ROIs_per_cluster_1, all_ROIs_per_cluster_2, group1_label, group2_label, metric_to_take);
    plot_scatter_comparison(all_ROIs_per_cluster_1, all_ROIs_per_cluster_2, group1_label, group2_label, metric_to_take);
end

function plot_bar_comparison(data1, data2, group1_label, group2_label, metric_to_take)
    Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]);
    
    % Define plotting constants
    barWidth = 0.75;
    gapWidth = 1;
    currentX = 1;
    
    clusters = fieldnames(data1);
    
    for i = 1:length(clusters)
        cluster_name = clusters{i};
        
        data_group1 = data1.(cluster_name);
        data_group2 = data2.(cluster_name);
        
        % Calculate statistics
        try
            [h, p, ci, stats] = ttest2(data_group1, data_group2);
        catch
            p = 0;
            stats.tstat = [];
        end
        
        % Calculate means and SEMs
        mean_1 = nanmean(data_group1);
        SEM_1 = std(data_group1)/sqrt(length(data_group1));
        mean_2 = nanmean(data_group2);
        SEM_2 = std(data_group2)/sqrt(length(data_group2));
        
        % Plot bars and error bars
        bar(currentX, mean_1, barWidth, 'b');
        hold on;
        errorbar(currentX, mean_1, SEM_1, 'k', 'LineStyle', 'none');
        bar(currentX + barWidth, mean_2, barWidth, 'r');
        errorbar(currentX + barWidth, mean_2, SEM_2, 'k', 'LineStyle', 'none');
        
        % Add significance markers
        try
            if p < 0.05 && p > 0.01
                text(currentX, max(mean_2, mean_1) + max(mean_2, mean_1) * 0.20, '*');
            elseif p < 0.01 && p > 0.001
                text(currentX - barWidth, max(mean_2, mean_1) + max(mean_2, mean_1) * 0.20, '**');
            elseif p < 0.001 && p > 0
                text(currentX - barWidth, max(mean_2, mean_1) + max(mean_2, mean_1) * 0.20, '***');
            end
        end
        
        currentX = currentX + barWidth * 2 + gapWidth;
    end
    
    % Add legend
    hG1 = bar(NaN, NaN, 'b');
    hG2 = bar(NaN, NaN, 'r');
    legend([hG1, hG2], {group1_label, group2_label});
    
    % Customize plot
    ylabel(['Mean ' metric_to_take]);
    xlabel('Clusters');
    title(['Mean ' metric_to_take ' of ROIs per Cluster for ' group1_label ' and ' group2_label]);
    set(gca, 'XTick', 1.5:barWidth*2+gapWidth:length(clusters)*(barWidth*2+gapWidth));
    set(gca, 'XTickLabel', clusters, 'TickLabelInterpreter', 'none');
    box off;
    set(gca, 'TickDir', 'out');
end

function plot_difference_chart(data1, data2, group1_label, group2_label, metric_to_take)
    Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]);
    
    barWidth = 0.85;
    clusters = fieldnames(data1);
    
    % Preallocate arrays
    differences = zeros(1, length(clusters));
    errors = zeros(1, length(clusters));
    p_values_array = zeros(1, length(clusters));
    
    % Calculate differences and statistics
    for i = 1:length(clusters)
        cluster_name = clusters{i};
        
        data_group1 = data1.(cluster_name);
        data_group2 = data2.(cluster_name);
        
        differences(i) = nanmean(data_group2) - nanmean(data_group1);
        errors(i) = sqrt((std(data_group2)^2/length(data_group2)) + (std(data_group1)^2/length(data_group1)));
        
        try
            [~, p] = ttest2(data_group1, data_group2);
            p_values_array(i) = p;
        catch
            p_values_array(i) = 1;
        end
    end
    
    % Create color coding based on p-values
    colors = zeros(length(clusters), 3);
    colors(p_values_array < 0.05 & p_values_array >= 0.01, :) = repmat([0.8 0.4 0], sum(p_values_array < 0.05 & p_values_array >= 0.01), 1);
    colors(p_values_array < 0.01 & p_values_array >= 0.001, :) = repmat([0.9 0.2 0], sum(p_values_array < 0.01 & p_values_array >= 0.001), 1);
    colors(p_values_array < 0.001, :) = repmat([1 0 0], sum(p_values_array < 0.001), 1);
    colors(p_values_array >= 0.05, :) = repmat([0.7 0.7 0.7], sum(p_values_array >= 0.05), 1);
    
    % Plot
    b = bar(differences, barWidth);
    b.FaceColor = 'flat';
    b.CData = colors;
    
    hold on;
    errorbar(1:length(clusters), differences, errors, 'k.', 'LineStyle', 'none');
    yline(0, 'k--', 'Alpha', 0.3);
    
    % Customize plot
    ylabel(['Difference in Mean ' metric_to_take ' (' group2_label ' - ' group1_label ')']);
    xlabel('Clusters');
    title(['Difference in ' metric_to_take ' between ' group2_label ' and ' group1_label ' Conditions']);
    set(gca, 'XTick', 1:length(clusters));
    set(gca, 'XTickLabel', clusters, 'TickLabelInterpreter', 'none');
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
end

function plot_heatmap(data1, data2, group1_label, group2_label, metric_to_take)
    Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]);
    
    clusters = fieldnames(data1);
    num_clusters = length(clusters);
    
    % Calculate means and p-values
    means_matrix = zeros(2, num_clusters);
    p_values_array = zeros(1, num_clusters);
    
    for i = 1:num_clusters
        cluster_name = clusters{i};
        
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
            text(i, 1.5, '***', 'HorizontalAlignment', 'center', 'Color', 'k');
        elseif p_values_array(i) < 0.01
            text(i, 1.5, '**', 'HorizontalAlignment', 'center', 'Color', 'k');
        elseif p_values_array(i) < 0.05
            text(i, 1.5, '*', 'HorizontalAlignment', 'center', 'Color', 'k');
        end
    end
    
    % Customize plot
    ylabel('Condition');
    xlabel('Clusters');
    title(['Mean ' metric_to_take ' Heatmap of ROIs per Cluster']);
    set(gca, 'YTick', 1:2);
    set(gca, 'YTickLabel', {group1_label, group2_label});
    set(gca, 'XTick', 1:num_clusters);
    set(gca, 'XTickLabel', clusters, 'TickLabelInterpreter', 'none');
end

function plot_scatter_comparison(data1, data2, group1_label, group2_label, metric_to_take)
    Fig_clusters_amplitude = figure('color', 'w', 'Position',[100 100 1500 700]);
    
    clusters = fieldnames(data1);
    num_clusters = length(clusters);
    
    % Calculate means and p-values
    means_1 = zeros(1, num_clusters);
    means_2 = zeros(1, num_clusters);
    p_values_array = zeros(1, num_clusters);
    
    for i = 1:num_clusters
        cluster_name = clusters{i};
        
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
        line([means_1(i) means_2(i)], [i i], 'Color', colors(i,:));
    end
    
    hold on;
    scatter(means_1, 1:num_clusters, 50, 'b', 'filled');
    scatter(means_2, 1:num_clusters, 50, 'r', 'filled');
    
    % Customize plot
    ylabel('Clusters');
    xlabel(['Mean ' metric_to_take]);
    title(['Comparison of Mean ' metric_to_take ' between Conditions']);
    set(gca, 'YTick', 1:num_clusters);
    set(gca, 'YTickLabel', clusters, 'TickLabelInterpreter', 'none');
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
end
