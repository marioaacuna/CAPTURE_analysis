    % Initialize arrays with proper structure using global_clusters
    num_global_clusters = length(global_clusters);
    H_weights_by_cluster = NaN(num_global_clusters, 0);  % Clusters x concatenated_neurons
    N_weights_by_cluster = NaN(num_global_clusters, 0);  % Clusters x concatenated_neurons
    
    neuron_count_H = 0;
    neuron_count_N = 0;
    
    % Get animal IDs
    H_animals = fieldnames(data_H);
    N_animals = fieldnames(data_N);
    
    % Process Sham animals
    for h = 1:length(H_animals)
        animal_ID = H_animals{h};
        animal_data = data_H.(animal_ID).ensemble_activity;
        
        % Get number of neurons for this animal
        num_neurons = size(animal_data(11).neural_weights, 1);
        
        % Initialize temporary matrix
        temp_weights = NaN(num_global_clusters, num_neurons);
        
        % Fill matrix with weights matching global clusters
        for local_cluster = 1:length(animal_data)
            if ~isempty(animal_data(local_cluster).neural_weights)
                % Find where this local cluster maps to in global_clusters
                global_idx = find(global_clusters == data_H.(animal_ID).unique_clusters(local_cluster));
                if ~isempty(global_idx)
                    temp_weights(global_idx, :) = animal_data(local_cluster).neural_weights';
                end
            end
        end
        
        % Concatenate
        H_weights_by_cluster = [H_weights_by_cluster, temp_weights];
        neuron_count_H = neuron_count_H + num_neurons;
        fprintf('Added %d neurons from Sham animal %s\n', num_neurons, animal_ID);
    end
    
    % Process Neuropathic animals
    for n = 1:length(N_animals)
        animal_ID = N_animals{n};
        animal_data = data_N.(animal_ID).ensemble_activity;
        
        % Get number of neurons for this animal
        num_neurons = size(animal_data(5).neural_weights, 1);
        
        % Initialize temporary matrix
        temp_weights = NaN(num_global_clusters, num_neurons);
        
        % Fill matrix with weights matching global clusters
        for local_cluster = 1:length(animal_data)
            if ~isempty(animal_data(local_cluster).neural_weights)
                % Find where this local cluster maps to in global_clusters
                global_idx = find(global_clusters == data_N.(animal_ID).unique_clusters(local_cluster));
                if ~isempty(global_idx)
                    temp_weights(global_idx, :) = animal_data(local_cluster).neural_weights';
                end
            end
        end
        
        % Concatenate
        N_weights_by_cluster = [N_weights_by_cluster, temp_weights];
        neuron_count_N = neuron_count_N + num_neurons;
        fprintf('Added %d neurons from Neuropathic animal %s\n', num_neurons, animal_ID);
    end
    
    %% Plot
    % Create figure
    figure('Position', [100 100 1200 600]);
    
    % Sort by cluster IDs
    [sorted_clusters, sort_idx] = sort(global_clusters);
    
    % Apply sorting to the weight matrices
    H_weights_sorted = H_weights_by_cluster(sort_idx, :);
    N_weights_sorted = N_weights_by_cluster(sort_idx, :);

    % Sort neurons by their mean weight across clusters
    mean_weights_H = nanmean(H_weights_sorted, 1);
    mean_weights_N = nanmean(N_weights_sorted, 1);
    [~, neuron_order_H] = sort(mean_weights_H, 'descend');
    [~, neuron_order_N] = sort(mean_weights_N, 'descend');

    % Reorder the matrices
    H_weights_sorted = H_weights_sorted(:, neuron_order_H);
    N_weights_sorted = N_weights_sorted(:, neuron_order_N);
    
    
    % Plot Sham condition
    subplot(1,2,1)
    imagesc(H_weights_sorted)
    colorbar
    title('Neural Ensembles Across Clusters - Sham')
    xlabel(['Neuron ID (Total: ' num2str(neuron_count_H) ')'])
    ylabel('Global Cluster ID')
    
    % Set y-axis ticks to show actual cluster IDs
    if length(global_clusters) > 20
        yticks(1:5:length(global_clusters))
        yticklabels(sorted_clusters(1:5:end))
    else
        yticks(1:length(global_clusters))
        yticklabels(sorted_clusters)
    end
    colormap('jet')
    
    % Plot Neuropathic condition
    subplot(1,2,2)
    imagesc(N_weights_sorted)
    colorbar
    title('Neural Ensembles Across Clusters - Neuropathic')
    xlabel(['Neuron ID (Total: ' num2str(neuron_count_N) ')'])
    ylabel('Global Cluster ID')
    
    % Set y-axis ticks to show actual cluster IDs
    if length(global_clusters) > 20
        yticks(1:5:length(global_clusters))
        yticklabels(sorted_clusters(1:5:end))
    else
        yticks(1:length(global_clusters))
        yticklabels(sorted_clusters)
    end
    colormap('jet')
%% Plot correlation and do analysis
    % 1. Correlation Analysis
    H_corr = corr(H_weights_sorted');
    N_corr = corr(N_weights_sorted');
    
    figure('Position', [100 100 1500 800]);
    
    % Plot correlations
    subplot(2,3,1)
    imagesc(H_corr)
    title('Cluster Correlations - Sham')
    colorbar
    xlabel('Sorted Cluster ID')
    ylabel('Sorted Cluster ID')
    % Add actual cluster IDs if not too many
    if length(sorted_clusters) <= 20
        xticks(1:length(sorted_clusters))
        yticks(1:length(sorted_clusters))
        xticklabels(sorted_clusters)
        yticklabels(sorted_clusters)
    end
    
    subplot(2,3,2)
    imagesc(N_corr)
    title('Cluster Correlations - Neuropathic')
    colorbar
    xlabel('Sorted Cluster ID')
    ylabel('Sorted Cluster ID')
    if length(sorted_clusters) <= 20
        xticks(1:length(sorted_clusters))
        yticks(1:length(sorted_clusters))
        xticklabels(sorted_clusters)
        yticklabels(sorted_clusters)
    end
    
    subplot(2,3,3)
    imagesc(H_corr - N_corr)
    title('Correlation Difference (Sham - Neuropathic)')
    colorbar
    xlabel('Sorted Cluster ID')
    ylabel('Sorted Cluster ID')
    if length(sorted_clusters) <= 20
        xticks(1:length(sorted_clusters))
        yticks(1:length(sorted_clusters))
        xticklabels(sorted_clusters)
        yticklabels(sorted_clusters)
    end
    
    % 2. Important Neurons Analysis
    % Calculate mean weights per neuron
    H_mean_weights = mean(H_weights_sorted, 1);
    N_mean_weights = mean(N_weights_sorted, 1);
    
    % Get top neurons (e.g., top 10%)
    threshold_percentile = 90;
    H_threshold = prctile(H_mean_weights, threshold_percentile);
    N_threshold = prctile(N_mean_weights, threshold_percentile);
    
    H_important = H_mean_weights > H_threshold;
    N_important = N_mean_weights > N_threshold;
    
    % Plot weight distributions
    subplot(2,3,4)
    histogram(H_mean_weights, 50, 'Normalization', 'probability', 'FaceColor', 'b', 'FaceAlpha', 0.5)
    hold on
    histogram(N_mean_weights, 50, 'Normalization', 'probability', 'FaceColor', 'r', 'FaceAlpha', 0.5)
    xline(H_threshold, 'b--')
    xline(N_threshold, 'r--')
    title('Distribution of Mean Neuronal Weights')
    legend('Sham', 'Neuropathic')
    xlabel('Mean Weight')
    ylabel('Probability')
    
    % 3. Cluster Activity Pattern Analysis
    % Calculate variance of weights per cluster
    H_cluster_var = var(H_weights_sorted, [], 2);
    N_cluster_var = var(N_weights_sorted, [], 2);
    
    subplot(2,3,5)
    plot(sorted_clusters, H_cluster_var, 'b', 'LineWidth', 2)
    hold on
    plot(sorted_clusters, N_cluster_var, 'r', 'LineWidth', 2)
    title('Weight Variance per Cluster')
    xlabel('Cluster ID')
    ylabel('Variance')
    legend('Sham', 'Neuropathic')
    
    % 4. Statistical Comparisons
    subplot(2,3,6)
    % Compare overall weight distributions
    [h,p] = kstest2(H_mean_weights, N_mean_weights);
    % 
    % boxplot([H_mean_weights'; N_mean_weights'], ...
    %         {'Sham', 'Neuropathic'})
    % title(sprintf('Weight Distributions\nKS-test p = %.3f', p))
    % ylabel('Mean Weight')
    
    % Print summary statistics
    fprintf('\nAnalysis Summary:\n')
    fprintf('Sham condition:\n')
    fprintf('  Number of important neurons: %d (%.1f%%)\n', ...
            sum(H_important), 100*sum(H_important)/length(H_important))
    fprintf('  Mean weight: %.3f ± %.3f\n', ...
            mean(H_mean_weights), std(H_mean_weights))
    
    fprintf('\nNeuropathic condition:\n')
    fprintf('  Number of important neurons: %d (%.1f%%)\n', ...
            sum(N_important), 100*sum(N_important)/length(N_important))
    fprintf('  Mean weight: %.3f ± %.3f\n', ...
            mean(N_mean_weights), std(N_mean_weights))
    
    % Additional cluster-specific analysis
    high_var_clusters_H = find(H_cluster_var > prctile(H_cluster_var, 90));
    high_var_clusters_N = find(N_cluster_var > prctile(N_cluster_var, 90));
    
    fprintf('\nClusters with high variability:\n')
    fprintf('Sham clusters: %s\n', num2str(sorted_clusters(high_var_clusters_H)'))
    fprintf('Neuropathic clusters: %s\n', num2str(sorted_clusters(high_var_clusters_N)'))


    %% distribution
    % Find active clusters (non-zero rows) in both conditions
    active_H = nansum(H_weights_sorted, 2) > 0;  % True for rows with any non-zero values
    active_N = nansum(N_weights_sorted, 2) > 0;
    
    % Find common active clusters
    common_clusters = active_H & active_N;
    common_cluster_ids = sorted_clusters(common_clusters);
    
    fprintf('Analysis of common clusters:\n')
    fprintf('Total clusters: %d\n', length(sorted_clusters))
    fprintf('Active in Sham: %d\n', sum(active_H))
    fprintf('Active in Neuropathic: %d\n', sum(active_N))
    fprintf('Common active clusters: %d\n', sum(common_clusters))
    
    % Continue analysis only with common clusters
    H_weights_common = H_weights_sorted(common_clusters, :);
    N_weights_common = N_weights_sorted(common_clusters, :);
    
    % Figure for distributions
    figure('Position', [100 100 1200 800]);
    
    % For each common cluster
    num_common = nansum(common_clusters);
    p_values = zeros(num_common, 1);
    
    % Store KS-test statistics
    for i = 1:num_common
        % Get weights for this cluster
        H_cluster_weights = H_weights_common(i, :);
        N_cluster_weights = N_weights_common(i, :);
        
        % Remove zeros (non-contributing neurons)
        H_cluster_weights = H_cluster_weights(H_cluster_weights > 0);
        N_cluster_weights = N_cluster_weights(N_cluster_weights > 0);
        
        % Perform KS test
        % [~, p_values(i)] = kstest2(H_cluster_weights, N_cluster_weights);
        [~, p_values(i)] = ttest2(H_cluster_weights, N_cluster_weights);
    end
    
    % Find significantly different clusters
    sig_clusters = find(p_values < 0.05);
    
    % Plot top different clusters
    [~, most_diff_idx] = sort(p_values);
    num_to_plot = min(6, length(sig_clusters));
    
    for plot_idx = 1:num_to_plot
        cluster_idx = most_diff_idx(plot_idx);
        
        subplot(2,ceil(min(6, length(sig_clusters))/2),plot_idx)
        % Get weights for this cluster
        H_cluster_weights = H_weights_common(cluster_idx, :);
        N_cluster_weights = N_weights_common(cluster_idx, :);
        
        % Remove zeros
        H_cluster_weights = H_cluster_weights(H_cluster_weights > 0);
        N_cluster_weights = N_cluster_weights(N_cluster_weights > 0);
        
        % Plot distributions
        histogram(H_cluster_weights, 50, 'Normalization', 'probability', ...
                 'FaceColor', 'b', 'FaceAlpha', 0.5)
        hold on
        histogram(N_cluster_weights, 50, 'Normalization', 'probability', ...
                 'FaceColor', 'r', 'FaceAlpha', 0.5)
        % Plot means as vertical lines
        mean_H = mean(H_cluster_weights);
        mean_N = mean(N_cluster_weights);
        
        xline(mean_H, 'b-', 'LineWidth', 2);
        xline(mean_N, 'r-', 'LineWidth', 2);
        
        
        title(sprintf('Cluster %d (p=%.3f)', common_cluster_ids(cluster_idx), p_values(cluster_idx)))
        xlabel('Weight Value')
        ylabel('Probability')
        legend('Sham', 'Neuropathic')
    end
    
    % Create summary figure
    figure('Position', [100 100 800 400]);
    
    % Plot p-values for common clusters
    subplot(1,2,1)
    plot(common_cluster_ids, -log10(p_values), 'k.', 'MarkerSize', 15)
    hold on
    yline(-log10(0.05), 'r--')
    xlabel('Cluster ID')
    ylabel('-log10(p-value)')
    title('Distribution Difference Significance')
    
    % Plot number of contributing neurons per cluster
    subplot(1,2,2)
    num_neurons_H = sum(H_weights_common > 0, 2);
    num_neurons_N = sum(N_weights_common > 0, 2);
    
    plot(common_cluster_ids, num_neurons_H, 'b.', 'MarkerSize', 15)
    hold on
    plot(common_cluster_ids, num_neurons_N, 'r.', 'MarkerSize', 15)
    xlabel('Cluster ID')
    ylabel('Number of Contributing Neurons')
    legend('Sham', 'Neuropathic')
    title('Neuronal Participation per Common Cluster')
    
    % Print summary
    fprintf('\nDistribution Analysis Summary (Common Clusters Only):\n')
    fprintf('Number of significantly different clusters: %d (%.1f%%)\n', ...
            length(sig_clusters), 100*length(sig_clusters)/num_common)
    
    fprintf('\nTop 5 most different clusters:\n')
    for i = 1:min(5, length(most_diff_idx))
        cluster_idx = most_diff_idx(i);
        fprintf('Cluster %d: p = %.3e\n', common_cluster_ids(cluster_idx), p_values(cluster_idx))
    end
