function analyze_group_PCA(group_analysis)
    % Create main figure
    figure('Position', [100 100 1200 800]);
    
    %% 1. Compare Explained Variance between conditions
    subplot(2,2,1)
    
    % Calculate mean and std of variance explained
    mean_var_H = mean(group_analysis.H.variance, 1);
    std_var_H = std(group_analysis.H.variance, [], 1);
    mean_var_N = mean(group_analysis.N.variance, 1);
    std_var_N = std(group_analysis.N.variance, [], 1);
    
    % Convert to cumulative sum for better visualization
    cum_mean_H = cumsum(mean_var_H);
    cum_mean_N = cumsum(mean_var_N);
    cum_std_H = cumsum(std_var_H);
    cum_std_N = cumsum(std_var_N);
    
    % Plot
    x = 1:length(mean_var_H);
    shadedErrorBar(x, cum_mean_H, cum_std_H, 'b', 0.5);
    hold on
    shadedErrorBar(x, cum_mean_N, cum_std_N, 'r', 0.5);
    
    legend('Sham', 'Neuropathic')
    title('Cumulative Variance Explained')
    xlabel('Principal Component')
    ylabel('Cumulative Variance (%)')
    
    % Statistical comparison of total variance explained
    [h,p] = ttest2(sum(group_analysis.H.variance,2), sum(group_analysis.N.variance,2));
    text(0.7*length(x), 30, sprintf('p = %.3f', p))
    
    %% 2. Neural Weights Distribution
    subplot(2,2,2)
    % Create boxplot instead of violin plot
    boxplot([group_analysis.H.weights(:), group_analysis.N.weights(:)], ...
            'Labels', {'Sham', 'Neuropathic'});
    ylabel('Neural Weights')
    title('Distribution of Neural Weights')
    
    % Add statistical comparison
    [~,p_weights] = ttest2(group_analysis.H.weights(:), group_analysis.N.weights(:));
    text(1.5, max(max(group_analysis.H.weights(:)), max(group_analysis.N.weights(:))), ...
        sprintf('p = %.3f', p_weights))
    
    %% 3. Component Structure Analysis
    subplot(2,2,3)
    % Calculate correlation between top components
    num_comp = min(10, size(group_analysis.H.components,2)); % Use top 10 components or less
    corr_matrix = zeros(num_comp);
    for i = 1:num_comp
        for j = 1:num_comp
            corr_matrix(i,j) = corr(group_analysis.H.components(:,i), ...
                                  group_analysis.N.components(:,j));
        end
    end
    
    imagesc(corr_matrix)
    colorbar
    colormap('jet')
    title('Component Correlation Between Conditions')
    xlabel('Neuropathic Components')
    ylabel('Sham Components')
    axis square
    
    %% 4. Dimensionality Analysis
    subplot(2,2,4)
    % Calculate effective dimensionality using participation ratio
    eff_dim_H = zeros(size(group_analysis.H.variance, 1), 1);
    eff_dim_N = zeros(size(group_analysis.N.variance, 1), 1);
    
    for i = 1:size(group_analysis.H.variance, 1)
        eff_dim_H(i) = (sum(group_analysis.H.variance(i,:)))^2 / ...
                        sum(group_analysis.H.variance(i,:).^2);
    end
    
    for i = 1:size(group_analysis.N.variance, 1)
        eff_dim_N(i) = (sum(group_analysis.N.variance(i,:)))^2 / ...
                        sum(group_analysis.N.variance(i,:).^2);
    end
    
    % Create boxplot of dimensionality
    boxplot([eff_dim_H; eff_dim_N], ...
            [ones(size(eff_dim_H)); 2*ones(size(eff_dim_N))], ...
            'Labels', {'Sham', 'Neuropathic'});
    ylabel('Effective Dimensionality')
    title('Neural Population Dimensionality')
    
    % Add statistical comparison
    [~,p_dim] = ttest2(eff_dim_H, eff_dim_N);
    text(1.5, max(max(eff_dim_H), max(eff_dim_N)), sprintf('p = %.3f', p_dim))
    
    % Print overall summary
    fprintf('\nStatistical Summary:\n')
    fprintf('Variance Explained: p = %.3f\n', p)
    fprintf('Neural Weights: p = %.3f\n', p_weights)
    fprintf('Dimensionality: p = %.3f\n', p_dim)
end