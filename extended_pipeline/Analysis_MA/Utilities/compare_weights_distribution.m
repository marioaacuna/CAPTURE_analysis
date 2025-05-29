
% For statistical comparison between conditions:
function compare_weights_distribution(group_analysis)
    % Calculate mean weights per neuron for each condition
    mean_weights_H = nanmean(group_analysis.H.weights, 1);
    mean_weights_N = nanmean(group_analysis.N.weights, 1);
    
    % Perform statistical test
    [h, p] = ttest2(group_analysis.H.weights, group_analysis.N.weights);
    
    % Plot comparison
    figure('Position', [100 100 800 400]);
    
    % Plot mean weights comparison
    subplot(1,2,1)
    scatter(mean_weights_H, mean_weights_N, '.')
    hold on
    plot([min([mean_weights_H mean_weights_N]) max([mean_weights_H mean_weights_N])], ...
         [min([mean_weights_H mean_weights_N]) max([mean_weights_H mean_weights_N])], 'k--')
    xlabel('Mean Weights - Sham')
    ylabel('Mean Weights - Neuropathic')
    title(sprintf('Weight Comparison\np = %.3f', p))
    axis square
    
    % Plot weight distributions
    subplot(1,2,2)
    boxplot([group_analysis.H.weights(:); group_analysis.N.weights(:)], ...
               [ones(size(group_analysis.H.weights(:))); 2*ones(size(group_analysis.N.weights(:)))]);
    xticklabels({'Sham', 'Neuropathic'})
    ylabel('Neural Weights')
    title('Weight Distributions')
end