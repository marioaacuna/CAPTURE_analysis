function visualize_group_heatmap(group_analysis)
    % Create figure
    figure('Position', [100 100 1200 800]);
    
    % Plot heatmap for Sham condition
    subplot(2,1,1)
    imagesc(group_analysis.H.weights)
    colorbar
    title('Neural Weights Distribution - Sham')
    xlabel('Neuron ID')
    ylabel('Sample')
    colormap('jet')  % You can change the colormap (e.g., 'parula', 'viridis', etc.)
    
    % Plot heatmap for Neuropathic condition
    subplot(2,1,2)
    imagesc(group_analysis.N.weights)
    colorbar
    title('Neural Weights Distribution - Neuropathic')
    xlabel('Neuron ID')
    ylabel('Sample')
    colormap('jet')
    
    % Adjust colorbar limits to be the same for both plots
    combined_weights = [group_analysis.H.weights; group_analysis.N.weights];
    clim([min(combined_weights(:)) max(combined_weights(:))])
end
