function sparsity_value = calculate_sparsity(transition_matrix)
    % Get the number of clusters (or behaviors)
    num_clusters = size(transition_matrix, 1);
    
    % Initialize sparsity_value
    sparsity_value = 0;
    
    % Calculate the frequency of each behavior (row sums of the transition matrix)
    behavior_frequency = sum(transition_matrix, 2);
    
    % Normalize to create probabilities
    behavior_frequency = behavior_frequency / sum(behavior_frequency);
    
    % Normalize the transition matrix to get probabilities
    transition_prob_matrix = bsxfun(@rdivide, transition_matrix, sum(transition_matrix, 2));
    
    % Calculate the sparsity
    for i = 1:num_clusters
        for j = 1:num_clusters
            sparsity_value = sparsity_value + (behavior_frequency(i) * (transition_prob_matrix(i, j)^2));
        end
    end
end
