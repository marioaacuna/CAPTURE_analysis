% This test runs after the hierarchystruct was generated.

length_of_your_data = length(hierarchystruct.clustered_behavior{2});
% get eigenvectors
[V, D] = eig(transition_matrix');

% Find Leading Eigenvalue and Eigenvector:
[max_eigenvalue, idx] = max(diag(D));
leading_eigenvector = V(:, idx);

% Check Perron-Frobenius Conditions:
if abs(max_eigenvalue - 1) < 1e-6
    fprintf('Leading eigenvalue is close to 1, data could be Markovian.\n');
else
    fprintf('Leading eigenvalue deviates from 1, data is likely not Markovian.\n');
end

% Simulate Markov Chain for Comparison: 
% Simulate Markov Chain
simulated_data = simulate_markov_chain(transition_matrix, length_of_your_data);  % You'll have to write this function

% Calculate transition matrix for simulated data
simulated_transition_matrix = calculate_transition_matrix(simulated_data, size(transition_matrix,1));  % You'll have to write this function

% Eigenvalue decomposition for simulated data
[V_sim, D_sim] = eig(simulated_transition_matrix');

% compaire Eigenvalues
fprintf('Largest Eigenvalue for observed data: %f\n', max_eigenvalue);
fprintf('Largest Eigenvalue for simulated data: %f\n', max(diag(D_sim)));

%% Calculate sparcity
% Assume transition_matrix is your transition matrix
sparsity_value = calculate_sparsity(transition_matrix);
fprintf('Sparsity val: %f \n',sparsity_value)

%% HELPER functions
function simulated_chain = simulate_markov_chain(transition_matrix, chain_length)
    % Initialize the simulated chain
    simulated_chain = zeros(1, chain_length);
    
    % Number of states
    num_states = size(transition_matrix, 1);
    
    % Start with a random state
    simulated_chain(1) = randi(num_states);
    
    % Loop to generate the Markov chain
    for t = 2:chain_length
        current_state = simulated_chain(t - 1);
        
        % Transition probabilities for the current state
        transition_prob = transition_matrix(current_state, :);
        
        % Generate the next state based on these probabilities
        next_state = randsample(1:num_states, 1, true, transition_prob);
        
        % Assign the next state to the chain
        simulated_chain(t) = next_state;
    end
end

% function transition_matrix = calculate_transition_matrix(sequence, num_states)
%     % Initialize the transition matrix with zeros
%     transition_matrix = zeros(num_states, num_states);
% 
%     % Count the number of transitions between each state
%     for i = 1:length(sequence) - 1
%         current_state = sequence(i);
%         next_state = sequence(i + 1);
% 
%         transition_matrix(current_state, next_state) = transition_matrix(current_state, next_state) + 1;
%     end
% 
%     % Normalize each row to get probabilities
%     for i = 1:num_states
%         row_sum = sum(transition_matrix(i, :));
%         if row_sum > 0  % Avoid division by zero
%             transition_matrix(i, :) = transition_matrix(i, :) / row_sum;
%         end
%     end
% end
% 
% 
% function sparsity_value = calculate_sparsity(transition_matrix)
%     % Get the number of clusters (or behaviors)
%     num_clusters = size(transition_matrix, 1);
% 
%     % Initialize sparsity_value
%     sparsity_value = 0;
% 
%     % Calculate the frequency of each behavior (row sums of the transition matrix)
%     behavior_frequency = sum(transition_matrix, 2);
% 
%     % Normalize to create probabilities
%     behavior_frequency = behavior_frequency / sum(behavior_frequency);
% 
%     % Normalize the transition matrix to get probabilities
%     transition_prob_matrix = bsxfun(@rdivide, transition_matrix, sum(transition_matrix, 2));
% 
%     % Calculate the sparsity
%     for i = 1:num_clusters
%         for j = 1:num_clusters
%             sparsity_value = sparsity_value + (behavior_frequency(i) * (transition_prob_matrix(i, j)^2));
%         end
%     end
% end
