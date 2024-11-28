function transition_matrix = calculate_transition_matrix(cluster_sequence, unique_clusters)
    % Create a mapping from cluster ID to index
    cluster_to_index_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for i = 1:length(unique_clusters)
        cluster_to_index_map(unique_clusters(i)) = i;
    end
    
    % Initialize transition matrix with zeros
    num_clusters = length(unique_clusters);
    transition_matrix = zeros(num_clusters, num_clusters);
    
    % Count the transitions from each cluster to each other cluster
    for i = 1:length(cluster_sequence) - 1
        current_cluster = cluster_sequence(i);
        next_cluster = cluster_sequence(i + 1);
        
        % Map the cluster IDs to indices
        current_index = cluster_to_index_map(current_cluster);
        next_index = cluster_to_index_map(next_cluster);
        
        % Update the transition matrix
        transition_matrix(current_index, next_index) = transition_matrix(current_index, next_index) + 1;
    end
    
    % Normalize the transition matrix rows to sum to 1 (convert counts to probabilities)
    row_sums = sum(transition_matrix, 2);
    transition_matrix = bsxfun(@rdivide, transition_matrix, row_sums);
    
    % Replace NaNs with 0 (this can happen if a row sums to zero)
    transition_matrix(isnan(transition_matrix)) = 0;
end

% 
% 
% 
% % 
% % function transition_matrix = calculate_transition_matrix(sequence, num_states)
% %     % Initialize the transition matrix with zeros
% %     transition_matrix = zeros(num_states, num_states);
% % 
% %     % Count the number of transitions between each state
% %     for i = 1:length(sequence) - 1
% %         current_state = sequence(i);
% %         next_state = sequence(i + 1);
% % 
% %         transition_matrix(current_state, next_state) = transition_matrix(current_state, next_state) + 1;
% %     end
% % 
% function transition_matrix = calculate_transition_matrix(cluster_sequence)
%     % Get unique clusters to determine the size of the transition matrix
%     unique_clusters = unique(cluster_sequence);
%     num_clusters = length(unique_clusters);
% 
%     % Initialize transition matrix with zeros
%     transition_matrix = zeros(num_clusters, num_clusters);
% 
%     % Count the transitions from each cluster to each other cluster
%     for i = 1:length(cluster_sequence) - 1
%         current_cluster = cluster_sequence(i);
%         next_cluster = cluster_sequence(i + 1);
% 
%         % Update the transition matrix
%         transition_matrix(current_cluster, next_cluster) = transition_matrix(current_cluster, next_cluster) + 1;
%     end
% 
%     % Normalize the transition matrix rows to sum to 1 (convert counts to probabilities)
%     row_sums = sum(transition_matrix, 2);
%     transition_matrix = bsxfun(@rdivide, transition_matrix, row_sums);
% 
%     % Replace NaNs with 0 (this can happen if a row sums to zero)
%     transition_matrix(isnan(transition_matrix)) = 0;
% end
