 
function [ensemble_activity] = analyze_neural_ensembles_poses(traces, cluster_vector_ds)
unique_clusters = unique(cluster_vector_ds);
num_neurons = size(traces, 1);
num_components = 3; 
% Calculate covariance matrix for each pose cluster
ensemble_activity = struct();

for i = 1:length(unique_clusters)
    cluster_mask = cluster_vector_ds == unique_clusters(i);
    if sum(cluster_mask) <= 20%3
        % Store ensemble information
        ensemble_activity(i).principal_components = []; % Top 3 components
        ensemble_activity(i).explained_variance = [];
        ensemble_activity(i).neural_weights = []; % Weight of each neuron
        ensemble_activity(i).pca_score = [];

        continue
    end
    cluster_traces = traces(:, cluster_mask);

    % PCA to identify ensembles
    [coeff, score, latent] = pca(cluster_traces', 'NumComponents',num_components);

    % Store ensemble information
    ensemble_activity(i).principal_components = coeff(:,1:end); % Top 3 components
    ensemble_activity(i).explained_variance = latent;
    ensemble_activity(i).neural_weights = abs(coeff(:,1)); % Weight of each neuron
    ensemble_activity(i).pca_score = score;
end

end