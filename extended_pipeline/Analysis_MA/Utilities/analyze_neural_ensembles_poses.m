 
function [ensemble_activity] = analyze_neural_ensembles_poses(traces, cluster_vector_ds)
unique_clusters = unique(cluster_vector_ds);
num_neurons = size(traces, 1);
num_components = 3; 
% Calculate covariance matrix for each pose cluster
ensemble_activity = struct();

for i = 1:length(unique_clusters)
    cluster_mask = cluster_vector_ds == unique_clusters(i);
    if sum(cluster_mask) <= 25%3
        % Store ensemble information
        ensemble_activity(i).principal_components = []; % Top 3 components
        ensemble_activity(i).explained_variance = [];
        ensemble_activity(i).neural_weights = []; % Weight of each neuron
        ensemble_activity(i).pca_score = [];
        ensemble_activity(i).peaks = [];
        ensemble_activity(i).freqs = [];

        continue
    end
    cluster_traces = traces(:, cluster_mask);

    % PCA to identify ensembles
    [coeff, score, latent, ~, ex] = pca(cluster_traces', 'NumComponents',2);
    
    % Calculate the mean amplitude of events
    % iterate through cells
    peaks = [];
    freqs= [];
    for icell = 1:size(cluster_traces,1)
        this_trace = cluster_traces(icell,:);
        smoothed = smooth(this_trace);

        thr = std(smoothed) * 1.5;
        l = findpeaks(smoothed,thr);% This is coming from the Utility script
        pks = this_trace(l.loc);
        peaks(icell,1)= mean(pks);

        freq = length(l.loc) / (size(this_trace,2)/5);
        freqs(icell,1) = freq;

    end
    % Store ensemble information
    ensemble_activity(i).principal_components = coeff(:,1:end); % Top 3 components
    ensemble_activity(i).explained_variance = latent;
    ensemble_activity(i).neural_weights = abs(coeff(:,1)); % Weight of each neuron
    ensemble_activity(i).pca_score = score;

    % Traces parameters
    ensemble_activity(i).peaks = peaks;
    ensemble_activity(i).freqs = freqs;

    ensemble_activity(i).traces = cluster_traces;

end

end