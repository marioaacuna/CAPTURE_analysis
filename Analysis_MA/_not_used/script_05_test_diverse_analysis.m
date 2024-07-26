% Initialize containers for storing transition matrices and sparsity values
% this needs to be run after script 3.
transition_matrices_F = {};  % For Formalin
transition_matrices_S = {};  % For Saline
sparcity_values_F = [];  % For Formalin
sparcity_values_S = [];  % For Saline
entropy_values_F = [];  % For Formalin
entropy_values_S = [];  % For Saline
auto_corr_F = [];
auto_corr_S = []; 

clustered_behavior = 1 +  hierarchystruct.clustered_behavior{1};
cond_ids_unique = unique(cond_inds, 'stable');
%% Loop over each animal
for i = 1:length(cond_ids_unique)
    
    % Get the animal index and corresponding condition
    animal_idx = cond_ids_unique(i);
    animal_condition = conditions{animal_idx};
    
    % Skip if the condition is 'N' (Naive)
    if strcmp(animal_condition, 'N')
        continue;
    end
    
    % Extract the behavioral sequence for this animal
    behavior_sequence = clustered_behavior(cond_inds == animal_idx);
    
    % Calculate transition matrix for this animal
    unique_clusters = unique(behavior_sequence);
    transition_matrix = calculate_transition_matrix(behavior_sequence, unique_clusters); 
    
    % Calculate sparsity for this animal
    sparcity_value = calculate_sparsity(transition_matrix);  
    
    % Calculate entropy
    probability_distribution = histcounts(behavior_sequence, 'Normalization', 'probability');
    entropy_value = -sum(probability_distribution .* log2(probability_distribution + eps));
    
    % calculate Autocorrelation
    % auto_corr = xcorr(behavior_sequence - mean(behavior_sequence), 'coeff');

    % Store results based on the condition
    if strcmp(animal_condition, 'F')
        transition_matrices_F{end+1} = transition_matrix;
        sparcity_values_F(end+1) = sparcity_value;
        entropy_values_F(end+1) = entropy_value;
        % auto_corr_F(end+1,:) = auto_corr;
    elseif strcmp(animal_condition, 'S')
        transition_matrices_S{end+1} = transition_matrix;
        sparcity_values_S(:, end+1) = sparcity_value;
        entropy_values_S(:, end+1) = entropy_value;
        % auto_corr_S(end+1,:) = auto_corr;
    end
    
end


% Perform group-level analysis, such as averaging or statistical tests,
% on transition_matrices_F, transition_matrices_S, sparcity_values_F, and sparcity_values_S

%% 
% Simulated data: Replace these with your actual sequences and condition indices
behavior_sequences = clustered_behavior;
condition_indices = conditions;

[dtw_distances_F, dtw_distances_S] = calculate_dtw_distances(behavior_sequences, cond_inds, conditions);

%% apply multidimensional dtw for the moseq dat
moseq_data = analysisstruct.mocapstruct_reduced_agg{1}.markers_preproc;
[dtw_distances_F, dtw_distances_S] = calculate_moseq_dtw_multidimensional(moseq_data, cond_inds, conditions);
