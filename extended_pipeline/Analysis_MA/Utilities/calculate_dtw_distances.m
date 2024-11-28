function [dtw_distances_F, dtw_distances_S] = calculate_dtw_distances(behavior_sequences, cond_inds, conditions)
    % Get the unique list of animal IDs
    animal_ids = unique(cond_inds, 'stable');
    
    % Initialize empty arrays to store the DTW distances
    dtw_distances_F = [];
    dtw_distances_S = [];
    
    % Loop over each animal ID to get the corresponding behavioral sequence
    for i = 1:length(animal_ids)
        for j = i+1:length(animal_ids)
            animal_idx_i = animal_ids(i);
            animal_idx_j = animal_ids(j);
            
            % Extract sequences for each animal
            seq_i = behavior_sequences(cond_inds == animal_idx_i);
            seq_j = behavior_sequences(cond_inds == animal_idx_j);
            
            % Extract condition for each animal
            cond_i = conditions{animal_idx_i};
            cond_j = conditions{animal_idx_j};
            
            % Calculate DTW only if the animals are in the same condition
            if cond_i == 'F' && cond_j == 'F'
                dtw_distances_F = [dtw_distances_F, dtw(seq_i, seq_j)];
            elseif cond_i == 'S' && cond_j == 'S'
                dtw_distances_S = [dtw_distances_S, dtw(seq_i, seq_j)];
            end
        end
    end
end
