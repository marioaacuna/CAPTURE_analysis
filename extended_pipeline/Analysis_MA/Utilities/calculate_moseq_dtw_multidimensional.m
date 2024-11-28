% analysis dtw in mocap data
function [dtw_distances_F, dtw_distances_S] = calculate_moseq_dtw_multidimensional(moseq_data, cond_inds, conditions)
    % Initialize empty arrays to store the DTW distances
    dtw_distances_F = [];
    dtw_distances_S = [];

    % Get the unique list of animal IDs
    animal_ids = unique(cond_inds, 'stable');

    body_parts = fieldnames(moseq_data);

    % Loop over each animal ID to get the corresponding behavioral sequence
    for i = 1:length(animal_ids)
        for j = i+1:length(animal_ids)
            animal_idx_i = animal_ids(i);
            animal_idx_j = animal_ids(j);

            % Extract condition for each animal
            cond_i = conditions{animal_idx_i};
            cond_j = conditions{animal_idx_j};

            % Calculate DTW only if the animals are in the same condition
            if cond_i == cond_j
                dtw_distance = 0;
                for part = 1:length(body_parts)
                    % Extract 3D coordinates for the current body part
                    coords_i = moseq_data.(body_parts{part})(cond_inds == animal_idx_i, :);
                    coords_j = moseq_data.(body_parts{part})(cond_inds == animal_idx_j, :);

                    % Apply multi-dimensional DTW
                    dist = dtw_multidimensional(coords_i, coords_j);
                    dtw_distance = dtw_distance + dist;
                end

                % Store the calculated DTW distance based on the condition
                if cond_i == 'F'
                    dtw_distances_F = [dtw_distances_F, dtw_distance];
                elseif cond_i == 'S'
                    dtw_distances_S = [dtw_distances_S, dtw_distance];
                end
            end
        end
    end
end




% Multi-dimensional DTW function
function dist = dtw_multidimensional(ts1, ts2)
    % Assuming ts1 and ts2 are n-by-3 matrices (n points in 3D space)
    n = size(ts1, 1);
    m = size(ts2, 1);
    DTW = zeros(n, m);
    
    % Initialization
    DTW(1, 1) = norm(ts1(1, :) - ts2(1, :));
    for i = 2:n
        DTW(i, 1) = norm(ts1(i, :) - ts2(1, :)) + DTW(i-1, 1);
    end
    for j = 2:m
        DTW(1, j) = norm(ts1(1, :) - ts2(j, :)) + DTW(1, j-1);
    end
    
    % Main loop
    for i = 2:n
        for j = 2:m
            cost = norm(ts1(i, :) - ts2(j, :));
            DTW(i, j) = cost + min([DTW(i-1, j), DTW(i, j-1), DTW(i-1, j-1)]);
        end
    end
    dist = DTW(n, m);
end
