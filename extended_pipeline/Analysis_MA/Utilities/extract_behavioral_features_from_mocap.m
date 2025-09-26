function [behavior_states, features] = extract_behavioral_features_from_mocap(mocapstruct, animal_condition_id, animal_condition_identifier, behaviors)
    
    % Set default behaviors if not provided
    if nargin < 4 || isempty(behaviors)
        behaviors = {'rearing', 'grooming', 'walking', 'quiet'};
    end
    
    % Ensure behaviors is a cell array
    if ischar(behaviors)
        behaviors = {behaviors};
    end
    
    % Find frames belonging to this specific animal/condition
    frame_indices = find(strcmp(animal_condition_identifier, animal_condition_id));
    
    if isempty(frame_indices)
        error('No frames found for animal_condition_id: %s', animal_condition_id);
    end
    
    % Extract marker data for this specific animal/condition
    % Static posture analysis (aligned/centered data)
    markers = struct();
    marker_fields = fieldnames(mocapstruct.markers_aligned_preproc);
    for i = 1:length(marker_fields)
        markers.(marker_fields{i}) = mocapstruct.markers_aligned_preproc.(marker_fields{i})(frame_indices, :);
    end
    
    % Movement analysis (non-aligned data with movement preserved)
    markers_moving = struct();
    marker_fields_moving = fieldnames(mocapstruct.markers_preproc);
    for i = 1:length(marker_fields_moving)
        markers_moving.(marker_fields_moving{i}) = mocapstruct.markers_preproc.(marker_fields_moving{i})(frame_indices, :);
    end
    
    % Pre-calculate body length for normalization using aligned frames
    body_length = mean(vecnorm(markers.SpineF - markers.Tail_base_, 2, 2));
    
    n_frames = size(markers_moving.SpineM, 1);
    
    %% HIERARCHICAL BEHAVIORAL ANALYSIS - CORRECTED APPROACH
    
    % Initialize behavior states only for requested behaviors
    behavior_states = struct();
    for i = 1:length(behaviors)
        behavior_states.(behaviors{i}) = false(n_frames, 1);
    end
    
    % Initialize features structure
    features = struct();
    
    %% STEP 1: REARING DETECTION (Independent - Z-axis behavior)
    if ismember('rearing', behaviors)
        % Rearing is primarily vertical behavior, independent of XY movement
        % Uses non-aligned data to capture true elevation
        tail_base_z = markers_moving.Tail_base_(:,3);
        spineM_z = markers_moving.SpineM(:,3);
        spineF_z = markers_moving.SpineF(:,3);
        
        % Condition 1: Paws are higher than the tail base
        paws_up = (markers_moving.ForepawL(:,3) > tail_base_z) & ...
                  (markers_moving.ForepawR(:,3) > tail_base_z);
                  
        % Condition 2: Snout is higher than tail base and mid-spine
        snout_up = (markers_moving.Snout(:,3) > tail_base_z) & ...
                   (markers_moving.Snout(:,3) > spineM_z);
                   
        % Condition 3: Ears are higher than tail base and mid-spine
        ears_up = (markers_moving.EarL(:,3) > tail_base_z) & ...
                  (markers_moving.EarR(:,3) > tail_base_z) & ...
                  (markers_moving.EarL(:,3) > spineM_z) & ...
                  (markers_moving.EarR(:,3) > spineM_z);
        
        % Condition 4: Absolute height threshold to filter false positives
        median_spineF_height = median(spineF_z);
        height_threshold = median_spineF_height + 1.5 * mad(spineF_z);
        spine_elevated = spineF_z > height_threshold;
        
        % Combine all conditions to define rearing state
        behavior_states.rearing = paws_up & snout_up & ears_up & spine_elevated;
    end
    
    %% STEP 2: ENHANCED XY MOVEMENT ANALYSIS (for remaining frames)
    if ismember('walking', behaviors)
        % Calculate horizontal body movement for non-rearing frames using center of mass
        mean_pos = mocapstruct.aligned_mean_position(frame_indices,1:2);
        mean_pos = [movmean(mean_pos(:,1),20), movmean(mean_pos(:,2), 20)];
        % Calculate kinematic features
        % 1. Velocity (displacement)
        velocity_xy = [0; sqrt(sum(diff(mean_pos).^2,2))];
        
        % 2. Acceleration (change in velocity magnitude)
        acceleration_magnitude = [0; 0; abs(diff(diff(velocity_xy)))];
        
        % 3. Smooth signals to reduce noise
        win_velocity = 2; % Short window for velocity to be responsive
        win_acceleration = 5; % Slightly longer window for acceleration to reduce noise
        
        velocity_smooth = movmean(velocity_xy, win_velocity);
        acceleration_smooth = movmean(acceleration_magnitude, win_acceleration);
        
        % 4. Calculate directional consistency (to distinguish walking from jitter)
        % Jitter = high acceleration, random directions
        % Walking = sustained direction, lower but consistent acceleration
        
        % Calculate direction changes (angular velocity)
        pos_diff = diff(mean_pos);
        pos_diff(vecnorm(pos_diff, 2, 2) < 0.1, :) = 0; % Ignore tiny movements
        
        direction_changes = zeros(size(velocity_xy));
        for i = 3:length(velocity_xy)
            if vecnorm(pos_diff(i-2,:)) > 0.1 && vecnorm(pos_diff(i-1,:)) > 0.1
                % Calculate angle between consecutive movement vectors
                v1 = pos_diff(i-2,:);
                v2 = pos_diff(i-1,:);
                cos_angle = dot(v1,v2) / (norm(v1) * norm(v2));
                cos_angle = max(-1, min(1, cos_angle)); % Clamp to [-1,1] for numerical stability
                direction_changes(i) = abs(acos(cos_angle));
            end
        end
        
        direction_consistency = movmean(direction_changes, 7);
        
        % 5. WALKING DETECTION with multiple criteria
        
        % Criterion 1: Velocity-based (for obvious walking)
        velocity_threshold = prctile(velocity_smooth, 70); % 70th percentile
        velocity_walking = velocity_smooth > velocity_threshold;
        
        % Criterion 2: Sustained movement (for slow walking)
        % Look for periods of consistent above-baseline movement
        baseline_velocity = prctile(velocity_smooth, 55); % Lower baseline for slow walks
        sustained_movement = movmean(velocity_smooth > baseline_velocity, 9) > 0.6; % 60% of 9-frame window
        
        % Criterion 3: Acceleration-based filtering to remove jitter
        % Jitter typically has high acceleration with random directions
        % Walking has moderate acceleration with consistent direction
        acceleration_threshold = prctile(acceleration_smooth, 10); % High acceleration threshold
        high_direction_change = direction_consistency > pi/3; % > 60 degrees average direction change
        
        % Identify jitter: high acceleration + high direction changes + low sustained velocity
        jitter_frames = (acceleration_smooth > acceleration_threshold) & ...
            (velocity_smooth < prctile(velocity_smooth, 33));
                       %high_direction_change & ...
                       
        
        % Criterion 4: Slow but directional movement
        % Detect slow walks: moderate velocity + low direction changes + sustained movement
        slow_walk_frames = (velocity_smooth > prctile(velocity_smooth, 40)) & ...
                          (direction_consistency < pi/6) & ... % < 30 degrees direction change
                          sustained_movement;
        
        % Combine all walking criteria and exclude jitter
        walking_combined = (velocity_walking | sustained_movement | slow_walk_frames) & ~jitter_frames;
        
        % Apply temporal smoothing to reduce isolated frame noise
        walking_smoothed = movmean(double(walking_combined), 5) > 0.4; % 40% of 5-frame window
        
        % Final walking state (exclude rearing frames if rearing is being analyzed)
        if ismember('rearing', behaviors)
            behavior_states.walking = walking_smoothed & ~behavior_states.rearing;
        else
            behavior_states.walking = walking_smoothed;
        end
        
        % Store movement features with enhanced debugging info
        features.velocity_xy = velocity_smooth;
        features.acceleration_magnitude = acceleration_smooth;
        features.direction_consistency = direction_consistency;
        features.velocity_threshold = velocity_threshold;
        features.baseline_velocity = baseline_velocity;
        features.jitter_frames = jitter_frames;
        features.slow_walk_frames = slow_walk_frames;
        features.walking_combined = walking_combined;
        features.walking_smoothed = walking_smoothed;
    end
    
    
    %% STEP 3: LEFT PAW LICKING DETECTION (for non-rearing frames)
    if ismember('left_paw_licking', behaviors)
        % Enhanced left paw attention detection using sophisticated algorithm
        pain_params = struct();
        pain_params.velocity_percentile = 50;
        pain_params.min_bout_duration = 0.01;
        pain_params.sampling_rate = 100;
        pain_params.smoothing_window = 5;
        pain_params.z_smoothing_window = 5;
        pain_params.z_threshold_percentile = 99.3;
        pain_params.direction_threshold = 75;
        pain_params.cycle_window = [-0.5 0.5];
        pain_params.n_normalized_points = 200;
        pain_params.interp_method = 'pchip';
        pain_params.markers_to_study = {'WristL', 'WristR', 'KneeL', 'KneeR','AnkleL', 'AnkleR', 'HindpawL', 'HindpawR'};
        
        % Apply sophisticated pain detection algorithm
        [pain_frames_strict, pain_metrics] = detectPainPhenotypes(markers_moving, markers, pain_params);
        
        % Add simple proximity-based detection for potential missed frames
        com_pos = (markers_moving.SpineF + markers_moving.SpineM + markers_moving.Tail_base_) / 3;
        velocity = [0; vecnorm(diff(com_pos), 2, 2)];
        smooth_velocity = movmean(velocity, 5);
        low_movement_frames = smooth_velocity < prctile(smooth_velocity, 50);
        
        % Secondary pain detection criteria
        snout_hindpawL_dist = vecnorm(markers.Snout - markers.HindpawL, 2, 2);
        close_to_left_hindpaw = snout_hindpawL_dist < (body_length * 0.4);
        
        % Left-leaning posture (from aligned data)
        snout_x = markers.Snout(:,1);
        left_lean = snout_x < -5; % More than 5mm left lean
        
        % Check relative positioning
        snout_y = markers.Snout(:,2);
        hindpaw_y = markers.HindpawL(:,2);
        near_left_side = abs(snout_y - hindpaw_y) < (body_length * 0.2);
        
        % Combine detection methods
        pain_frames_secondary = low_movement_frames & close_to_left_hindpaw & ...
                               left_lean & near_left_side;
        
        % Combine strict and secondary detection (logical OR), but exclude rearing if analyzed
        if ismember('rearing', behaviors)
            pain_frames_combined = (logical(pain_frames_strict) | pain_frames_secondary) & ~behavior_states.rearing;
        else
            pain_frames_combined = logical(pain_frames_strict) | pain_frames_secondary;
        end
        behavior_states.left_paw_licking = pain_frames_combined;
        
        % Store pain metrics as features
        features.pain_metrics = pain_metrics;
        features.pain_frames_strict = double(pain_frames_strict);
        features.pain_frames_secondary = double(pain_frames_secondary);
    end
    
    %% QUIET AND GROOMING DETECTION
    % Initialize variables for quiet detection if needed by any behavior
    if ismember('quiet', behaviors) || ismember('grooming', behaviors)
        % Need pain_metrics for is_still calculation
        if ~ismember('left_paw_licking', behaviors)
            % If left_paw_licking wasn't analyzed, we need to compute pain_metrics for is_still
            pain_params = struct();
            pain_params.velocity_percentile = 50;
            pain_params.min_bout_duration = 0.01;
            pain_params.sampling_rate = 100;
            pain_params.smoothing_window = 5;
            pain_params.z_smoothing_window = 5;
            pain_params.z_threshold_percentile = 99.3;
            pain_params.direction_threshold = 75;
            pain_params.cycle_window = [-0.5 0.5];
            pain_params.n_normalized_points = 200;
            pain_params.interp_method = 'pchip';
            pain_params.markers_to_study = {'WristL', 'WristR', 'KneeL', 'KneeR','AnkleL', 'AnkleR', 'HindpawL', 'HindpawR'};
            
            [~, pain_metrics] = detectPainPhenotypes(markers_moving, markers, pain_params);
        end
        
        % Calculate quiet frames
        if ismember('rearing', behaviors)
            is_quiet = pain_metrics.is_still & ~behavior_states.rearing;
        else
            is_quiet = pain_metrics.is_still;
        end
    end
    
    %% GROOMING DETECTION
    if ismember('grooming', behaviors)
        snout_pos = markers_moving.Snout(:, :);
        forepawL_pos = markers_moving.ForepawL(:, :);
        forepawR_pos = markers_moving.ForepawR(:, :);
        hindpawL_pos = markers_moving.HindpawL(:, :);

        snout_velocity = [0;vecnorm(diff(snout_pos), 2, 2)];

         % Proximity to FORE paws vs HIND paws
         dist_to_forepawL = vecnorm(snout_pos - forepawL_pos, 2, 2);
         dist_to_forepawR = vecnorm(snout_pos - forepawR_pos, 2, 2);
         dist_to_hindpawL = vecnorm(snout_pos - hindpawL_pos, 2, 2);

         min_forepaw_dist = min(dist_to_forepawL, dist_to_forepawR);

        % Grooming criteria - must be closer to FOREpaws than HINDpaws
            close_to_forepaws = (min_forepaw_dist) < (body_length * 0.7); % Within 70% of body length
            closer_to_fore_than_hind = (min_forepaw_dist) < mean(dist_to_hindpawL); % Closer to fore than hind
            high_snout_activity = (snout_velocity) > prctile(snout_velocity, 90);
                                 
        is_grooming_bin = close_to_forepaws & closer_to_fore_than_hind & high_snout_activity;
        behavior_states.grooming = is_grooming_bin & is_quiet;
    end

    % %% STEP 4: GROOMING AND QUIET DETECTION (for remaining frames)
    % % Use time-bin analysis for stationary behaviors (non-rearing, non-walking, non-licking)
    % 
    % % Define analysis parameters
    % params_behavior = struct();
    % params_behavior.sampling_rate = 100; % Hz
    % params_behavior.time_bin_size = 1.0; % seconds
    % params_behavior.bin_frames = round(params_behavior.time_bin_size * params_behavior.sampling_rate);
    % params_behavior.overlap = 0.5; % 50% overlap between bins
    % params_behavior.step_size = round(params_behavior.bin_frames * (1 - params_behavior.overlap));
    % 
    % % Analyze in time bins
    % for start_frame = 1:params_behavior.step_size:(n_frames - params_behavior.bin_frames + 1)
    %     end_frame = min(start_frame + params_behavior.bin_frames - 1, n_frames);
    %     bin_indices = start_frame:end_frame;
    % 
    %     % Skip if already classified as rearing, walking, or left paw licking
    %     if any(behavior_states.rearing(bin_indices)) || ...
    %        any(behavior_states.walking(bin_indices)) || ...
    %        any(behavior_states.left_paw_licking(bin_indices))
    %         continue;
    %     end
    % 
    %     % --- GROOMING DETECTION ---
    %     % Fine motor activity: snout-FOREPAW proximity + repetitive movements
    %     snout_pos = markers_moving.Snout(bin_indices, :);
    %     forepawL_pos = markers_moving.ForepawL(bin_indices, :);
    %     forepawR_pos = markers_moving.ForepawR(bin_indices, :);
    %     hindpawL_pos = markers_moving.HindpawL(bin_indices, :);
    % 
    %     % Proximity to FORE paws vs HIND paws
    %     dist_to_forepawL = vecnorm(snout_pos - forepawL_pos, 2, 2);
    %     dist_to_forepawR = vecnorm(snout_pos - forepawR_pos, 2, 2);
    %     dist_to_hindpawL = vecnorm(snout_pos - hindpawL_pos, 2, 2);
    % 
    %     min_forepaw_dist = min(dist_to_forepawL, dist_to_forepawR);
    % 
    %     % High-frequency snout movement (grooming signature)
    %     snout_velocity = vecnorm(diff(snout_pos), 2, 2);
    % 
    %     % Grooming criteria - must be closer to FOREpaws than HINDpaws
    %     close_to_forepaws = mean(min_forepaw_dist) < (body_length * 0.3); % Within 30% of body length
    %     closer_to_fore_than_hind = mean(min_forepaw_dist) < mean(dist_to_hindpawL); % Closer to fore than hind
    %     high_snout_activity = mean(snout_velocity) > prctile(snout_velocity, 70) && ...
    %                          std(snout_velocity) > prctile(snout_velocity, 60); % Variable, active movement
    % 
    %     is_grooming_bin = close_to_forepaws && closer_to_fore_than_hind && high_snout_activity;
    % 
    %     if is_grooming_bin
    %         behavior_states.grooming(bin_indices) = true;
    %         continue; % Don't classify as quiet if grooming
    %     end
    % 
    %     % --- QUIET DETECTION ---
    %     % Overall body stillness + low posture + no transitions
    %     spineM_pos = markers_moving.SpineM(bin_indices, :);
    %     spineF_pos = markers_moving.SpineF(bin_indices, :);
    % 
    %     % Body center of mass movement
    %     com_movement = vecnorm(diff(spineM_pos), 2, 2);
    %     body_stillness = mean(com_movement) < .5; % < 5mm average movement
    % 
    %     % Low acceleration (no sudden movements)
    %     body_accel = vecnorm(diff(diff(spineM_pos)), 2, 2);
    %     low_acceleration = mean(body_accel) < 2; % < 2mm/frame²
    % 
    %     % Check for velocity transitions (sign of behavioral change)
    %     velocity_consistency = std(com_movement) < 2; % Low variability in movement
    % 
    %     % Low posture (not elevated)
    %     low_posture = mean(spineF_pos(:,3)) < (median(markers_moving.SpineF(:,3)) + 10); % Within 10mm of median height
    % 
    %     % Minimal snout movement (no active sniffing/exploration)
    %     snout_stillness = mean(snout_velocity) < 3; % < 3mm/frame average
    % 
    %     % Check temporal context - avoid transitions
    %     % Look at frames before and after this bin
    %     context_before = max(1, start_frame - 25):start_frame-1; % 0.25s before
    %     context_after = end_frame+1:min(n_frames, end_frame + 25); % 0.25s after
    % 
    %     transitioning = false;
    %     if ~isempty(context_before) && ~isempty(context_after)
    %         % Check if there's a significant change in movement before/after
    %         movement_before = vecnorm(diff(markers_moving.SpineM(context_before, :)), 2, 2);
    %         movement_after = vecnorm(diff(markers_moving.SpineM(context_after, :)), 2, 2);
    %         current_movement = com_movement;
    % 
    %         % Detect transitions by comparing movement levels
    %         movement_change_before = abs(mean(movement_before) - mean(current_movement)) > 3;
    %         movement_change_after = abs(mean(movement_after) - mean(current_movement)) > 3;
    %         transitioning = movement_change_before || movement_change_after;
    %     end
    % 
    %     % Additional check: detect if rearing or walking occurs nearby
    %     nearby_active_behavior = false;
    %     context_window = max(1, start_frame - 50):min(n_frames, end_frame + 50); % 0.5s window
    %     if any(behavior_states.rearing(context_window)) || any(behavior_states.walking(context_window))
    %         nearby_active_behavior = true;
    %     end
    % 
    %     % Combine quiet criteria - must be still AND not transitioning
    %     is_quiet_bin = body_stillness && low_acceleration && low_posture && ...
    %                   snout_stillness && velocity_consistency && ...
    %                   ~transitioning && ~nearby_active_behavior;
    % 
    %     if is_quiet_bin
    %         behavior_states.quiet(bin_indices) = true;
    %     end
    % end

    
    %% FINALIZE BEHAVIORAL STATES WITH HIERARCHICAL LOGIC
    % Apply hierarchical logic based on analyzed behaviors
    
    % Define quiet state if analyzed
    if ismember('quiet', behaviors)
        quiet_base = is_quiet;
        % Exclude other behaviors from quiet
        if ismember('grooming', behaviors)
            quiet_base = quiet_base & ~behavior_states.grooming;
        end
        if ismember('left_paw_licking', behaviors)
            quiet_base = quiet_base & ~behavior_states.left_paw_licking;
        end
        behavior_states.quiet = quiet_base;
    end
    
    % Adjust walking to exclude quiet frames if both are analyzed
    if ismember('walking', behaviors) && ismember('quiet', behaviors)
        behavior_states.walking = behavior_states.walking & ~is_quiet;
    end
    
    
    %% FINALIZE FEATURES AND ENSURE MUTUAL EXCLUSIVITY
    % 
    % % Store additional features for analysis
    % features.body_length = body_length;
    % 
    % % Calculate continuous features for backward compatibility
    % paw_head_dist_L = vecnorm(markers.ForepawL - markers.Snout, 2, 2);
    % paw_head_dist_R = vecnorm(markers.ForepawR - markers.Snout, 2, 2);
    % features.grooming_proximity = 1 ./ min(paw_head_dist_L, paw_head_dist_R);
    % 
    % dist_snout_pawL = vecnorm(markers.Snout - markers.HindpawL, 2, 2);
    % features.left_paw_proximity = 1 ./ dist_snout_pawL;
    % 
    % % Enforce hierarchical mutual exclusivity
    % % Priority: Rearing > Walking > Left Paw Licking > Grooming > Quiet
    % 
    % % First level: Calculate moving state for summary
    % is_moving = behavior_states.rearing | behavior_states.walking;
    % 
    % % Enforce priorities
    % behavior_states.walking = behavior_states.walking & ~behavior_states.rearing;
    % behavior_states.grooming = behavior_states.grooming & ~behavior_states.left_paw_licking & ~behavior_states.rearing & ~behavior_states.walking;
    % behavior_states.quiet = behavior_states.quiet & ~behavior_states.left_paw_licking & ~behavior_states.grooming & ~behavior_states.rearing & ~behavior_states.walking;
    % 
    % % Calculate total classified frames
    % total_classified = behavior_states.rearing | behavior_states.walking | ...
    %                   behavior_states.left_paw_licking | behavior_states.grooming | ...
    %                   behavior_states.quiet;
    % 
    % % Add summary statistics
    % features.behavior_summary = struct();
    % features.behavior_summary.total_frames = n_frames;
    % features.behavior_summary.moving_frames = sum(is_moving);
    % features.behavior_summary.rearing_frames = sum(behavior_states.rearing);
    % features.behavior_summary.walking_frames = sum(behavior_states.walking);
    % features.behavior_summary.left_paw_licking_frames = sum(behavior_states.left_paw_licking);
    % features.behavior_summary.grooming_frames = sum(behavior_states.grooming);
    % features.behavior_summary.quiet_frames = sum(behavior_states.quiet);
    % features.behavior_summary.unclassified_frames = n_frames - sum(total_classified);
    % 
    % % Add percentages
    % features.behavior_summary.percent_moving = (features.behavior_summary.moving_frames / n_frames) * 100;
    % features.behavior_summary.percent_rearing = (features.behavior_summary.rearing_frames / n_frames) * 100;
    % features.behavior_summary.percent_walking = (features.behavior_summary.walking_frames / n_frames) * 100;
    % features.behavior_summary.percent_left_paw_licking = (features.behavior_summary.left_paw_licking_frames / n_frames) * 100;
    % features.behavior_summary.percent_grooming = (features.behavior_summary.grooming_frames / n_frames) * 100;
    % features.behavior_summary.percent_quiet = (features.behavior_summary.quiet_frames / n_frames) * 100;
    % features.behavior_summary.percent_unclassified = (features.behavior_summary.unclassified_frames / n_frames) * 100;

end