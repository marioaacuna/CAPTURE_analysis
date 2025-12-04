function [pain_frames, metrics, confidence_scores] = detectPainPhenotypes_v2(non_aligned_mocap, aligned_mocap, params)
    % Refined multi-modal pain detection with bilateral features and adaptive thresholds
    
    num_frames = size(non_aligned_mocap.SpineM, 1);
    
    % Step 1: Enhanced movement-based frame filtering with centered difference
    spine_vel = calculateVelocityRobust(non_aligned_mocap.SpineM, params);
    movement_threshold = prctile(spine_vel, params.velocity_percentile); % Configurable percentile
    candidate_frames = spine_vel < movement_threshold;
    
    % Step 2: Calculate baseline statistics for adaptive thresholds
    baseline_stats = calculateBaselineStats(aligned_mocap, candidate_frames, params);
    
    % Initialize outputs
    pain_frames = false(num_frames, 1);
    confidence_scores = zeros(num_frames, 1);
    
    % Initialize comprehensive metrics
    metrics = initializeMetrics(num_frames);
    
    % Step 3: Extract and smooth temporal features
    raw_features = extractAllFrameFeatures(aligned_mocap);
    % smoothed_features = temporalFeatureSmoothing(raw_features, params);
    smoothed_features = raw_features;
    % Step 4: Vectorized multi-modal pain behavior detection
    % Calculate probabilities for all frames at once
    % prob_paw_licking = assessPawLickingBilateral(smoothed_features, baseline_stats);
    prob_paw_licking = assessPawLickingLeftOnly(smoothed_features, baseline_stats);
    prob_guarding = assessGuardingAdaptive(smoothed_features, baseline_stats);
    prob_asymmetry = assessAsymmetryAdaptive(smoothed_features, baseline_stats);
    prob_hunched = assessHunchedPostureAdaptive(smoothed_features, baseline_stats);
    prob_grooming = assessGroomingBilateral(smoothed_features, baseline_stats);
    
    % Combine probabilities using weighted ensemble
    %weights = [0.50, 0.25, 0.10, 0.15, 0.0];
    %combined_prob = weights(1) * prob_paw_licking + ...
    %                weights(2) * prob_guarding + ...
    %                weights(3) * prob_asymmetry + ...
    %                weights(4) * prob_hunched + ...
    %                weights(5) * prob_grooming;
    

    % % Combine probabilities by taking the maximum score across key pain phenotypes.
    % % This ensures that a strong signal from any single behavior is not diluted.
    % % Paw licking, guarding, asymmetry, and hunching are considered primary pain indicators.
    % pain_phenotype_probs = [prob_paw_licking, prob_guarding, prob_asymmetry, prob_hunched];
    % combined_prob = max(pain_phenotype_probs, [], 2);
    
    % Stage 1: Primary behavior detection
    primary_detection = max(prob_paw_licking, prob_hunched); % Your strongest signals
    
    % Stage 2: Confirmation from secondary behaviors  
    secondary_evidence = mean([prob_guarding, prob_asymmetry],2);
    
    % Combine with interaction
    combined_prob = primary_detection * 0.7 + (0.3 * secondary_evidence);


    % Apply only to candidate frames
    confidence_scores(candidate_frames) = combined_prob(candidate_frames);


    % Apply only to candidate frames
    confidence_scores(candidate_frames) = combined_prob(candidate_frames);
    pain_frames(candidate_frames) = combined_prob(candidate_frames) > params.pain_threshold;
    
    % Store all metrics at once
    metrics = updateMetricsVectorized(metrics, smoothed_features, ...
                                    [prob_paw_licking, prob_guarding, prob_asymmetry, prob_hunched, prob_grooming], ...
                                    candidate_frames);
    
    % Step 5: Adaptive temporal smoothing
    % pain_frames = adaptiveTemporalSmoothing(pain_frames, params);
    
end

function speed_smooth = calculateVelocityRobust(positions, params)
    % Calculate velocity using centered difference with NaN protection
    
    num_frames = size(positions, 1);
    vel_xy = zeros(num_frames, 2);
    
    % Centered difference for interior points
    for i = 2:num_frames-1
        vel_xy(i, :) = (positions(i+1, 1:2) - positions(i-1, 1:2)) * params.sampling_rate / 2;
    end
    
    % Forward/backward difference for endpoints
    vel_xy(1, :) = (positions(2, 1:2) - positions(1, 1:2)) * params.sampling_rate;
    vel_xy(end, :) = (positions(end, 1:2) - positions(end-1, 1:2)) * params.sampling_rate;
    
    % Calculate speed and handle NaN values
    speed = sqrt(sum(vel_xy.^2, 2));
    speed(isnan(speed)) = 0; % Replace NaN with 0
    
    % Smooth with robust method
    speed_smooth = movmedian(speed, params.smoothing_window);
end

function baseline_stats = calculateBaselineStats(aligned_mocap, candidate_frames, params)
    % Calculate adaptive baseline statistics from calm periods
    
    % Use bottom 25% of candidate frames as baseline (calmest periods)
    baseline_indices = find(candidate_frames);
    n_baseline = round(length(baseline_indices) * 0.25);
    baseline_frames = baseline_indices(1:n_baseline);
    
    baseline_stats = struct();
    
    % Calculate baseline means and stds for key metrics
    snout_paw_L_dists = arrayfun(@(f) norm(aligned_mocap.Snout(f,:) - aligned_mocap.HindpawL(f,:)), baseline_frames);
    snout_paw_R_dists = arrayfun(@(f) norm(aligned_mocap.Snout(f,:) - aligned_mocap.HindpawR(f,:)), baseline_frames);
    
    baseline_stats.snout_paw_L_mean = mean(snout_paw_L_dists);
    baseline_stats.snout_paw_L_std = std(snout_paw_L_dists);
    baseline_stats.snout_paw_R_mean = mean(snout_paw_R_dists);
    baseline_stats.snout_paw_R_std = std(snout_paw_R_dists);
    
    % Limb angle baselines
    left_angles = arrayfun(@(f) calculateLimbAngle(aligned_mocap.KneeL(f,:), aligned_mocap.AnkleL(f,:), aligned_mocap.HindpawL(f,:)), baseline_frames);
    right_angles = arrayfun(@(f) calculateLimbAngle(aligned_mocap.KneeR(f,:), aligned_mocap.AnkleR(f,:), aligned_mocap.HindpawR(f,:)), baseline_frames);
    
    baseline_stats.limb_angle_L_mean = mean(left_angles);
    baseline_stats.limb_angle_L_std = std(left_angles);
    baseline_stats.limb_angle_R_mean = mean(right_angles);
    baseline_stats.limb_angle_R_std = std(right_angles);
    
    % Asymmetry baselines
    asymmetries = abs(left_angles - right_angles);
    baseline_stats.asymmetry_mean = mean(asymmetries);
    baseline_stats.asymmetry_std = std(asymmetries);
    
    % Trunk curvature baseline
    curvatures = arrayfun(@(f) calculateTrunkCurvature(aligned_mocap, f), baseline_frames);
    baseline_stats.curvature_mean = mean(curvatures);
    baseline_stats.curvature_std = std(curvatures);
end

function all_features = extractAllFrameFeatures(aligned_mocap)
    % Extract features for all frames at once for temporal smoothing
    
    num_frames = size(aligned_mocap.Snout, 1);
    all_features = struct();
    
    % Preallocate arrays
    all_features.snout_to_hindpaw_L = zeros(num_frames, 1);
    all_features.snout_to_hindpaw_R = zeros(num_frames, 1);
    all_features.forepaw_to_hindpaw_L = zeros(num_frames, 1);
    all_features.forepaw_to_hindpaw_R = zeros(num_frames, 1);
    all_features.trunk_curvature = zeros(num_frames, 1);
    all_features.lateral_deviation = zeros(num_frames, 1);
    all_features.spine_elevation = zeros(num_frames, 1);
    all_features.left_hindlimb_angle = zeros(num_frames, 1);
    all_features.right_hindlimb_angle = zeros(num_frames, 1);
    all_features.limb_asymmetry = zeros(num_frames, 1);
    all_features.paw_height_asymmetry = zeros(num_frames, 1);
    all_features.head_orientation = zeros(num_frames, 1);
    all_features.body_axis_deviation = zeros(num_frames, 1);
    
    % Extract features for all frames
    for frame = 1:num_frames
        all_features.snout_to_hindpaw_L(frame) = norm(aligned_mocap.Snout(frame,:) - aligned_mocap.HindpawL(frame,:));
        all_features.snout_to_hindpaw_R(frame) = norm(aligned_mocap.Snout(frame,:) - aligned_mocap.HindpawR(frame,:));
        all_features.forepaw_to_hindpaw_L(frame) = norm(aligned_mocap.ForepawL(frame,:) - aligned_mocap.HindpawL(frame,:));
        all_features.forepaw_to_hindpaw_R(frame) = norm(aligned_mocap.ForepawR(frame,:) - aligned_mocap.HindpawR(frame,:));
        
        all_features.trunk_curvature(frame) = calculateTrunkCurvature(aligned_mocap, frame);
        all_features.lateral_deviation(frame) = calculateLateralDeviation(aligned_mocap, frame);
        all_features.spine_elevation(frame) = calculateSpineElevation(aligned_mocap, frame);
        
        all_features.left_hindlimb_angle(frame) = calculateLimbAngle(aligned_mocap.KneeL(frame,:), ...
                                                                    aligned_mocap.AnkleL(frame,:), ...
                                                                    aligned_mocap.HindpawL(frame,:));
        all_features.right_hindlimb_angle(frame) = calculateLimbAngle(aligned_mocap.KneeR(frame,:), ...
                                                                     aligned_mocap.AnkleR(frame,:), ...
                                                                     aligned_mocap.HindpawR(frame,:));
        
        all_features.limb_asymmetry(frame) = abs(all_features.left_hindlimb_angle(frame) - all_features.right_hindlimb_angle(frame));
        all_features.paw_height_asymmetry(frame) = abs(aligned_mocap.HindpawL(frame,3) - aligned_mocap.HindpawR(frame,3));
        
        all_features.head_orientation(frame) = calculateHeadOrientation(aligned_mocap, frame);
        all_features.body_axis_deviation(frame) = calculateBodyAxisDeviation(aligned_mocap, frame);
    end
end

function smoothed_features = temporalFeatureSmoothing(raw_features, params)
    % Apply temporal smoothing to reduce measurement jitter
    
    smoothed_features = raw_features;
    field_names = fieldnames(raw_features);
    
    for i = 1:length(field_names)
        field = field_names{i};
        % Use rolling median for robust smoothing
        smoothed_features.(field) = movmedian(raw_features.(field), params.feature_smoothing_window);
    end
end

% function prob = assessPawLickingBilateral(features, baseline_stats)
%     % Vectorized bilateral assessment of paw licking behavior
% 
%     % Z-scores for both paws
%     z_score_L = (features.snout_to_hindpaw_L - baseline_stats.snout_paw_L_mean) / baseline_stats.snout_paw_L_std;
%     z_score_R = (features.snout_to_hindpaw_R - baseline_stats.snout_paw_R_mean) / baseline_stats.snout_paw_R_std;
% 
%     % Sigmoid scoring for both sides
%     dist_score_L = sigmoid(-z_score_L, -2, 0.5);
%     dist_score_R = sigmoid(-z_score_R, -2, 0.5);
% 
%     % Head orientation score
%     orientation_score = sigmoid(-features.head_orientation, -15, 5);
% 
%     % Proximity scores for both sides
%     proximity_score_L = sigmoid(-features.forepaw_to_hindpaw_L + 20, 0, 3);
%     proximity_score_R = sigmoid(-features.forepaw_to_hindpaw_R + 20, 0, 3);
% 
%     % Take maximum across sides (element-wise)
%     max_dist_score = max(dist_score_L, dist_score_R);
%     max_proximity_score = max(proximity_score_L, proximity_score_R);
% 
%     % Weighted combination
%     prob = 0.5 * max_dist_score + 0.3 * orientation_score + 0.2 * max_proximity_score;
%     prob = max(0, min(1, prob));
% end


function prob = assessPawLickingLeftOnly(features, baseline_stats)
    % Z-score for left paw distance (keep this - it works)
    z_score_L = (features.snout_to_hindpaw_L - baseline_stats.snout_paw_L_mean) / baseline_stats.snout_paw_L_std;
    dist_score_L = sigmoid(-z_score_L, 3, 0.3); 

    % Make other components more stringent
    proximity_score = sigmoid(-features.forepaw_to_hindpaw_L + 30, 0, .1); % Tighter threshold
    % lateral_lean_score = sigmoid(features.lateral_deviation, 0, .1); % Require significant lean
    orientation_score = sigmoid(-features.head_orientation, 45, .1); % Require clear head turn
    
    % Use multiplicative logic instead of additive (requires ALL conditions)
    prob = dist_score_L * 0.6 + ...
              (dist_score_L .* proximity_score .* orientation_score) * 0.4; %           % (dist_score_L .* proximity_score .* lateral_lean_score .* orientation_score) * 0.4;

 
    prob = max(0, min(1, prob));
end

% function prob = assessPawLickingLeftOnly(features, baseline_stats)
%     % Vectorized assessment of left paw licking behavior (injured side)
% 
%     % Z-score for left paw only
%     z_score_L = (features.snout_to_hindpaw_L - baseline_stats.snout_paw_L_mean) / baseline_stats.snout_paw_L_std;
% 
%     % Sigmoid scoring for left side
%     dist_score_L = sigmoid(-z_score_L, -2, 0.5);
% 
%     % Head orientation score (negative = left turn toward injured paw)
%     orientation_score = sigmoid(-features.head_orientation, -15, 5);
% 
%     % Proximity score for left side only
%     proximity_score_L = sigmoid(-features.forepaw_to_hindpaw_L + 20, 0, 3);
% 
%     % Weighted combination - focus on injured left paw
%     prob = 0.5 * dist_score_L + 0.3 * orientation_score + 0.2 * proximity_score_L;
%     prob = max(0, min(1, prob));
% end

% function prob = assessGuardingAdaptive(features, baseline_stats)
%     % Vectorized adaptive assessment of protective guarding behavior
% 
%   z_score_L = (features.left_hindlimb_angle - baseline_stats.limb_angle_L_mean) / baseline_stats.limb_angle_L_std;
%     z_score_R = (features.right_hindlimb_angle - baseline_stats.limb_angle_R_mean) / baseline_stats.limb_angle_R_std;
% 
%     limb_flexion_score = max(sigmoid(z_score_L, 1.5, 0.5), sigmoid(z_score_R, 1.5, 0.5));
% 
%     % Fix elevation scoring - use stability around baseline
%     elevation_z = zscore(features.spine_elevation);
%     elevation_stability_score = sigmoid((elevation_z), 0.5, 0.3); % Higher score for low deviation
% 
%     % Asymmetry score
%     asymmetry_z = (features.limb_asymmetry - baseline_stats.asymmetry_mean) / baseline_stats.asymmetry_std;
%     asymmetry_score = sigmoid(asymmetry_z, 2, 0.5);
% 
%     prob = 0.4 * limb_flexion_score + 0.4 * elevation_stability_score + 0.2 * asymmetry_score;
%     prob = max(0, min(1, prob));
% end

function prob = assessGuardingAdaptive(features, baseline_stats)
    % Classical guarding assessment for left hindpaw injury
    
    % 1. Weight shifting - body shifted away from injured left paw
    lateral_shift_z = zscore(features.lateral_deviation);
    weight_shift_score = sigmoid(lateral_shift_z, 0.5, 0.5); % Positive = shift right (away from left injury)
    
    % 2. Protective flexion - left hindlimb more flexed than baseline
    left_limb_z = (features.left_hindlimb_angle - baseline_stats.limb_angle_L_mean) / baseline_stats.limb_angle_L_std;
    protective_flexion_score = sigmoid(left_limb_z, 2.0, 0.15); % Higher flexion angle
    
    % 3. Limb elevation - left hindpaw lifted to avoid ground contact
    paw_height_asym_z = zscore(features.paw_height_asymmetry);
    limb_elevation_score = sigmoid(paw_height_asym_z, 1.0, 0.5); % Left paw higher than right
    
    % 4. Asymmetric loading - overall limb asymmetry increased
    asymmetry_z = (features.limb_asymmetry - baseline_stats.asymmetry_mean) / baseline_stats.asymmetry_std;
    asymmetric_loading_score = sigmoid(asymmetry_z, 1.5, 0.5); % Increased asymmetry
    
    % Weighted combination focusing on protective behaviors
    prob = 0.3 * weight_shift_score + ...
           0.3 * protective_flexion_score + ...
           0.20 * limb_elevation_score + ...
           0.20 * asymmetric_loading_score;
    
    prob = max(0, min(1, prob));
end


function prob = assessAsymmetryAdaptive(features, baseline_stats)
    % Vectorized adaptive assessment of postural asymmetry
    
    % Z-score for limb asymmetry
    asymmetry_z = (features.limb_asymmetry - baseline_stats.asymmetry_mean) / baseline_stats.asymmetry_std;
    limb_asym_score = sigmoid(asymmetry_z, 2, 0.5);
    
    lat_dev_z = zscore(features.lateral_deviation);
    lateral_asym_score = sigmoid(abs(lat_dev_z), 2, 3);
    height_asym_score = sigmoid(features.paw_height_asymmetry, 2, 2);
    
    prob = 0.5 * limb_asym_score + 0.3 * lateral_asym_score + 0.2 * height_asym_score;
    prob = max(0, min(1, prob));
end

function prob = assessHunchedPostureAdaptive(features, baseline_stats)
    % Vectorized adaptive assessment of hunched/arched back posture
    
    % Z-score for trunk curvature
    curvature_z = (features.trunk_curvature - baseline_stats.curvature_mean) / baseline_stats.curvature_std;
    curvature_score = sigmoid(-curvature_z, 0.5, 0.5);

    elevation_score = sigmoid(-zscore(features.spine_elevation), 0.2, 0.5);
    
    prob = 0.7 * curvature_score + 0.3 * elevation_score;
    prob = max(0, min(1, prob));
end

function prob = assessGroomingBilateral(features, baseline_stats)
    % Vectorized bilateral assessment of excessive grooming behavior
    
    % Take minimum distance (closest paw) element-wise
    min_proximity = min(features.snout_to_hindpaw_L, features.snout_to_hindpaw_R);
    
    proximity_score = sigmoid(-min_proximity + 18, 0, 4);
    
    prob = proximity_score * 0.8;
    prob = max(0, min(1, prob));
end

function smoothed_frames = adaptiveTemporalSmoothing(pain_frames, params)
    % Adaptive temporal smoothing based on frame rate
    
    % Scale window size based on sampling rate
    base_window_ms = 100; % 100ms base window
    window_frames = round(base_window_ms * params.sampling_rate / 1000);
    window_frames = max(3, min(window_frames, params.max_smoothing_window));
    
    % Morphological opening: remove isolated positives
    smoothed_frames = imopen(pain_frames, ones(window_frames, 1));
    
    % Fill small gaps (half window size)
    gap_fill_window = max(1, round(window_frames / 2));
    smoothed_frames = imclose(smoothed_frames, ones(gap_fill_window, 1));
end

function y = sigmoid(x, midpoint, steepness)
    % Sigmoid function for soft thresholding
    y = 1 ./ (1 + exp(-(x - midpoint) / steepness));
end

function metrics = initializeMetrics(num_frames)
    metrics = struct( ...
        'snout_paw_dist_L', zeros(num_frames,1), ...
        'snout_paw_dist_R', zeros(num_frames,1), ...
        'forepaw_hindpaw_dist_L', zeros(num_frames,1), ...
        'forepaw_hindpaw_dist_R', zeros(num_frames,1), ...
        'trunk_curvature', zeros(num_frames,1), ...
        'lateral_deviation', zeros(num_frames,1), ...
        'spine_elevation', zeros(num_frames,1), ...
        'left_hindlimb_angle', zeros(num_frames,1), ...
        'right_hindlimb_angle', zeros(num_frames,1), ...
        'limb_asymmetry', zeros(num_frames,1), ...
        'paw_height_asymmetry', zeros(num_frames,1), ...
        'head_orientation', zeros(num_frames,1), ...
        'body_axis_deviation', zeros(num_frames,1), ...
        'prob_paw_licking', zeros(num_frames,1), ...
        'prob_guarding', zeros(num_frames,1), ...
        'prob_asymmetry', zeros(num_frames,1), ...
        'prob_hunched', zeros(num_frames,1), ...
        'prob_grooming', zeros(num_frames,1), ...
        'is_processed', false(num_frames,1));
end

% Include all the helper functions from previous version
function curvature = calculateTrunkCurvature(aligned_mocap, frame)
    p1 = aligned_mocap.SpineF(frame,:);
    p2 = aligned_mocap.SpineM(frame,:);
    p3 = aligned_mocap.Tail_base_(frame,:);
    
    v1 = p2 - p1;
    v2 = p3 - p2;
    
    cos_angle = dot(v1, v2) / (norm(v1) * norm(v2));
    cos_angle = max(-1, min(1, cos_angle));
    
    angle = acosd(cos_angle);
    curvature = 180 - angle;
end

function deviation = calculateLateralDeviation(aligned_mocap, frame)
    spine_y = aligned_mocap.SpineF(frame, 2);
    deviation = spine_y;
end

function elevation = calculateSpineElevation(aligned_mocap, frame)
    spine_z = aligned_mocap.SpineM(frame, 3);
    paw_z_avg = mean([aligned_mocap.HindpawL(frame, 3), ...
                      aligned_mocap.HindpawR(frame, 3), ...
                      aligned_mocap.ForepawL(frame, 3), ...
                      aligned_mocap.ForepawR(frame, 3)]);
    elevation = spine_z - paw_z_avg;
end

function angle = calculateLimbAngle(knee, ankle, paw)
    v1 = knee - ankle;
    v2 = paw - ankle;
    cos_angle = dot(v1, v2) / (norm(v1) * norm(v2));
    cos_angle = max(-1, min(1, cos_angle));
    angle = acosd(cos_angle);
end

function orientation = calculateHeadOrientation(aligned_mocap, frame)
    head_vec = aligned_mocap.Snout(frame,:) - aligned_mocap.SpineF(frame,:);
    body_vec = aligned_mocap.SpineF(frame,:) - aligned_mocap.Tail_base_(frame,:);
    
    head_xy = head_vec(1:2);
    body_xy = body_vec(1:2);
    
    cross_prod = head_xy(1) * body_xy(2) - head_xy(2) * body_xy(1);
    dot_prod = dot(head_xy, body_xy);
    
    orientation = atan2d(cross_prod, dot_prod);
end

function deviation = calculateBodyAxisDeviation(aligned_mocap, frame)
    body_vec = aligned_mocap.SpineF(frame,:) - aligned_mocap.Tail_base_(frame,:);
    forward_vec = [1, 0, 0];
    
    body_xy = body_vec(1:2) / norm(body_vec(1:2));
    
    cos_angle = dot(body_xy, forward_vec(1:2));
    cos_angle = max(-1, min(1, cos_angle));
    
    deviation = acosd(cos_angle);
end

function metrics = updateMetricsVectorized(metrics, behavioral_features, probabilities, candidate_frames)
    % Vectorized update of metrics structure
    
    % Store behavioral features for all frames
    metrics.snout_paw_dist_L = behavioral_features.snout_to_hindpaw_L;
    metrics.snout_paw_dist_R = behavioral_features.snout_to_hindpaw_R;
    metrics.forepaw_hindpaw_dist_L = behavioral_features.forepaw_to_hindpaw_L;
    metrics.forepaw_hindpaw_dist_R = behavioral_features.forepaw_to_hindpaw_R;
    metrics.trunk_curvature = behavioral_features.trunk_curvature;
    metrics.lateral_deviation = behavioral_features.lateral_deviation;
    metrics.spine_elevation = behavioral_features.spine_elevation;
    metrics.left_hindlimb_angle = behavioral_features.left_hindlimb_angle;
    metrics.right_hindlimb_angle = behavioral_features.right_hindlimb_angle;
    metrics.limb_asymmetry = behavioral_features.limb_asymmetry;
    metrics.paw_height_asymmetry = behavioral_features.paw_height_asymmetry;
    metrics.head_orientation = behavioral_features.head_orientation;
    metrics.body_axis_deviation = behavioral_features.body_axis_deviation;
    
    % Store probability scores for all frames
    metrics.prob_paw_licking = probabilities(:,1);
    metrics.prob_guarding = probabilities(:,2);
    metrics.prob_asymmetry = probabilities(:,3);
    metrics.prob_hunched = probabilities(:,4);
    metrics.prob_grooming = probabilities(:,5);
    
    % Mark processed frames
    metrics.is_processed(candidate_frames) = true;
end