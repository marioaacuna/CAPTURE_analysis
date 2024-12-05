function [walking_bouts, metrics] = detectWalkingFrom2D(spineM, snout, params)
    % Added snout marker for direction detection
    if nargin < 3
        params = struct();
        params.sampling_rate = 100; 
        params.velocity_percentile = 75; 
        params.min_bout_duration = 0.1; 
        params.smoothing_window = 5;
        params.z_threshold_percentile = 85;
        params.z_smoothing_window = 10;
        params.direction_threshold = 60;    % Max angle deviation from heading (degrees)
    end
    
    % Extract coordinates
    xy_pos = spineM(:,1:2);
    z_pos = spineM(:,3);
    snout_xy = snout(:,1:2);
    
    % Calculate heading direction (from spineM to snout)
    heading_vector = snout_xy - xy_pos;
    heading_angle = atan2d(heading_vector(:,2), heading_vector(:,1));
    
    % Calculate movement direction
    movement_vector = diff(xy_pos);
    movement_vector = [zeros(1,2); movement_vector];  % Add initial zero vector
    movement_angle = atan2d(movement_vector(:,2), movement_vector(:,1));
    
    % Calculate angular difference between heading and movement
    angle_diff = abs(angdiff(deg2rad(heading_angle), deg2rad(movement_angle)));
    angle_diff_deg = rad2deg(angle_diff);
    
    % Identify forward movement (within threshold of heading direction)
    forward_movement = angle_diff_deg < params.direction_threshold;
    
    % Z movement processing (as before)
    z_vel = [0; diff(z_pos)] * params.sampling_rate;
    z_vel_smooth = movmean(abs(z_vel), params.z_smoothing_window);
    z_threshold = prctile(z_vel_smooth, params.z_threshold_percentile);
    rearing_periods = z_vel_smooth > z_threshold;
    rearing_periods = imdilate(rearing_periods, ones(params.sampling_rate, 1));
    
    % XY velocity processing
    vel_xy = diff(xy_pos) * params.sampling_rate;
    vel_xy = [zeros(1,2); vel_xy];
    speed = sqrt(sum(vel_xy.^2, 2));
    speed_smooth = movmean(speed, params.smoothing_window);
    velocity_threshold = prctile(speed_smooth, params.velocity_percentile);
    
    % Combine all criteria for walking detection
    walking_bouts = (speed_smooth > velocity_threshold) & ...  % Moving fast enough
                   ~rearing_periods & ...                      % Not rearing
                   forward_movement;                           % Moving forward
    
    % Filter short bouts
    walking_bouts = filterShortBouts(walking_bouts, params.min_bout_duration * params.sampling_rate);
    
    % Calculate metrics
    if nargout > 1
        metrics = struct();
        metrics.mean_speed = mean(speed_smooth(walking_bouts));
        metrics.max_speed = max(speed_smooth(walking_bouts));
        metrics.velocity_threshold = velocity_threshold;
        metrics.z_threshold = z_threshold;
        metrics.total_distance = sum(speed_smooth(walking_bouts))/params.sampling_rate;
        metrics.speed_smooth = speed_smooth;
        metrics.z_velocity = z_vel_smooth;
        metrics.rearing_periods = rearing_periods;
        metrics.heading_angle = heading_angle;
        metrics.movement_angle = movement_angle;
        metrics.angle_difference = angle_diff_deg;
        metrics.forward_movement = forward_movement;
        
        CC = bwconncomp(walking_bouts);
        bout_durations = cellfun(@length, CC.PixelIdxList)/params.sampling_rate;
        metrics.num_bouts = CC.NumObjects;
        metrics.mean_bout_duration = mean(bout_durations);
        metrics.bout_durations = bout_durations;
    end
end

function filtered = filterShortBouts(signal, min_duration)
    CC = bwconncomp(signal);
    for i = 1:CC.NumObjects
        if length(CC.PixelIdxList{i}) < min_duration
            signal(CC.PixelIdxList{i}) = false;
        end
    end
    filtered = signal;
end