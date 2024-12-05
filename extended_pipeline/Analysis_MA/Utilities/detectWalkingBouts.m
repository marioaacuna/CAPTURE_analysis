function [walking_bouts, metrics] = detectWalkingBouts(leg_angles, params)
    % Set default parameters
    if nargin < 2
        params = struct();
        params.sampling_rate = 100;
        params.hip_range = [-20, 20];    
        params.knee_range = [40, 100];   
        params.ankle_range = [50, 100];  
        params.window_size = 0.1;        % 100ms for quick steps
        params.min_bout_duration = 0.1;  % seconds
    end
    
    % Calculate velocities
    vel = diff(leg_angles) * params.sampling_rate;
    vel = [zeros(1,3); vel];
    total_ang_vel = sqrt(sum(vel.^2, 2));
    
    % Check angle ranges - using corrected indexing
    in_range = leg_angles(:,3) >= params.hip_range(1) & leg_angles(:,3) <= params.hip_range(2) & ...
               leg_angles(:,1) >= params.knee_range(1) & leg_angles(:,1) <= params.knee_range(2) & ...
               leg_angles(:,2) >= params.ankle_range(1) & leg_angles(:,2) <= params.ankle_range(2);
    
    % Set up for phase analysis
    window_size = round(params.window_size * params.sampling_rate);
    overlap = round(window_size * 0.5);
    
    % Calculate phase relationships between joints
    walking_scores = zeros(size(leg_angles,1), 1);
    for t = 1:size(leg_angles,1)-window_size
        window_data = leg_angles(t:t+window_size-1, :);
        
        % Cross-correlation between joints
        hip_knee_xcorr = xcorr(window_data(:,1), window_data(:,2), 'coeff');
        knee_ankle_xcorr = xcorr(window_data(:,2), window_data(:,3), 'coeff');
        
        % Get phase relationships
        [~, hip_knee_lag] = max(abs(hip_knee_xcorr));
        [~, knee_ankle_lag] = max(abs(knee_ankle_xcorr));
        
        % Score based on expected walking phase relationships
        phase_score = gaussmf(hip_knee_lag, [window_size/4, window_size]) * ...
                     gaussmf(knee_ankle_lag, [window_size/4, window_size]);
        
        % Check for alternating movement in window
        hip_vel_window = vel(t:t+window_size-1, 1);
        knee_vel_window = vel(t:t+window_size-1, 2);
        ankle_vel_window = vel(t:t+window_size-1, 3);
        
        alternating_pattern = any(diff(sign(hip_vel_window)) ~= 0) && ...
                            any(diff(sign(knee_vel_window)) ~= 0) && ...
                            any(diff(sign(ankle_vel_window)) ~= 0) && ...
                            mean(abs(hip_vel_window)) < mean(abs(knee_vel_window));
                            
        walking_scores(t) = phase_score * alternating_pattern;
    end
    
    % Pad the end
    walking_scores(end-window_size+1:end) = walking_scores(end-window_size);
    
    % Combine all criteria
    walking_bouts = in_range & ...           % Angles within range
                   walking_scores > 0.1 & ... % Good phase relationships
                   total_ang_vel > 7;        % Some movement happening
    
    % Filter short bouts
    walking_bouts = filterShortBouts(walking_bouts, params.min_bout_duration * params.sampling_rate);
    
    % Calculate metrics if needed
    if nargout > 1
        metrics = struct();
        metrics.walking_scores = walking_scores;
        metrics.in_range = in_range;
        metrics.total_velocity = total_ang_vel;
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