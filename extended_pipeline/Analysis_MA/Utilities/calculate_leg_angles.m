function angles = calculate_leg_angles(knee, ankle, paw)
    % this function calculates the angles between the segments of the leg
    % knee, ankle and paw are the 3D coordinates of the markers
    % output is a matrix with the angles between the segments
    % where the first column is the knee angle and the second column is the ankle angle
    thigh_vec = knee;  % Already relative to SpineM (origin)
    shank_vec = ankle - knee;
    foot_vec = paw - ankle;
    
    % Calculate cross products for all frames
    cross_knee = cross(thigh_vec, shank_vec, 2);  % '2' specifies operation along 2nd dimension
    cross_ankle = cross(-shank_vec, foot_vec, 2);  % Note the negative shank_vec for correct ankle direction
    
    % Get directions from y-component (sagittal plane)
    knee_direction = sign(cross_knee(:,2));
    ankle_direction = sign(cross_ankle(:,2));
    
    % Calculate dot products and norms for all frames
    knee_dots = sum(-thigh_vec .* shank_vec, 2);
    ankle_dots = sum(shank_vec .* foot_vec, 2);
    
    thigh_norms = sqrt(sum(thigh_vec.^2, 2));
    shank_norms = sqrt(sum(shank_vec.^2, 2));
    foot_norms = sqrt(sum(foot_vec.^2, 2));
    
    % Calculate angles
    knee_angles = knee_direction .* acosd(knee_dots ./ (thigh_norms .* shank_norms));
    ankle_angles = ankle_direction .* acosd(ankle_dots ./ (shank_norms .* foot_norms));
    
    angles = [knee_angles, ankle_angles];
end