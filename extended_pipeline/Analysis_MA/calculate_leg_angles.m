
function leg_angles = calculate_leg_angles(spineM, spineF, knee, ankle, paw)
    % Calculate angles between markers in 3D space
    % spineM, spineF, knee, ankle, paw are Nx3 matrices where N is the number of frames

    % Calculate vectors
    thigh_vector = knee - spineM;
    shank_vector = ankle - knee;
    foot_vector = paw - ankle;
    hip_vector = spineF - spineM;

    % Calculate angles
    knee_angle = acosd(dot(thigh_vector, shank_vector, 2) ./ (vecnorm(thigh_vector, 2, 2) .* vecnorm(shank_vector, 2, 2)));
    ankle_angle = acosd(dot(shank_vector, foot_vector, 2) ./ (vecnorm(shank_vector, 2, 2) .* vecnorm(foot_vector, 2, 2)));
    hip_angle = acosd(dot(hip_vector, thigh_vector, 2) ./ (vecnorm(hip_vector, 2, 2) .* vecnorm(thigh_vector, 2, 2)));

    % Combine angles into a single matrix
    leg_angles = [knee_angle, ankle_angle, hip_angle];
end