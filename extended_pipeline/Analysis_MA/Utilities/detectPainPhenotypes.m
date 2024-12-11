function [pain_frames, metrics] = detectPainPhenotypes(non_aligned_mocap, aligned_mocap, params)
    % Input: both non_aligned and aligned mocap data structures
    % Output: pain_frames - logical array indicating frames with pain behavior
    %         metrics - structure with detailed measurements
    
    % Step 1: Identify still frames based on non-aligned SpineM velocity
    vel_threshold = 0.05; % adjust based on your data
    spine_vel = calculateVelocity(non_aligned_mocap.SpineM, params);
    still_frames = spine_vel < median(spine_vel);  

    % Initialize output
    pain_frames = false(size(still_frames));

    num_frames =size(still_frames,1);
    % Initialize metrics structure with arrays for each measurement
    metrics = struct(...
        'snout_paw_dist', zeros(num_frames, 1), ...
        'L_hind_fore_paw_dist', zeros(num_frames, 1), ...
        'trunk_angle_coronal', zeros(num_frames, 1), ...
        'lateral_asymmetry', zeros(num_frames, 1), ...
        'leg_asymmetry', zeros(num_frames, 1), ...
        'arm_asymmetry', zeros(num_frames, 1), ...
        'neck_angle', zeros(num_frames, 1), ...
        'right_paw_height', zeros(num_frames, 1), ...
        'left_paw_height', zeros(num_frames, 1), ...
        'is_still', false(num_frames, 1), ...
        'paw_licking_detected', false(num_frames, 1), ...
        'timestamp', zeros(num_frames, 1) ...
    );


    % Only analyze still frames using aligned data
    for frame = find(still_frames)'
        % Calculate key metrics for each frame using aligned data
        
        % 1. Paw Licking Detection
        % Distance between snout and left hindpaw in aligned space
        snout_paw_dist = norm(aligned_mocap.Snout(frame,:) - aligned_mocap.HindpawL(frame,:));
        L_hind_fore_paw_dist = norm(aligned_mocap.HindpawL(frame,:) - aligned_mocap.ForepawL(frame,:));
        
        % Calculate trunk angles in coronal plane (XY plane in aligned data)
        % Since data is centered on SpineM, we can directly use relative positions
        spineF_rel = aligned_mocap.SpineF(frame,:) - aligned_mocap.SpineM(frame,:);
        tailbase_rel = aligned_mocap.Tail_base_(frame,:) - aligned_mocap.SpineM(frame,:);
        
        % Trunk angle in coronal plane (using X and Y components)
        trunk_angle_coronal = calculateCoronalAngle(aligned_mocap.Snout(frame,2),spineF_rel, tailbase_rel);
        
        % Check right paw floor contact using aligned Z coordinate
        right_paw_height = aligned_mocap.HindpawR(frame,3);
        floor_threshold = 0.02; % adjust based on your data
        
        % 2. Posture Analysis in aligned space
        % All vectors are already relative to SpineM
        neck_angle = calculateNeckAngle(aligned_mocap.Snout(frame,:), ...
                                      aligned_mocap.SpineF(frame,:), ...
                                      aligned_mocap.SpineM(frame,:));
                                      
        [left_leg_angle, right_leg_angle] = calculateLegAngles(aligned_mocap, frame);
        [left_arm_angle, right_arm_angle] = calculateArmAngles(aligned_mocap, frame);
        
        % Calculate lateral deviation in coronal plane
        left_paw_lateral = aligned_mocap.HindpawL(frame,2); % Y coordinate
        right_paw_lateral = aligned_mocap.HindpawR(frame,2);
        lateral_asymmetry = left_paw_lateral - right_paw_lateral;
        
        % Define pain criteria using aligned metrics
        is_paw_licking =   snout_paw_dist < 10 && ...
                           aligned_mocap.Snout(frame,2) > 10 &&...; % leaning left
                           L_hind_fore_paw_dist < 15 &&...
                           aligned_mocap.HindpawL(frame,1) > 0;
                           
                        % trunk_angle_coronal < 90 ; ... % negative angle indicates left lean
                        % right_paw_height < floor_threshold;
        % if is_paw_licking, keyboard, end
        % Posture asymmetry check
        leg_asymmetry = abs(left_leg_angle - right_leg_angle);
        arm_asymmetry = abs(left_arm_angle - right_arm_angle);
        
        % Combine criteria for final pain detection
        pain_frames(frame) = is_paw_licking && ...
                            (leg_asymmetry > 10) && ... % adjust threshold
                            (lateral_asymmetry > 15) &&... % adjust threshold
                            trunk_angle_coronal > 15;

        % Store metrics for this frame
        metrics.snout_paw_dist(frame) = snout_paw_dist;
        metrics.L_hind_fore_paw_dist(frame) = L_hind_fore_paw_dist;

        metrics.trunk_angle_coronal(frame) = trunk_angle_coronal;
        metrics.lateral_asymmetry(frame) = lateral_asymmetry;
        metrics.leg_asymmetry(frame) = leg_asymmetry;
        metrics.arm_asymmetry(frame) = arm_asymmetry;
        metrics.neck_angle(frame) = neck_angle;
        metrics.paw_licking_detected(frame) = is_paw_licking;
    end
   
end

% Helper functions
function speed_smooth = calculateVelocity(positions, params)
    % Calculate velocity using central difference
    % XY velocity processing
    vel_xy = diff(positions(:,1:2)) * params.sampling_rate;
    vel_xy = [zeros(1,2); vel_xy];
    speed = sqrt(sum(vel_xy.^2, 2));
    speed_smooth = movmean(speed, params.smoothing_window);

end

function angle = calculateCoronalAngle(snouty,vector1, vector2)
    % Calculate angle in coronal (XY) plane
    % Project vectors onto XY plane
    v1_coronal = vector1(1:2);
    v2_coronal = vector2(1:2);

    % Calculate angle
    angle = atan2d(norm(cross([v1_coronal, 0], [v2_coronal, 0])), ...
                   dot(v1_coronal, v2_coronal));
    
    angle = 180-angle;
    % Determine sign (negative for left lean)
    if snouty < 0  % Y component negative indicates left lean
        angle = -angle;
    end
end

% function angle = calculateCoronalAngle(vector1, vector2)
%     % For spine vectors in YZ plane (coronal view)
%     % Y is medial-lateral 
%     % Z is dorsal-ventral
% 
%     % Project vectors onto YZ plane
%     v1_coronal = vector1([2 3]);  % [Y Z]
%     v2_coronal = vector2([2 3]);  % [Y Z]
% 
%     % Normalize vectors
%     v1_norm = v1_coronal / norm(v1_coronal);
%     v2_norm = v2_coronal / norm(v2_coronal);
% 
%     % Calculate angle using dot product and cross product
%     angle = atan2d(cross([v1_norm, 0], [v2_norm, 0]), dot(v1_norm, v2_norm));
% 
%     % Keep angle in [-90, 90] range
%     if angle > 90
%         angle = 180 - angle;
%     elseif angle < -90
%         angle = -180 - angle;
%     end
% end

function angle = calculateNeckAngle(snout, spine_f, spine_m)
    % Calculate neck angle using aligned coordinates
    neck_vector = snout - spine_f;
    spine_vector = spine_f - spine_m;
    
    % Project vectors onto sagittal (XZ) plane for primary angle
    neck_sagittal = neck_vector([1 3]);
    spine_sagittal = spine_vector([1 3]);
    
    angle = atan2d(norm(cross([neck_sagittal, 0], [spine_sagittal, 0])), ...
                   dot(neck_sagittal, spine_sagittal));
end

function [left_angle, right_angle] = calculateLegAngles(aligned_mocap, frame)
    % Calculate leg angles in aligned space
    left_angle = calculateJointAngle(aligned_mocap.KneeL(frame,:), ...
                                   aligned_mocap.AnkleL(frame,:), ...
                                   aligned_mocap.HindpawL(frame,:));
    right_angle = calculateJointAngle(aligned_mocap.KneeR(frame,:), ...
                                    aligned_mocap.AnkleR(frame,:), ...
                                    aligned_mocap.HindpawR(frame,:));
end

function [left_angle, right_angle] = calculateArmAngles(aligned_mocap, frame)
    % Calculate arm angles in aligned space
    left_angle = calculateJointAngle(aligned_mocap.ShoulderL(frame,:), ...
                                   aligned_mocap.ElbowL(frame,:), ...
                                   aligned_mocap.ForepawL(frame,:));
    right_angle = calculateJointAngle(aligned_mocap.ShoulderR(frame,:), ...
                                    aligned_mocap.ElbowR(frame,:), ...
                                    aligned_mocap.ForepawR(frame,:));
end

function angle = calculateJointAngle(p1, p2, p3)
    % Calculate angle between three points in aligned space
    v1 = p1 - p2;
    v2 = p3 - p2;
    angle = atan2d(norm(cross(v1, v2)), dot(v1, v2));
end