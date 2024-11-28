function anglestruct = load_mario_mouse22_anglestruct_v2()

%% XYZ reference frame angles
% Sagittal angles (side view, z-y plane)
anglestruct.saggital_names = {'head_sagg', 'neck_sagg', 'spine_sagg', 'tail_sagg'};
anglestruct.saggital_pairs = {[2,3], [3,4], [4,5], [5,22]};

% Transverse angles (overhead view, x-y plane)
anglestruct.transverse_names = {'head_trans', 'neck_trans', 'spine_trans', 'KneeL_trans', 'KneeR_trans', 'shoulderL_trans', 'shoulderR_trans'};
anglestruct.transverse_pairs = {[2,3], [3,4], [4,5], [5,6], [5,7], [4,8], [4,9]};

% Coronal angles (front view, x-z plane)
anglestruct.coronal_names = {'head_coronal', 'KneeL_coronal', 'KneeR_coronal', 'shoulderL_coronal', 'shoulderR_coronal'};
anglestruct.coronal_pairs = {[3,4], [5,18], [5,19], [4,8], [4,9]};

% % Coronal angles with elbow
anglestruct.coronal_names = {'head_coronal', 'KneeL_coronal', 'KneeR_coronal', 'shoulderL_coronal', 'shoulderR_coronal', 'elbowL_coronal', 'elbowR_coronal'};
anglestruct.coronal_pairs = {[3,4], [5,18], [5,19], [4,8], [4,9], [10,12], [11,13]};

%% Segment pairs for all angles
anglestruct.segment_pairs = {
    {'EarL', 'EarR'},   ...     % 1
    {'EarR', 'Snout'},     ...  % 2
    {'Snout', 'SpineF'},     ...% 3
    {'SpineF', 'SpineM'},...    % 4
    {'SpineM', 'Tail_base_'},... % 5
    {'Tail_base_', 'KneeL'},... % 6
    {'Tail_base_', 'KneeR'},... % 7
    {'SpineF', 'ShoulderL'},... % 8
    {'SpineF', 'ShoulderR'},... % 9
    {'ShoulderL', 'ElbowL'},... % 10
    {'ShoulderR', 'ElbowR'},... % 11
    {'ElbowL', 'WristL'},...    % 12
    {'ElbowR', 'WristR'},...    % 13
    {'WristL', 'ForepawL'},...  % 14
    {'WristR', 'ForepawR'},...  % 15
    {'ShoulderL', 'ElbowL'},... % 16 (duplicate of 10, but needed for indexing)
    {'ShoulderR', 'ElbowR'},... % 17 (duplicate of 11, but needed for indexing)
    {'KneeL', 'AnkleL'},...     % 18
    {'KneeR', 'AnkleR'},...     % 19
    {'AnkleL', 'HindpawL'},...  % 20
    {'AnkleR', 'HindpawR'},...  % 21
    {'Tail_base_', 'Tail_mid_'} % 22
};
%% Planar trios for complex angle calculations
% Hip angles
anglestruct.planar_trios{1}.plane = {{'SpineM', 'Tail_base_'}, {'zvector', 'zvector'}};
anglestruct.planar_trios{1}.vector = {'Tail_base_', 'KneeL'};
anglestruct.planar_trios{1}.name1 = 'hipl_pitch';
anglestruct.planar_trios{1}.name2 = 'hipl_yaw';
anglestruct.planar_trios{1}.namesuse = [1 2];

anglestruct.planar_trios{2}.plane = {{'SpineM', 'Tail_base_'}, {'zvector', 'zvector'}};
anglestruct.planar_trios{2}.vector = {'Tail_base_', 'KneeR'};
anglestruct.planar_trios{2}.name1 = 'hipr_pitch';
anglestruct.planar_trios{2}.name2 = 'hipr_yaw';
anglestruct.planar_trios{2}.namesuse = [1 2];

% Knee angles
anglestruct.planar_trios{3}.plane = {{'Tail_base_', 'SpineM'}, {'KneeR', 'KneeL'}};
anglestruct.planar_trios{3}.vector = {'KneeL', 'AnkleL'};
anglestruct.planar_trios{3}.name1 = 'kneel_yaw';
anglestruct.planar_trios{3}.name2 = 'kneel_pitch';
anglestruct.planar_trios{3}.namesuse = [1 2];

anglestruct.planar_trios{4}.plane = {{'Tail_base_', 'SpineM'}, {'KneeL', 'KneeR'}};
anglestruct.planar_trios{4}.vector = {'KneeR', 'AnkleR'};
anglestruct.planar_trios{4}.name1 = 'kneer_yaw';
anglestruct.planar_trios{4}.name2 = 'kneer_pitch';
anglestruct.planar_trios{4}.namesuse = [1 2];

% Ankle angles (formerly shin angles)
anglestruct.planar_trios{5}.plane = {{'KneeR', 'KneeL'}, {'KneeL', 'AnkleL'}};
anglestruct.planar_trios{5}.vector = {'HindpawL', 'AnkleL'};
anglestruct.planar_trios{5}.name1 = 'anklel_yaw';
anglestruct.planar_trios{5}.name2 = 'anklel_pitch';
anglestruct.planar_trios{5}.namesuse = [1 2];

anglestruct.planar_trios{6}.plane = {{'KneeL', 'KneeR'}, {'KneeR', 'AnkleR'}};
anglestruct.planar_trios{6}.vector = {'HindpawR', 'AnkleR'};
anglestruct.planar_trios{6}.name1 = 'ankler_yaw';
anglestruct.planar_trios{6}.name2 = 'ankler_pitch';
anglestruct.planar_trios{6}.namesuse = [1 2];

% Shoulder angles (formerly arm angles)
anglestruct.planar_trios{7}.plane = {{'SpineF', 'SpineM'}, {'ShoulderR', 'ShoulderL'}};
anglestruct.planar_trios{7}.vector = {'ShoulderL', 'ElbowL'};
anglestruct.planar_trios{7}.name1 = 'shoulderl_yaw';
anglestruct.planar_trios{7}.name2 = 'shoulderl_pitch';
anglestruct.planar_trios{7}.namesuse = [1 2];

anglestruct.planar_trios{8}.plane = {{'SpineF', 'SpineM'}, {'ShoulderL', 'ShoulderR'}};
anglestruct.planar_trios{8}.vector = {'ShoulderR', 'ElbowR'};
anglestruct.planar_trios{8}.name1 = 'shoulderr_yaw';
anglestruct.planar_trios{8}.name2 = 'shoulderr_pitch';
anglestruct.planar_trios{8}.namesuse = [1 2];

% Elbow angles (revised as suggested)
anglestruct.planar_trios{9}.plane = {{'ShoulderL', 'ElbowL'}, {'ShoulderR', 'ShoulderL'}};
anglestruct.planar_trios{9}.vector = {'ElbowL', 'WristL'};
anglestruct.planar_trios{9}.name1 = 'elbowl_pitch';
anglestruct.planar_trios{9}.name2 = 'elbowl_yaw';
anglestruct.planar_trios{9}.namesuse = [1 2];

anglestruct.planar_trios{10}.plane = {{'ShoulderR', 'ElbowR'}, {'ShoulderL', 'ShoulderR'}};
anglestruct.planar_trios{10}.vector = {'ElbowR', 'WristR'};
anglestruct.planar_trios{10}.name1 = 'elbowr_pitch';
anglestruct.planar_trios{10}.name2 = 'elbowr_yaw';
anglestruct.planar_trios{10}.namesuse = [1 2];

%% Angles to include in analysis
anglestruct.include_angles = [
    anglestruct.saggital_names, anglestruct.transverse_names, anglestruct.coronal_names, ...
    'hipl_pitch', 'hipl_yaw', 'hipr_pitch', 'hipr_yaw', ...
    'kneel_pitch', 'kneel_yaw', 'kneer_pitch', 'kneer_yaw', ...
    'anklel_pitch', 'anklel_yaw', 'ankler_pitch', 'ankler_yaw', ...
    'shoulderl_pitch', 'shoulderl_yaw', 'shoulderr_pitch', 'shoulderr_yaw', ...
    'elbowl_pitch', 'elbowl_yaw', 'elbowr_pitch', 'elbowr_yaw'
];

end