function ML_features = compute_joint_angles_demo(mocapstruct,linkname)
ML_features = struct();

%% get joint angle features
fprintf('%% computing joint angles %% \n');

%% Define indices for different planes
saggital_inds = [2,3];
coronal_inds = [1,3];
transverse_inds = [1,2];
allangles_inds = [1,2,3];

% Load angle structure based on linkname
switch linkname
    case 'rats'
        anglestruct = load_default_anglestruct();
    case 'bird'
        anglestruct = load_bird_anglestruct();
    case 'mouse'
        anglestruct = load_mouse_anglestruct();
    case 'kylemouse'
        anglestruct = load_mouse_kyle_anglestruct();
    case {'mario_mouse', 'mario_mouse_14_pts'}
        anglestruct = load_mario_mouse_anglestruct();
    case 'mario_mouse22'
        anglestruct = load_mario_mouse22_anglestruct_v2();
end

segment_pairs = anglestruct.segment_pairs;
coronal_pairs = anglestruct.coronal_pairs;
saggital_pairs = anglestruct.saggital_pairs;
transverse_pairs = anglestruct.transverse_pairs;
planar_trios = anglestruct.planar_trios;

%% Load variables from anglestruct
assigns = structvars(anglestruct);
for kk = 1:size(assigns,1)
    eval(assigns(kk,:))
end

%% Pre-allocate structures - SINGLE PRECISION
jointangle_struct = struct();
all_seglengths = cell(1,numel(segment_pairs));
all_segments = cell(1,numel(segment_pairs));
transverse_seglengths = cell(1,numel(segment_pairs));

ML_features.segment_pairs = segment_pairs;
ML_features.include_angles = include_angles;

%% SAGGITAL ANGLES - SINGLE PRECISION
for ll = 1:numel(saggital_pairs)
    % Extract vectors as single precision
    vec1 = single(mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(1)}{1}) - ...
                  mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(1)}{2}));
    
    vec2 = single(mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(2)}{1}) - ...
                  mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(2)}{2}));
    
    % Compute angles with single precision
    vec1_proj = single(vec1(:,saggital_inds));
    vec2_proj = single(vec2(:,saggital_inds));
    
    dot_prod = single(sum(vec1_proj .* vec2_proj, 2));
    norm1 = single(sqrt(sum(vec1_proj.^2, 2)));
    norm2 = single(sqrt(sum(vec2_proj.^2, 2)));
    
    jointangle_struct.(saggital_names{ll}) = single(acosd(dot_prod ./ (norm1 .* norm2)));
    
    % Clear temporary variables immediately
    clear vec1 vec2 vec1_proj vec2_proj dot_prod norm1 norm2
end

%% CORONAL ANGLES - SINGLE PRECISION
for ll = 1:numel(coronal_pairs)
    vec1 = single(mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(1)}{1}) - ...
                  mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(1)}{2}));
    
    vec2 = single(mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(2)}{1}) - ...
                  mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(2)}{2}));
    
    vec1_proj = single(vec1(:,coronal_inds));
    vec2_proj = single(vec2(:,coronal_inds));
    
    dot_prod = single(sum(vec1_proj .* vec2_proj, 2));
    norm1 = single(sqrt(sum(vec1_proj.^2, 2)));
    norm2 = single(sqrt(sum(vec2_proj.^2, 2)));
    
    jointangle_struct.(coronal_names{ll}) = single(acosd(dot_prod ./ (norm1 .* norm2)));
    
    clear vec1 vec2 vec1_proj vec2_proj dot_prod norm1 norm2
end

%% TRANSVERSE ANGLES - SINGLE PRECISION
for ll = 1:numel(transverse_pairs)
    vec1 = single(mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(1)}{1}) - ...
                  mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(1)}{2}));
    
    vec2 = single(mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(2)}{1}) - ...
                  mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(2)}{2}));
    
    vec1_proj = single(vec1(:,transverse_inds));
    vec2_proj = single(vec2(:,transverse_inds));
    
    dot_prod = single(sum(vec1_proj .* vec2_proj, 2));
    norm1 = single(sqrt(sum(vec1_proj.^2, 2)));
    norm2 = single(sqrt(sum(vec2_proj.^2, 2)));
    
    jointangle_struct.(transverse_names{ll}) = single(acosd(dot_prod ./ (norm1 .* norm2)));
    
    clear vec1 vec2 vec1_proj vec2_proj dot_prod norm1 norm2
end

%% TRANSVERSE SEGMENT LENGTHS - SINGLE PRECISION
for ll = 1:numel(transverse_pairs)
    vec1 = single(mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(1)}{1}) - ...
                  mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(1)}{2}));
    
    vec1_proj = single(vec1(:,transverse_inds));
    transverse_seglengths{ll} = single(sqrt(sum(vec1_proj.^2, 2)));
    
    clear vec1 vec1_proj
end

%% PLANAR ANGLES - SINGLE PRECISION
for kk = 1:numel(planar_trios)
    if numel(planar_trios{kk}.namesuse) == 2
        [angle1, angle2] = get_planar_jointangles(mocapstruct, planar_trios{kk}.plane, planar_trios{kk}.vector);
        jointangle_struct.(planar_trios{kk}.name1) = single(angle1);
        jointangle_struct.(planar_trios{kk}.name2) = single(angle2);
        clear angle1 angle2
    elseif find(planar_trios{kk}.namesuse == 1)
        [angle1, ~] = get_planar_jointangles(mocapstruct, planar_trios{kk}.plane, planar_trios{kk}.vector);
        jointangle_struct.(planar_trios{kk}.name1) = single(angle1);
        clear angle1
    elseif find(planar_trios{kk}.namesuse == 2)
        [~, angle2] = get_planar_jointangles(mocapstruct, planar_trios{kk}.plane, planar_trios{kk}.vector);
        jointangle_struct.(planar_trios{kk}.name2) = single(angle2);
        clear angle2
    end
end

%% SEGMENT LENGTHS AND VECTORS - SINGLE PRECISION
fprintf('getting segment lengths \n')
for ll = 1:numel(segment_pairs)
    vec1 = single(mocapstruct.markers_aligned_preproc.(segment_pairs{ll}{1}) - ...
                  mocapstruct.markers_aligned_preproc.(segment_pairs{ll}{2}));
    
    vec1_all = single(vec1(:,allangles_inds));
    all_seglengths{ll} = single(sqrt(sum(vec1_all.^2, 2)));
    all_segments{ll} = single(vec1);
    
    clear vec1 vec1_all
end

%% CLEAN UP JOINT ANGLES - SINGLE PRECISION
jointangle_struct = structfun(@(x) single(real(x)), jointangle_struct, 'UniformOutput', false);
fname = fieldnames(jointangle_struct);
for lk = 1:numel(fname)
    temp_angles = jointangle_struct.(fname{lk});
    temp_angles(isnan(temp_angles)) = 0;
    jointangle_struct.(fname{lk}) = single(temp_angles);
    clear temp_angles
end

%% COMPUTE MEANS AND STORE RESULTS - SINGLE PRECISION
ML_features.joint_angles_mean = single(real(structfun(@nanmean, jointangle_struct)));
ML_features.jointangle_struct = jointangle_struct;
ML_features.all_seglengths = all_seglengths;
ML_features.all_segments = all_segments;
ML_features.transverse_seglengths = transverse_seglengths;

% Final cleanup
clear anglestruct segment_pairs coronal_pairs saggital_pairs transverse_pairs planar_trios

end