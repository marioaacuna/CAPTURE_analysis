function ML_features = create_behavioral_features(mocapstruct,coeff_file,overwrite_coeff, linkname, MLmatobjfile)
% Make features for tsne - Memory Optimized Version
% ---------------------------
% (C) Jesse D Marshall 2020
%     Harvard University
% Modified for memory optimization using matfile objects

fprintf('Starting memory-optimized behavioral feature computation...\n');

global GC
%% compute the joint angles - returns matfile-based structure
compute_joint_angles_demo(mocapstruct,linkname, MLmatobjfile);

% Check if we got a matfile-based structure
%if isfield(ML_features, 'is_matfile_based') && ML_features.is_matfile_based
    fprintf('Using matfile-based ML_features for memory optimization.\n');
    % matfile_path = ML_features.matfile_path;

    % % Create matfile object for the pipeline
    % ml_matfile = matfile(matfile_path, 'Writable', true);

    %% compute the principal components of the joint angles
    fprintf('Computing appendage PCs with matfile optimization...\n');
    ML_features = compute_appendage_pc_demos_matfile(mocapstruct, coeff_file, overwrite_coeff, MLmatobjfile);

    %% compute the wavelet transform
    fprintf('Computing wavelet transform with matfile optimization...\n');
    tic
    ML_features = compute_wl_transform_features_demo_matfile_new(coeff_file, overwrite_coeff,MLmatobjfile);
    toc

    %% Return matfile-based ML_features structure for downstream use
    fprintf('Returning matfile-based ML_features structure...\n');
    %ML_features.is_matfile_based = true;
    %ML_features.matfile_path = matfile_path;
    %ML_features.matfile_obj = ml_matfile;  % Keep the matfile object available
    % try
    %     ML_features.mocapstruct = ml_matfile.mocapstruct;
    % catch
    %     warning('Could not load mocapstruct from matfile');
    % end
    ML_features = load(MLmatobjfile);
    %fprintf('Matfile-based ML_features ready. Use matfile_path or matfile_obj to access data.\n');
    %fprintf('Matfile location: %s\n', matfile_path);
    %fprintf('Keeping matfile for downstream processing: %s\n', matfile_path);

%else
%    % Fallback to original method if matfile approach failed
%    warning('Falling back to original memory-intensive approach.');
%
%    %% compute the principal components of the joint angles and
%    ML_features = compute_appendage_pc_demos(mocapstruct,ML_features,coeff_file,overwrite_coeff);
%    %% compute the wavelet transform
%    tic
%    ML_features = compute_wl_transform_features_demo(mocapstruct,ML_features,coeff_file,overwrite_coeff);
%    toc
%end

%% add window/vel/old features
% TODO: Add additional feature computation here

fprintf('Behavioral feature computation completed.\n');

%% code for visualization
