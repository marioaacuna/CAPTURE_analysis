function ML_features = compute_wl_transform_features_demo_matfile(coeffstruct_in, overwrite_coeff,MLmatobjfile)
% TRUE Memory-optimized version - loads only specific fields when needed
global GC

fprintf('Starting wavelet transform computation with TRUE memory optimization...\n');

% Load only the metadata we need
fprintf('Loading only appendage_anglegps from matfile...\n');
appendage_anglegps = load(ML_matobjfile, 'appendage_anglegps');
if isstruct(appendage_anglegps)
    appendage_anglegps = appendage_anglegps.appendage_anglegps;
end
%appendage_anglegps = ml_matfile.appendage_anglegps;

if exist(coeffstruct_in,'file')
    try
        coeffstruct = load(coeffstruct_in);
    catch ME
        coeffstruct = load(coeffstruct_in);
    end
else
    coeffstruct = struct();
end

%% spectrogram parameters
opts.fps = GC.upsampling_to;%300./1;

opts.clustering_window = opts.fps./2;
opts.clustering_overlap = opts.fps./4;
opts.numclusters = 100;
opts.lambda = 0.1;
opts.num = 1;

opts.whiten = 0;
opts.frameNormalize = 0;
opts.clustermethod = 'GMM';
opts.ds = 1;
opts.samprate = GC.upsampling_to;
opts.params = struct;
opts.params.samplingFreq = GC.upsampling_to;
opts.params.numPeriods = 25;
opts.params.minF = 1;
opts.params.maxF = 25;
spacing = 6;
hipass_val = 0.5;

opts.pcuse = 20;
opts.numclusters = 100;
opts.lambda = 0.1;
num_spectrogram_pcs = 15;

% params.fps = 300;
params.fps = GC.upsampling_to;
params.difforder = 10;
params.medfiltorder = 3;
params.gaussorder = 2.5;

dHipass = designfilt('highpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', hipass_val/(params.fps/2), ...
    'DesignMethod', 'butter');
[f1_hipass,f2_hipass] = tf(dHipass);

for kk = 8
    fprintf('starting group %i joint angles \n',kk)

    % Load only the specific data we need for this iteration
    field_base = sprintf('appendage_%d_', kk);
    fprintf('Loading only %s from matfile...\n', strcat(field_base, 'joint_angles_pcs_hipass'));
    appendage_joint_angles_pcs_hipass_kk = ml_matfile.(strcat(field_base, 'joint_angles_pcs_hipass'));

    if size(appendage_joint_angles_pcs_hipass_kk,1) > 10*1

        %% JOINT ANGLES WAVELETS - SINGLE PRECISION
        agg_features_wl = single([]);
        num_components = size(appendage_joint_angles_pcs_hipass_kk,2);

        for ll = 1:num_components
            % Clean data - force single precision
            frame_data = single(appendage_joint_angles_pcs_hipass_kk(:,ll));
            frame_data(isnan(frame_data) | isinf(frame_data)) = 0;

            % Filter and process
            frame_fragment = single(filtfilt(f1_hipass, f2_hipass, frame_data));
            frame_fragment(isnan(frame_fragment) | isinf(frame_fragment)) = 0;
            frame_fragment = single(real(frame_fragment));

            opts.samprate = opts.fps./spacing;
            [~, w_map, ~] = return_wavelets(frame_fragment(1:spacing:end,1), 1:size(frame_fragment(1:spacing:end,1),1), opts);

            w_map = single(w_map + 3);
            w_map(w_map < 0) = 0;
            agg_features_wl = single(cat(2, agg_features_wl, w_map));

            % Clear large temporary variables immediately
            clear frame_data frame_fragment w_map
        end

        %% Load/compute coefficients - SINGLE PRECISION
        coeffname_appendage = strcat('COEFFS_appendages_wl',num2str(kk));
        explainedname_appendage = strcat('EXPLAINED_appendages_wl',num2str(kk));

        if (~isfield(coeffstruct,coeffname_appendage) || overwrite_coeff)
            fprintf('OVERWRITING APPENDAGES WL \n')
                        % Memory-optimized PCA on wavelet features
            [COEFFS_feat_wl_appendages_kk, ~, ~, ~, explained_wl_appendages_kk] = pca(agg_features_wl(1:1:end,:), 'Economy', true);
            COEFFS_feat_wl_appendages_kk = single(COEFFS_feat_wl_appendages_kk);
            explained_wl_appendages_kk = single(explained_wl_appendages_kk);
            coeffstruct.(coeffname_appendage) = COEFFS_feat_wl_appendages_kk;
            coeffstruct.(explainedname_appendage) = explained_wl_appendages_kk;
        else
            COEFFS_feat_wl_appendages_kk = single(coeffstruct.(coeffname_appendage));
            explained_wl_appendages_kk = single(coeffstruct.(explainedname_appendage));
        end

        % Center data and compute scores - SINGLE PRECISION
        agg_mean = single(mean(agg_features_wl, 1));
        agg_features_wl = single(bsxfun(@minus, agg_features_wl, agg_mean));
        dyadic_spectrograms_score_wl_appendages_kk = single(agg_features_wl * COEFFS_feat_wl_appendages_kk);

        % Clear large array immediately
        clear agg_features_wl agg_mean

        %% Replicate elements - SINGLE PRECISION
        replication_factor_wl = spacing;
        num_pcs_to_keep = min(num_spectrogram_pcs, size(dyadic_spectrograms_score_wl_appendages_kk,2));
        dyadic_spectrograms_score_wl_appendages_kk = single(repelem( ...
            dyadic_spectrograms_score_wl_appendages_kk(:,1:num_pcs_to_keep), replication_factor_wl, 1));

        %% Size adjustment - SINGLE PRECISION
        target_size = size(appendage_joint_angles_pcs_hipass_kk, 1);
        current_size = size(dyadic_spectrograms_score_wl_appendages_kk, 1);

        if current_size < target_size
            padding = single(zeros(target_size - current_size, size(dyadic_spectrograms_score_wl_appendages_kk,2)));
            dyadic_spectrograms_score_wl_appendages_kk = single(cat(1, dyadic_spectrograms_score_wl_appendages_kk, padding));
            clear padding
        else
            dyadic_spectrograms_score_wl_appendages_kk((end-(current_size-target_size)):end,:) = [];
        end

        % Add final row
        final_row = single(zeros(1, size(dyadic_spectrograms_score_wl_appendages_kk,2)));
        dyadic_spectrograms_score_wl_appendages_kk = single(cat(1, dyadic_spectrograms_score_wl_appendages_kk, final_row));
        clear final_row

        %% EUCLIDEAN WAVELETS - SINGLE PRECISION (same pattern)
        fprintf('starting group %f euclidean \n',kk)

        % Load euclidean data from matfile and clear the joint angles data
        clear appendage_joint_angles_pcs_hipass_kk  % Free memory
        fprintf('Loading only %s from matfile...\n', strcat(field_base, 'pca_score_euc_hipassclip'));
        appendage_pca_score_euc_hipassclip_kk = ml_matfile.(strcat(field_base, 'pca_score_euc_hipassclip'));

        agg_features_wl_euc = single([]);

        for ll = 1:size(appendage_pca_score_euc_hipassclip_kk,2)
            % Clean data
            frame_data = single(appendage_pca_score_euc_hipassclip_kk(:,ll));
            frame_data(isnan(frame_data) | isinf(frame_data)) = 0;

            frame_fragment = single(filtfilt(f1_hipass, f2_hipass, frame_data));
            frame_fragment(isnan(frame_fragment) | isinf(frame_fragment)) = 0;
            frame_fragment = single(real(frame_fragment));

            opts.samprate = opts.fps./spacing;
            [~, w_map, ~] = return_wavelets(frame_fragment(1:6:end,1), 1:size(frame_fragment(1:6:end,1),1), opts);

            w_map = single(w_map + 3);
            w_map(w_map < 0) = 0;
            agg_features_wl_euc = single(cat(2, agg_features_wl_euc, w_map));

            clear frame_data frame_fragment w_map
        end

        %% Coefficients for euclidean
        coeffname_appendage = strcat('COEFFS_appendages_wl_euc',num2str(kk));
        explainedname_appendage = strcat('EXPLAINED_appendages_wl_euc',num2str(kk));

        if (~isfield(coeffstruct,coeffname_appendage) || overwrite_coeff)
            fprintf('OVERWRITING APPENDAGES WL EUC \n')
            [COEFFS_feat_wl_appendages_euc_kk, ~, ~, ~, explained_wl_appendages_euc_kk] = pca(agg_features_wl_euc(1:1:end,:), 'Economy', true);
            COEFFS_feat_wl_appendages_euc_kk = single(COEFFS_feat_wl_appendages_euc_kk);
            explained_wl_appendages_euc_kk = single(explained_wl_appendages_euc_kk);
            coeffstruct.(coeffname_appendage) = COEFFS_feat_wl_appendages_euc_kk;
            coeffstruct.(explainedname_appendage) = explained_wl_appendages_euc_kk;
        else
            COEFFS_feat_wl_appendages_euc_kk = single(coeffstruct.(coeffname_appendage));
            explained_wl_appendages_euc_kk = single(coeffstruct.(explainedname_appendage));
        end

        agg_mean_euc = single(mean(agg_features_wl_euc, 1));
        agg_features_wl_euc = single(bsxfun(@minus, agg_features_wl_euc, agg_mean_euc));
        dyadic_spectrograms_score_wl_appendages_euc_kk = single(agg_features_wl_euc * COEFFS_feat_wl_appendages_euc_kk);

        clear agg_features_wl_euc agg_mean_euc

        %% Replicate and adjust size for euclidean
        num_pcs_to_keep_euc = min(num_spectrogram_pcs, size(dyadic_spectrograms_score_wl_appendages_euc_kk,2));
        dyadic_spectrograms_score_wl_appendages_euc_kk = single(repelem( ...
            dyadic_spectrograms_score_wl_appendages_euc_kk(:,1:num_pcs_to_keep_euc), replication_factor_wl, 1));

        current_size_euc = size(dyadic_spectrograms_score_wl_appendages_euc_kk, 1);
        if current_size_euc < target_size
            padding_euc = single(zeros(target_size - current_size_euc, size(dyadic_spectrograms_score_wl_appendages_euc_kk,2)));
            dyadic_spectrograms_score_wl_appendages_euc_kk = single(cat(1, dyadic_spectrograms_score_wl_appendages_euc_kk, padding_euc));
            clear padding_euc
        else
            dyadic_spectrograms_score_wl_appendages_euc_kk((end-(current_size_euc-target_size)):end,:) = [];
        end

        final_row_euc = single(zeros(1, size(dyadic_spectrograms_score_wl_appendages_euc_kk,2)));
        dyadic_spectrograms_score_wl_appendages_euc_kk = single(cat(1, dyadic_spectrograms_score_wl_appendages_euc_kk, final_row_euc));
        clear final_row_euc

        % Clear euclidean data from memory
        clear appendage_pca_score_euc_hipassclip_kk

        %% Store results back to matfile - SINGLE PRECISION
        ml_matfile.ML_features.COEFFS_feat_wl_appendages_euc(1,kk) = {COEFFS_feat_wl_appendages_euc_kk};
        ml_matfile.ML_features.dyadic_spectrograms_score_wl_appendages_euc(1,kk) = {dyadic_spectrograms_score_wl_appendages_euc_kk};
        ml_matfile.ML_features.explained_wl_appendages_euc(1,kk) = {explained_wl_appendages_euc_kk};

        ml_matfile.ML_features.COEFFS_feat_wl_appendages(1,kk) = {COEFFS_feat_wl_appendages_kk};
        ml_matfile.ML_features.dyadic_spectrograms_score_wl_appendages(1,kk) = {dyadic_spectrograms_score_wl_appendages_kk};
        ml_matfile.ML_features.explained_wl_appendages(1,kk) = {explained_wl_appendages_kk};

        % Clear processed variables
        clear COEFFS_feat_wl_appendages_kk explained_wl_appendages_kk dyadic_spectrograms_score_wl_appendages_kk
        clear COEFFS_feat_wl_appendages_euc_kk explained_wl_appendages_euc_kk dyadic_spectrograms_score_wl_appendages_euc_kk

    else
        %% Empty case - store empty arrays to matfile
        fprintf('Storing empty arrays for group %d...\n', kk);

        ml_matfile.ML_features.COEFFS_feat_wl_appendages_euc(1,kk) = {[]};
        ml_matfile.ML_features.dyadic_spectrograms_score_wl_appendages_euc(1,kk) = {[]};
        ml_matfile.ML_features.explained_wl_appendages_euc(1,kk) = {[]};

        ml_matfile.ML_features.COEFFS_feat_wl_appendages(1,kk) = {[]};
        ml_matfile.ML_features.dyadic_spectrograms_score_wl_appendages(1,kk) = {[]};
        ml_matfile.ML_features.explained_wl_appendages(1,kk) = {[]};
    end
end

%% Save coefficients
fprintf('saving appendage coefficients !!! WAVELET !!!\n')
try
    save(coeffstruct_in,'-struct','coeffstruct','-v7.3')
catch ME
    save(coeffstruct_in,'-struct','coeffstruct','-v7.3')
end


% Return empty ML_features since everything is stored in matfile
ML_features = struct();
ML_features.status = 'completed_wavelet_processing';

fprintf('Wavelet transform computation completed with matfile optimization.\n');

end
