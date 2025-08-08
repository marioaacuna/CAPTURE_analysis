function ML_features = compute_wl_transform_features_demo(~,ML_features,coeffstruct_in,overwrite_coeff)
global GC

% Memory optimization: Save ML_features to temporary matfile and work with matfile object
fprintf('Saving ML_features to temporary file to reduce memory usage...\n');
temp_ml_file = [tempname, '_ML_features.mat'];
save(temp_ml_file, 'ML_features', '-v7.3');

% Create matfile object for efficient memory access
ml_matfile = matfile(temp_ml_file, 'Writable', true);

% Load only the minimal required fields into memory
appendage_anglegps = ml_matfile.ML_features.appendage_anglegps;
fprintf('ML_features saved to temp file. Working with matfile object to minimize memory usage.\n');

% Clear the large ML_features structure from memory
clear ML_features;

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

%% Pre-allocate with reasonable sizes - SINGLE PRECISION
COEFFS_feat_wl_appendages = cell(1,numel(appendage_anglegps));
dyadic_spectrograms_score_wl_appendages = cell(1,numel(appendage_anglegps));
explained_wl_appendages = cell(1,numel(appendage_anglegps));

COEFFS_feat_wl_appendages_euc = cell(1,numel(appendage_anglegps));
dyadic_spectrograms_score_wl_appendages_euc = cell(1,numel(appendage_anglegps));
explained_wl_appendages_euc = cell(1,numel(appendage_anglegps));

dHipass = designfilt('highpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', hipass_val/(params.fps/2), ...
    'DesignMethod', 'butter');
[f1_hipass,f2_hipass] = tf(dHipass);

for kk = 8
    fprintf('starting group %i joint angles \n',kk)
    
    % Load only the required data for this iteration from matfile
    fprintf('Loading appendage_joint_angles_pcs_hipass{%d} from matfile...\n', kk);
    appendage_joint_angles_pcs_hipass_kk = ml_matfile.ML_features.appendage_joint_angles_pcs_hipass(1,kk);
    appendage_joint_angles_pcs_hipass_kk = appendage_joint_angles_pcs_hipass_kk{1};
    
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
            [COEFFS_feat_wl_appendages{kk}, ~, ~, ~, explained_wl_appendages{kk}] = pca(single(squeeze(agg_features_wl)));
            COEFFS_feat_wl_appendages{kk} = single(COEFFS_feat_wl_appendages{kk});
            explained_wl_appendages{kk} = single(explained_wl_appendages{kk});
            coeffstruct.(coeffname_appendage) = COEFFS_feat_wl_appendages{kk};
            coeffstruct.(explainedname_appendage) = explained_wl_appendages{kk};
        else
            COEFFS_feat_wl_appendages{kk} = single(coeffstruct.(coeffname_appendage));
            explained_wl_appendages{kk} = single(coeffstruct.(explainedname_appendage));
        end
        
        % Center data and compute scores - SINGLE PRECISION
        agg_mean = single(mean(agg_features_wl, 1));
        agg_features_wl = single(bsxfun(@minus, agg_features_wl, agg_mean));
        dyadic_spectrograms_score_wl_appendages{kk} = single(agg_features_wl * COEFFS_feat_wl_appendages{kk});
        
        % Clear large array immediately
        clear agg_features_wl agg_mean
        
        %% Replicate elements - SINGLE PRECISION
        replication_factor_wl = spacing;
        num_pcs_to_keep = min(num_spectrogram_pcs, size(dyadic_spectrograms_score_wl_appendages{kk},2));
        dyadic_spectrograms_score_wl_appendages{kk} = single(repelem( ...
            dyadic_spectrograms_score_wl_appendages{kk}(:,1:num_pcs_to_keep), replication_factor_wl, 1));
        
        %% Size adjustment - SINGLE PRECISION
        target_size = size(appendage_joint_angles_pcs_hipass_kk, 1);
        current_size = size(dyadic_spectrograms_score_wl_appendages{kk}, 1);
        
        if current_size < target_size
            padding = single(zeros(target_size - current_size, size(dyadic_spectrograms_score_wl_appendages{kk},2)));
            dyadic_spectrograms_score_wl_appendages{kk} = single(cat(1, dyadic_spectrograms_score_wl_appendages{kk}, padding));
            clear padding
        else
            dyadic_spectrograms_score_wl_appendages{kk}((end-(current_size-target_size)):end,:) = [];
        end
        
        % Add final row
        final_row = single(zeros(1, size(dyadic_spectrograms_score_wl_appendages{kk},2)));
        dyadic_spectrograms_score_wl_appendages{kk} = single(cat(1, dyadic_spectrograms_score_wl_appendages{kk}, final_row));
        clear final_row
        
        %% EUCLIDEAN WAVELETS - SINGLE PRECISION (same pattern)
        fprintf('starting group %f euclidean \n',kk)
        
        % Load euclidean data from matfile and clear the joint angles data
        clear appendage_joint_angles_pcs_hipass_kk  % Free memory
        fprintf('Loading appendage_pca_score_euc_hipassclip{%d} from matfile...\n', kk);
        appendage_pca_score_euc_hipassclip_kk = ml_matfile.ML_features.appendage_pca_score_euc_hipassclip(1,kk);
        appendage_pca_score_euc_hipassclip_kk = appendage_pca_score_euc_hipassclip_kk{1};
        
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
            [COEFFS_feat_wl_appendages_euc{kk}, ~, ~, ~, explained_wl_appendages_euc{kk}] = pca(single(squeeze(agg_features_wl_euc)));
            COEFFS_feat_wl_appendages_euc{kk} = single(COEFFS_feat_wl_appendages_euc{kk});
            explained_wl_appendages_euc{kk} = single(explained_wl_appendages_euc{kk});
            coeffstruct.(coeffname_appendage) = COEFFS_feat_wl_appendages_euc{kk};
            coeffstruct.(explainedname_appendage) = explained_wl_appendages_euc{kk};
        else
            COEFFS_feat_wl_appendages_euc{kk} = single(coeffstruct.(coeffname_appendage));
            explained_wl_appendages_euc{kk} = single(coeffstruct.(explainedname_appendage));
        end
        
        agg_mean_euc = single(mean(agg_features_wl_euc, 1));
        agg_features_wl_euc = single(bsxfun(@minus, agg_features_wl_euc, agg_mean_euc));
        dyadic_spectrograms_score_wl_appendages_euc{kk} = single(agg_features_wl_euc * COEFFS_feat_wl_appendages_euc{kk});
        
        clear agg_features_wl_euc agg_mean_euc
        
        %% Replicate and adjust size for euclidean
        num_pcs_to_keep_euc = min(num_spectrogram_pcs, size(dyadic_spectrograms_score_wl_appendages_euc{kk},2));
        dyadic_spectrograms_score_wl_appendages_euc{kk} = single(repelem( ...
            dyadic_spectrograms_score_wl_appendages_euc{kk}(:,1:num_pcs_to_keep_euc), replication_factor_wl, 1));
        
        current_size_euc = size(dyadic_spectrograms_score_wl_appendages_euc{kk}, 1);
        if current_size_euc < target_size
            padding_euc = single(zeros(target_size - current_size_euc, size(dyadic_spectrograms_score_wl_appendages_euc{kk},2)));
            dyadic_spectrograms_score_wl_appendages_euc{kk} = single(cat(1, dyadic_spectrograms_score_wl_appendages_euc{kk}, padding_euc));
            clear padding_euc
        else
            dyadic_spectrograms_score_wl_appendages_euc{kk}((end-(current_size_euc-target_size)):end,:) = [];
        end
        
        final_row_euc = single(zeros(1, size(dyadic_spectrograms_score_wl_appendages_euc{kk},2)));
        dyadic_spectrograms_score_wl_appendages_euc{kk} = single(cat(1, dyadic_spectrograms_score_wl_appendages_euc{kk}, final_row_euc));
        clear final_row_euc
        
        % Clear euclidean data from memory
        clear appendage_pca_score_euc_hipassclip_kk
        
        %% Store results back to matfile - SINGLE PRECISION
        ml_matfile.ML_features.COEFFS_feat_wl_appendages_euc(1,kk) = {COEFFS_feat_wl_appendages_euc{kk}};
        ml_matfile.ML_features.dyadic_spectrograms_score_wl_appendages_euc(1,kk) = {dyadic_spectrograms_score_wl_appendages_euc{kk}};
        ml_matfile.ML_features.explained_wl_appendages_euc(1,kk) = {explained_wl_appendages_euc{kk}};
        
        ml_matfile.ML_features.COEFFS_feat_wl_appendages(1,kk) = {COEFFS_feat_wl_appendages{kk}};
        ml_matfile.ML_features.dyadic_spectrograms_score_wl_appendages(1,kk) = {dyadic_spectrograms_score_wl_appendages{kk}};
        ml_matfile.ML_features.explained_wl_appendages(1,kk) = {explained_wl_appendages{kk}};
        
    else
        %% Empty case - store empty arrays to matfile
        ml_matfile.ML_features.COEFFS_feat_wl_appendages_euc(1,kk) = {[]};
        ml_matfile.ML_features.dyadic_spectrograms_score_wl_appendages_euc(1,kk) = {[]};
        ml_matfile.ML_features.explained_wl_appendages_euc(1,kk) = {[]};
        
        ml_matfile.ML_features.COEFFS_feat_wl_appendages(1,kk) = {[]};
        ml_matfile.ML_features.dyadic_spectrograms_score_wl_appendages(1,kk) = {[]};
        ml_matfile.ML_features.explained_wl_appendages(1,kk) = {[]};
    end
end

% Clear all temporary cell arrays
clear COEFFS_feat_wl_appendages dyadic_spectrograms_score_wl_appendages explained_wl_appendages
clear COEFFS_feat_wl_appendages_euc dyadic_spectrograms_score_wl_appendages_euc explained_wl_appendages_euc

%% Save coefficients
fprintf('saving appendage coefficients !!! WAVELET !!!\n')
try
    save(coeffstruct_in,'-struct','coeffstruct','-v7.3')
catch ME
    save(coeffstruct_in,'-struct','coeffstruct','-v7.3')
end

%% Load final ML_features from matfile and cleanup
fprintf('Loading final ML_features from matfile and cleaning up...\n');
ML_features = ml_matfile.ML_features;

% % Clean up temporary file
% clear ml_matfile;
% if exist(temp_ml_file, 'file')
%     delete(temp_ml_file);
%     fprintf('Temporary matfile deleted.\n');
% end

end