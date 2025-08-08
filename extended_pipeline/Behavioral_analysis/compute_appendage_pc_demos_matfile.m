function ML_features = compute_appendage_pc_demos_matfile(mocapstruct,  coeffstruct_in, overwrite_coeff, MLmatobjfile)
% Memory-optimized version of compute_appendage_pc_demos using matfile objects
% Only loads specific fields when needed to minimize memory usage

%% load in the coefficients
if exist(coeffstruct_in,'file') && overwrite_coeff == 0
    try
        coeffstruct = load(coeffstruct_in);
    catch ME
        coeffstruct = load(coeffstruct_in);
    end
else
    coeffstruct = struct();
end

%% specify the specific angles for the different appendages
appendage_names = {'Head','axial','LArm','Rarm','LLeg','RLeg','','globall','trunk'};
%ml_matfile = load(MLmatobjfile);

% Load only the jointangle_struct field we need
fprintf('Loading only jointangle_struct from matfile...\n');
load(MLmatobjfile,"jointangle_struct");

% Load only the segment data we need
fprintf('Loading only all_seglengths and all_segments from matfile...\n');
load(MLmatobjfile, "all_seglengths");
load(MLmatobjfile, "all_segments");


appendage_anglegps{8} = fieldnames(jointangle_struct);
appendage_segvals{8} = [1:numel(all_seglengths)];

%% angles to include for mai tsne
appendage_gps = 8;
for kk = appendage_gps
    frames_appendage_gps{kk} = mocapstruct.modular_cluster_properties.clipped_index{8};
end

%% hipass and clip all of the joint angles
params.fps = 300;
params.difforder = 10;
params.medfiltorder = 3;
params.gaussorder = 2.5;

save_coeffs = 0;

fprintf('starting PCA over appendages \n')
for kk = appendage_gps
    fprintf('group %f \n',kk);
    % init again structure cleaning up memory
    ml_matfile = struct();

    fieldnames_here = appendage_anglegps{kk};

    %% FORCE SINGLE PRECISION throughout
    appendage_anglevals_kk = single(zeros(numel(jointangle_struct.(fieldnames_here{1})(frames_appendage_gps{kk})), ...
        numel(fieldnames_here)));

    % Get the angles - ensure single precision
    for zz = 1:numel(fieldnames_here)
        appendage_anglevals_kk(:,zz) = single(jointangle_struct.(fieldnames_here{zz})(frames_appendage_gps{kk}));
    end

    meanval_ja_kk = single(nanmean(appendage_anglevals_kk,1));
    appendage_anglevals_kk = single(bsxfun(@minus,appendage_anglevals_kk, meanval_ja_kk));

    %% load coeffs
    coeffname_appendage = strcat('COEFFS_appendages',num2str(kk));
    explainedname_appendage = strcat('EXPLAINED_appendages',num2str(kk));
    if (~isfield(coeffstruct,coeffname_appendage) || overwrite_coeff)
        [COEFFS_appendages_kk, ~, ~, ~,appendage_explained_kk] = pca(single(squeeze(appendage_anglevals_kk)));
        COEFFS_appendages_kk = single(COEFFS_appendages_kk);
        appendage_explained_kk = single(appendage_explained_kk);
        coeffstruct.(coeffname_appendage) = COEFFS_appendages_kk;
        coeffstruct.(explainedname_appendage) = appendage_explained_kk;
        save_coeffs = 1;
    else
        COEFFS_appendages_kk = single(coeffstruct.(coeffname_appendage));
        appendage_explained_kk = single(coeffstruct.(explainedname_appendage));
    end

    appendage_dyadic_spectrograms_kk = single(appendage_anglevals_kk * COEFFS_appendages_kk);

    %% do for segments - SINGLE PRECISION
    appendage_lengths_kk = single(zeros(numel(jointangle_struct.(fieldnames_here{1})(frames_appendage_gps{kk})), ...
        numel(appendage_segvals{kk})));

    for ll = 1:numel(appendage_segvals{kk})
        appendage_lengths_kk(:,ll) = single(all_seglengths{appendage_segvals{kk}(ll)}(frames_appendage_gps{kk}));
    end

    meanval_lengths = single(nanmean(appendage_lengths_kk,1));
    appendage_lengths_kk = single(bsxfun(@minus,appendage_lengths_kk, meanval_lengths));

    %% load coeffs for lengths
    coeffname_appendage_lengths = strcat('COEFFS_appendages_lengths',num2str(kk));
    explainedname_appendage_lengths = strcat('EXPLAINED_appendages_lengths',num2str(kk));

    if (~isfield(coeffstruct,explainedname_appendage_lengths) || overwrite_coeff)
        [COEFFS_appendages_lengths_kk, ~, ~, ~,appendage_explained_lengths_kk] = pca(single(squeeze(appendage_lengths_kk)));
        COEFFS_appendages_lengths_kk = single(COEFFS_appendages_lengths_kk);
        appendage_explained_lengths_kk = single(appendage_explained_lengths_kk);
        coeffstruct.(coeffname_appendage_lengths) = COEFFS_appendages_lengths_kk;
        coeffstruct.(explainedname_appendage_lengths) = appendage_explained_lengths_kk;
    else
        COEFFS_appendages_lengths_kk = single(coeffstruct.(coeffname_appendage_lengths));
        appendage_explained_lengths_kk = single(coeffstruct.(explainedname_appendage_lengths));
    end
    appendage_dyadic_spectrograms_lengths_kk = single(appendage_lengths_kk * COEFFS_appendages_lengths_kk);

    %% get the PCS of the euclidean distances - SINGLE PRECISION
    appendage_euc_vecs_kk = single([]);
    for ll = 1:numel(appendage_segvals{kk})
        temp_segment = single(all_segments{appendage_segvals{kk}(ll)}(frames_appendage_gps{kk}',:));
        appendage_euc_vecs_kk = single(cat(2, appendage_euc_vecs_kk, temp_segment));
        clear temp_segment  % Clear immediately
    end

    meanval_euc_kk = single(nanmean(appendage_euc_vecs_kk,1));
    appendage_euc_vecs_kk = single(bsxfun(@minus,appendage_euc_vecs_kk, meanval_euc_kk));

    %% load coeffs for euclidean
    coeffname_appendage_euc = strcat('COEFFS_appendages_euc',num2str(kk));
    explainedname_appendage_euc = strcat('EXPLAINED_appendages_euc',num2str(kk));

    if (~isfield(coeffstruct,coeffname_appendage_euc) || overwrite_coeff)
        [COEFFS_appendages_euc_kk, ~, ~, ~,appendage_explained_euc_kk] = pca(single(squeeze(appendage_euc_vecs_kk)));
        COEFFS_appendages_euc_kk = single(COEFFS_appendages_euc_kk);
        appendage_explained_euc_kk = single(appendage_explained_euc_kk);
        coeffstruct.(coeffname_appendage_euc) = COEFFS_appendages_euc_kk;
        coeffstruct.(explainedname_appendage_euc) = appendage_explained_euc_kk;
    else
        COEFFS_appendages_euc_kk = single(coeffstruct.(coeffname_appendage_euc));
        appendage_explained_euc_kk = single(coeffstruct.(explainedname_appendage_euc));
    end
    appendage_dyadic_spectrograms_euc_kk = single(appendage_euc_vecs_kk * COEFFS_appendages_euc_kk);

    % Store results directly to individual matfile fields - TRUE MEMORY OPTIMIZATION
    fprintf('Saving results to individual matfile fields...\n');

    % Store results in cell arrays indexed by kk

    ml_matfile.appendage_pca_score{kk} = appendage_dyadic_spectrograms_kk;
    ml_matfile.appendage_pca_explained{kk} = appendage_explained_kk;
    ml_matfile.appendage_coeffs_euc{kk} = COEFFS_appendages_euc_kk;
    ml_matfile.appendage_pca_score_euc{kk} = appendage_dyadic_spectrograms_euc_kk;
    ml_matfile.appendage_pca_explained_euc{kk} = appendage_explained_euc_kk;
    ml_matfile.appendage_coeffs_lengths{kk} = COEFFS_appendages_lengths_kk;
    ml_matfile.appendage_pca_score_lengths{kk} = appendage_dyadic_spectrograms_lengths_kk;
    ml_matfile.appendage_pca_explained_lengths{kk} = appendage_explained_lengths_kk;
    % Store mean values for later use
    save(MLmatobjfile, '-struct','ml_matfile', '-append')
    % CRITICAL: Clear large temporary variables immediately
    ml_matfile = rmfield(ml_matfile,fieldnames(ml_matfile)); %save up memory
    ml_matfile.appendage_coeffs{kk} = COEFFS_appendages_kk;
    ml_matfile.meanval_ja{kk} = meanval_ja_kk;
    ml_matfile.meanval_euc{kk} = meanval_euc_kk;

     clear appendage_euc_vecs_kk appendage_anglevals_kk appendage_lengths_kk
    clear appendage_dyadic_spectrograms_kk appendage_dyadic_spectrograms_lengths_kk appendage_dyadic_spectrograms_euc_kk
    clear COEFFS_appendages_kk COEFFS_appendages_lengths_kk COEFFS_appendages_euc_kk
    clear meanval_lengths
end


% Save metadata to matfile
fprintf('Saving metadata to matfile...\n');
ml_matfile.appendage_names = appendage_names;
ml_matfile.appendage_anglegps = appendage_anglegps;
ml_matfile.appendage_segvals = appendage_segvals;
ml_matfile.frames_appendage_gps = frames_appendage_gps;
ml_matfile.appendage_gps = appendage_gps;
% Append the metadata to the matfile
save(MLmatobjfile, '-struct','ml_matfile', '-append')

%% save coeffs
if save_coeffs
    fprintf('saving appendage coefficients \n')
    try
        save(coeffstruct_in,'-struct','coeffstruct','-v7.3')
    catch ME
        save(coeffstruct_in,'-struct','coeffstruct','-v7.3')
    end
end

%% Apply PCA to smoothed dynamics - SINGLE PRECISION
for kk = appendage_gps
    fprintf('hipass clipping group %f JOINT ANGLE PCs \n',kk);

    % Load only the specific data we need from matfile
    fieldnames_here = appendage_anglegps{kk};
    meanval_ja_kk = ml_matfile.meanval_ja{kk};
    appendage_coeffs_kk = ml_matfile.appendage_coeffs{kk};

    temp_anglevals = single(zeros(numel(jointangle_struct.(fieldnames_here{1})), numel(appendage_anglegps{kk})));
    for zz = 1:numel(fieldnames_here)
        temp_anglevals(:,zz) = single(jointangle_struct.(fieldnames_here{zz}));
    end

    temp_anglevals = single(bsxfun(@minus, temp_anglevals, meanval_ja_kk) * appendage_coeffs_kk);
    temp_anglevals(isnan(temp_anglevals) | isinf(temp_anglevals)) = 0;

    appendage_joint_angles_pcs_hipass_kk = hipass_clip_cell(temp_anglevals, frames_appendage_gps{kk}, params);

    % Save result in cell array
    ml_matfile.appendage_joint_angles_pcs_hipass{kk} = appendage_joint_angles_pcs_hipass_kk;

    clear temp_anglevals appendage_joint_angles_pcs_hipass_kk meanval_ja_kk appendage_coeffs_kk % Clear immediately
end

% load coeffs
load(MLmatobjfile, "appendage_coeffs_euc");
ml_matfile.appendage_coeffs_euc = appendage_coeffs_euc;

%% Loop over segment vectors - SINGLE PRECISION
for kk = appendage_gps
    fprintf('hipass clipping group %f EUCLIDEAN \n',kk);

    % Load only the specific data we need from matfile
    % field_base = sprintf('appendage_%d_', kk);
    meanval_euc_kk = ml_matfile.meanval_euc{kk};
    appendage_coeffs_euc_kk = ml_matfile.appendage_coeffs_euc{kk};%(strcat(field_base, 'coeffs_euc'));

    temp_euc_vecs = single([]);
    for ll = 1:numel(appendage_segvals{kk})
        temp_segment = single(all_segments{appendage_segvals{kk}(ll)}(:,:));
        temp_euc_vecs = single(cat(2, temp_euc_vecs, temp_segment));
        clear temp_segment
    end

    temp_euc_vecs = single((temp_euc_vecs - meanval_euc_kk) * appendage_coeffs_euc_kk);
    appendage_pca_score_euc_hipassclip_kk = hipass_clip_cell(temp_euc_vecs, frames_appendage_gps{kk}, params);

    % Save result directly to matfile
    ml_matfile.pca_score_euc_hipassclip{kk} = appendage_pca_score_euc_hipassclip_kk;

    clear temp_euc_vecs appendage_pca_score_euc_hipassclip_kk meanval_euc_kk appendage_coeffs_euc_kk
end
save(MLmatobjfile,'-struct','ml_matfile','-append');
% Clear large data arrays we no longer need
clear jointangle_struct all_seglengths all_segments

% Return empty ML_features since everything is stored in matfile
ML_features = struct();
ML_features.status = 'completed_matfile_processing';

fprintf('Appendage PC computation completed with matfile optimization.\n');

end
