function ML_features = compute_appendage_pc_demos(mocapstruct,ML_features,coeffstruct_in,overwrite_coeff)

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
ML_features.appendage_names = appendage_names;

appendage_anglegps{8} = fieldnames(ML_features.jointangle_struct);
appendage_segvals{8} = [1:numel(ML_features.all_seglengths)];

ML_features.appendage_anglegps = appendage_anglegps;
ML_features.appendage_segvals = appendage_segvals;

%% angles to include for mai tsne
appendage_gps = 8;
for kk = appendage_gps
    frames_appendage_gps{kk} = mocapstruct.modular_cluster_properties.clipped_index{8};
end
ML_features.frames_appendage_gps = frames_appendage_gps;
ML_features.appendage_gps = appendage_gps;

%% hipass and clip all of the joint angles
params.fps = 300;
params.difforder = 10;
params.medfiltorder = 3;
params.gaussorder = 2.5;

%% Pre-allocate with single precision
appendage_anglevals = cell(1,numel(appendage_names));
appendage_explained = cell(1,numel(appendage_names));
appendage_dyadic_spectrograms = cell(1,numel(appendage_names));
COEFFS_appendages = cell(1,numel(appendage_names));

appendage_lengths = cell(1,numel(appendage_names));
appendage_explained_lengths = cell(1,numel(appendage_names));
appendage_dyadic_spectrograms_lengths = cell(1,numel(appendage_names));
COEFFS_appendages_lengths = cell(1,numel(appendage_names));
COEFFS_appendages_euc = cell(1,numel(appendage_names));

save_coeffs = 0;

fprintf('starting PCA over appendages \n')
for kk = appendage_gps
    fprintf('group %f \n',kk);
    
    fieldnames_here = appendage_anglegps{kk};
    
    %% FORCE SINGLE PRECISION throughout
    appendage_anglevals{kk} = single(zeros(numel(ML_features.jointangle_struct.(fieldnames_here{1})(frames_appendage_gps{kk})), ...
        numel(fieldnames_here)));
    
    % Get the angles - ensure single precision
    for zz = 1:numel(fieldnames_here)
        appendage_anglevals{kk}(:,zz) = single(ML_features.jointangle_struct.(fieldnames_here{zz})(frames_appendage_gps{kk}));
    end
    
    meanval_ja{kk} = single(nanmean(appendage_anglevals{kk},1));
    appendage_anglevals{kk} = single(bsxfun(@minus,appendage_anglevals{kk}, meanval_ja{kk}));
    
    %% load coeffs
    coeffname_appendage = strcat('COEFFS_appendages',num2str(kk));
    explainedname_appendage = strcat('EXPLAINED_appendages',num2str(kk));
    if (~isfield(coeffstruct,coeffname_appendage) || overwrite_coeff)
        [COEFFS_appendages{kk}, ~, ~, ~,appendage_explained{kk}] = pca(single(squeeze(appendage_anglevals{kk})));
        COEFFS_appendages{kk} = single(COEFFS_appendages{kk});
        appendage_explained{kk} = single(appendage_explained{kk});
        coeffstruct.(coeffname_appendage) = COEFFS_appendages{kk};
        coeffstruct.(explainedname_appendage) = appendage_explained{kk};
        save_coeffs = 1;
    else
        COEFFS_appendages{kk} = single(coeffstruct.(coeffname_appendage));
        appendage_explained{kk} = single(coeffstruct.(explainedname_appendage));
    end
    
    appendage_dyadic_spectrograms{kk} = single(appendage_anglevals{kk} * COEFFS_appendages{kk});
    
    %% do for segments - SINGLE PRECISION
    appendage_lengths{kk} = single(zeros(numel(ML_features.jointangle_struct.(fieldnames_here{1})(frames_appendage_gps{kk})), ...
        numel(appendage_segvals{kk})));
    
    for ll = 1:numel(appendage_segvals{kk})
        appendage_lengths{kk}(:,ll) = single(ML_features.all_seglengths{appendage_segvals{kk}(ll)}(frames_appendage_gps{kk}));
    end
    
    meanval_lengths = single(nanmean(appendage_lengths{kk},1));
    appendage_lengths{kk} = single(bsxfun(@minus,appendage_lengths{kk}, meanval_lengths));
    
    %% load coeffs for lengths
    coeffname_appendage_lengths = strcat('COEFFS_appendages_lengths',num2str(kk));
    explainedname_appendage_lengths = strcat('EXPLAINED_appendages_lengths',num2str(kk));
    
    if (~isfield(coeffstruct,explainedname_appendage_lengths) || overwrite_coeff)
        [COEFFS_appendages_lengths{kk}, ~, ~, ~,appendage_explained_lengths{kk}] = pca(single(squeeze(appendage_lengths{kk})));
        COEFFS_appendages_lengths{kk} = single(COEFFS_appendages_lengths{kk});
        appendage_explained_lengths{kk} = single(appendage_explained_lengths{kk});
        coeffstruct.(coeffname_appendage_lengths) = COEFFS_appendages_lengths{kk};
        coeffstruct.(explainedname_appendage_lengths) = appendage_explained_lengths{kk};
    else
        COEFFS_appendages_lengths{kk} = single(coeffstruct.(coeffname_appendage_lengths));
        appendage_explained_lengths{kk} = single(coeffstruct.(explainedname_appendage_lengths));
    end
    appendage_dyadic_spectrograms_lengths{kk} = single(appendage_lengths{kk} * COEFFS_appendages_lengths{kk});
    
    %% get the PCS of the euclidean distances - SINGLE PRECISION
    appendage_euc_vecs{kk} = single([]);
    for ll = 1:numel(appendage_segvals{kk})
        temp_segment = single(ML_features.all_segments{appendage_segvals{kk}(ll)}(frames_appendage_gps{kk}',:));
        appendage_euc_vecs{kk} = single(cat(2, appendage_euc_vecs{kk}, temp_segment));
        clear temp_segment  % Clear immediately
    end
    
    meanval_euc{kk} = single(nanmean(appendage_euc_vecs{kk},1));
    appendage_euc_vecs{kk} = single(bsxfun(@minus,appendage_euc_vecs{kk}, meanval_euc{kk}));
    
    %% load coeffs for euclidean
    coeffname_appendage_euc = strcat('COEFFS_appendages_euc',num2str(kk));
    explainedname_appendage_euc = strcat('EXPLAINED_appendages_euc',num2str(kk));
    
    if (~isfield(coeffstruct,coeffname_appendage_euc) || overwrite_coeff)
        [COEFFS_appendages_euc{kk}, ~, ~, ~,appendage_explained_euc{kk}] = pca(single(squeeze(appendage_euc_vecs{kk})));
        COEFFS_appendages_euc{kk} = single(COEFFS_appendages_euc{kk});
        appendage_explained_euc{kk} = single(appendage_explained_euc{kk});
        coeffstruct.(coeffname_appendage_euc) = COEFFS_appendages_euc{kk};
        coeffstruct.(explainedname_appendage_euc) = appendage_explained_euc{kk};
    else
        COEFFS_appendages_euc{kk} = single(coeffstruct.(coeffname_appendage_euc));
        appendage_explained_euc{kk} = single(coeffstruct.(explainedname_appendage_euc));
    end
    appendage_dyadic_spectrograms_euc{kk} = single(appendage_euc_vecs{kk} * COEFFS_appendages_euc{kk});
    
    % Store results in ML_features as single
    ML_features.appendage_coeffs{kk} = COEFFS_appendages{kk};
    ML_features.appendage_pca_score{kk} = appendage_dyadic_spectrograms{kk};
    ML_features.appendage_pca_explained{kk} = appendage_explained{kk};
    ML_features.appendage_coeffs_euc{kk} = COEFFS_appendages_euc{kk};
    ML_features.appendage_pca_score_euc{kk} = appendage_dyadic_spectrograms_euc{kk};
    ML_features.appendage_pca_explained_euc{kk} = appendage_explained_euc{kk};
    ML_features.appendage_coeffs_lengths{kk} = COEFFS_appendages_lengths{kk};
    ML_features.appendage_pca_score_lengths{kk} = appendage_dyadic_spectrograms_lengths{kk};
    ML_features.appendage_pca_explained_lengths{kk} = appendage_explained_lengths{kk};
    
    % CRITICAL: Clear large temporary variables immediately
    clear appendage_euc_vecs meanval_lengths
end

% Clear all temporary cell arrays
clear appendage_anglevals appendage_lengths appendage_dyadic_spectrograms
clear appendage_dyadic_spectrograms_lengths appendage_dyadic_spectrograms_euc
clear COEFFS_appendages COEFFS_appendages_lengths COEFFS_appendages_euc

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
    fieldnames_here = appendage_anglegps{kk};
    
    temp_anglevals = single(zeros(numel(ML_features.jointangle_struct.(fieldnames_here{1})), numel(appendage_anglegps{kk})));
    for zz = 1:numel(fieldnames_here)
        temp_anglevals(:,zz) = single(ML_features.jointangle_struct.(fieldnames_here{zz}));
    end
    
    temp_anglevals = single(bsxfun(@minus, temp_anglevals, meanval_ja{kk}) * ML_features.appendage_coeffs{kk});
    temp_anglevals(isnan(temp_anglevals) | isinf(temp_anglevals)) = 0;
    
    ML_features.appendage_joint_angles_pcs_hipass{kk} = hipass_clip_cell(temp_anglevals, frames_appendage_gps{kk}, params);
    clear temp_anglevals  % Clear immediately
end

%% Loop over segment vectors - SINGLE PRECISION
for kk = appendage_gps
    fprintf('hipass clipping group %f EUCLIDEAN \n',kk);
    
    temp_euc_vecs = single([]);
    for ll = 1:numel(appendage_segvals{kk})
        temp_segment = single(ML_features.all_segments{appendage_segvals{kk}(ll)}(:,:));
        temp_euc_vecs = single(cat(2, temp_euc_vecs, temp_segment));
        clear temp_segment
    end
    
    temp_euc_vecs = single((temp_euc_vecs - meanval_euc{kk}) * ML_features.appendage_coeffs_euc{kk});
    ML_features.appendage_pca_score_euc_hipassclip{kk} = hipass_clip_cell(temp_euc_vecs, frames_appendage_gps{kk}, params);
    clear temp_euc_vecs
end

%% Optional timescales processing - SINGLE PRECISION
do_timescales = 0;
if do_timescales
    timescales = [10,33,100];
    for zz = 1:numel(timescales)
        params.gaussorder = timescales(zz)./2;
        gfilter = single(fspecial('gaussian',[timescales(zz)*6 1], params.gaussorder));
        
        appendage_pca_score_smoothed = cell(1,numel(appendage_anglegps));
        appendage_pca_score_smoothed_lengths = cell(1,numel(appendage_anglegps));
        
        for kk = 1:numel(appendage_anglegps)
            appendage_pca_score_smoothed{kk} = single(convn(ML_features.appendage_pca_score{kk}, gfilter, 'same'));
            appendage_pca_score_smoothed_lengths{kk} = single(convn(ML_features.appendage_pca_score_lengths{kk}, gfilter, 'same'));
        end
        
        ML_features.(strcat('appendage_pca_score',num2str(timescales(zz)))) = appendage_pca_score_smoothed;
        ML_features.(strcat('appendage_pca_score_lengths',num2str(timescales(zz)))) = appendage_pca_score_smoothed_lengths;
        
        clear appendage_pca_score_smoothed appendage_pca_score_smoothed_lengths gfilter
    end
end

% Final cleanup
clear meanval_ja meanval_euc appendage_explained appendage_explained_lengths

end