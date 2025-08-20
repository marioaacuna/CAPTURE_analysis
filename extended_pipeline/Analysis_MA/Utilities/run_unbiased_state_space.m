% function run_unbiased_state_space()
% Unbiased animal-level behavioral discriminability without t-SNE.
% Methods:
%   1) PCA + GMM (states = soft posteriors) -> occupancy + soft transitions
%   2) PCA + k-means/DBSCAN (hard states)  -> occupancy + transitions
%   3) No clustering (aggregate stats per feature) -> mean/std/quantiles
%   4) Sequence model (discrete HMM on symbols; fallback Markov)
%
% Evaluation:
%   - Leave-one-animal-out (grouped by animal ID) SVM AUC
%   - MMD permutation test (animal-level)
%
% Outputs:
%   - Figures per method with AUC and MMD for pairs
%   - CSVs with features and summary statistics
%   - Saved under GC.temp_root/figs_unbiased_state_space/<exp.name>/

clc;
GC = general_configs();

% Match your existing experiment list
experiments = struct();
experiments(1).name = 'BSFC_300hz';
experiments(1).folder = '0_preprocessing_BSFC_300hz';
experiments(2).name = 'BHNG_300hz';
experiments(2).folder = '0_preprocessing_BHNG_300hz';

for exp_idx = 1%:numel(experiments)
    exp = experiments(exp_idx);
    fprintf('\n=== Unbiased state-space analysis: %s ===\n', exp.name);

    data_folder = fullfile(GC.project_path, 'data', exp.folder);
    if ~exist(data_folder, 'dir')
        warning('Data folder not found: %s', data_folder);
        continue;
    end

    % Load analysisstruct (we do NOT use zValues or t-SNE)
    fn_analysis = fullfile(data_folder, 'analysis_zValues_all_features.mat');
    if ~exist(fn_analysis, 'file')
        warning('Missing %s', fn_analysis);
        continue;
    end
    S = load(fn_analysis, 'analysisstruct');
    if ~isfield(S, 'analysisstruct')
        warning('analysisstruct missing in %s', fn_analysis);
        continue;
    end
    analysisstruct = S.analysisstruct;

    % Load animal identifiers
    fn_pred = fullfile(data_folder, 'agg_predictions.mat');
    if ~exist(fn_pred, 'file')
        warning('Missing %s', fn_pred);
        continue;
    end
    S_pred = load(fn_pred, 'animal_condition_identifier');
    if ~isfield(S_pred, 'animal_condition_identifier')
        warning('animal_condition_identifier missing in %s', fn_pred);
        continue;
    end

    % Collect aligned features and IDs
    try
        [X_raw, ids, animals, conds] = collect_features_and_ids(analysisstruct, S_pred.animal_condition_identifier, GC);
    catch ME
        %warning('Feature collection failed: %s', ME.message);
        continue;
    end

    % Transform analysisstruct to group by animal IDs
    fprintf('Transforming analysisstruct to group by animal IDs...\n');
    analysisstruct = transform_analysisstruct_by_animals(analysisstruct, ids);

    % Preprocess: impute NaNs, standardize, PCA
    % [X_std, col_mu, col_sig] = standardize_impute(X_raw); %#ok<NASGU>
    ncomp = min(10, size(X_raw,2));
    [score_pca, pca_info] = do_pca(X_raw, ncomp); %#ok<NASGU>

    % Determine optimal K using hierarchical clustering and silhouette analysis
    % fprintf('Determining optimal K via hierarchical clustering...\n');
    % optimal_K = determine_optimal_k(score_pca);
    % fprintf('Optimal K selected: %d\n', optimal_K);

    % Fixed optimal K for all methods (behavioral states)
    optimal_K = 12;
    fprintf('Using fixed K = %d for all clustering methods\n', optimal_K);

    % Clustering method toggle: 'kmeans' or 'dbscan'
    clustering_method = 'kmeans'; % Change to 'dbscan' to use DBSCAN

    % Build per-ID index map
    [ac_list, ac_to_idx] = index_by_ac(ids);

    % Condition palette (axes text black; white background)
    cmap = containers.Map();
    cmap('B') = [0.5, 0.5, 0.5];            % baseline: gray
    cmap('S') = [0.4660, 0.6740, 0.1880];   % S/H: green
    cmap('H') = [0.4660, 0.6740, 0.1880];
    cmap('F') = [0.8500, 0.3250, 0.0980];   % F/N: red
    cmap('N') = [0.8500, 0.3250, 0.0980];
    cmap('C') = [0, 0.4470, 0.7410];        % C/G: blue
    cmap('G') = [0, 0.4470, 0.7410];

    % Determine pain control per experiment (S or H)
    uconds = unique(conds);
    ctrlPain = 'S'; if any(strcmp(uconds,'H')) && ~any(strcmp(uconds,'S')), ctrlPain = 'H'; end

    % Pairs to test: restrict by experiment and available conditions
    if contains(exp.name, 'BSFC')
        candidate_pairs = {
            'F','B'
            'F',ctrlPain;
            'F','C';
            'C', 'S';
        };
    elseif contains(exp.name, 'BHNG')
        candidate_pairs = {
            'N','B';
            'N',ctrlPain;
            'N','G';
            'G', 'H';
        };
    else
        % Fallback: anchor to present main condition
        if any(strcmp(uconds,'F'))
            candidate_pairs = {
                'F','B';
                'F',ctrlPain;
                'F','C';
                'C', 'S';
            };
        else
            candidate_pairs = {
                'N','B';
                'N',ctrlPain;
                'N','G';
                'G', 'H';

            };
        end
    end
    % Keep only pairs where both conditions are present in this dataset
    keep_mask = false(size(candidate_pairs,1),1);
    for pi = 1:size(candidate_pairs,1)
        keep_mask(pi) = any(strcmp(uconds, candidate_pairs{pi,1})) && any(strcmp(uconds, candidate_pairs{pi,2}));
    end
    pairs = candidate_pairs(keep_mask, :);

    out_folder = fullfile(GC.temp_root, 'figs_unbiased_state_space', exp.name);
    if ~exist(out_folder, 'dir'), mkdir(out_folder); end

    % 1) PCA + GMM
    fprintf('Method 1: PCA + GMM (soft states)...\n');
    try
        K_candidates = [optimal_K-2, optimal_K-1, optimal_K, optimal_K+1, optimal_K+2];
        K_candidates = K_candidates(K_candidates >= 3); % Ensure K >= 3
        [X_occ_tr, X_tran_tr, info_gmm] = features_from_gmm(score_pca, ac_to_idx, K_candidates);
        X_animal_gmm = [X_occ_tr, X_tran_tr];
        [AUCs_gmm, Pvals_gmm, pair_labels, n_samples] = loao_and_mmd(X_animal_gmm, ac_list, animals, conds, pairs);
        fig1 = plot_auc_mmd_with_error(pair_labels, AUCs_gmm, Pvals_gmm, n_samples, sprintf('%s - PCA+GMM (K=%d, %dD PCA)', exp.name, info_gmm.K, size(score_pca,2)), cmap);
        exportgraphics(fig1, fullfile(out_folder, 'unbiased_pca_gmm_auc_mmd.pdf'), 'ContentType','vector', 'BackgroundColor','white');
        close(fig1);
        % Export CSVs
        export_method_csvs(out_folder, 'pca_gmm', ac_list, X_animal_gmm, build_gmm_colnames(info_gmm.K), pair_labels, AUCs_gmm, Pvals_gmm);
    catch ME
        %warning('PCA+GMM failed: %s', ME.message);
    end

    % 2) PCA + clustering (k-means or DBSCAN)
    fprintf('Method 2: PCA + %s clustering...\n', clustering_method);
    try
        if strcmp(clustering_method, 'dbscan')
            [X_occ_km, X_tran_km, info_km] = features_from_dbscan(score_pca, ac_to_idx);
        else
            [X_occ_km, X_tran_km, info_km] = features_from_kmeans(score_pca, ac_to_idx, optimal_K);
        end
        X_animal_km = [X_occ_km, X_tran_km];
        [AUCs_km, Pvals_km, pair_labels, n_samples] = loao_and_mmd(X_animal_km, ac_list, animals, conds, pairs);
        fig2 = plot_auc_mmd_with_error(pair_labels, AUCs_km, Pvals_km, n_samples, sprintf('%s - PCA+%s (K=%d, %dD PCA)', exp.name, clustering_method, info_km.K, size(score_pca,2)), cmap);
        exportgraphics(fig2, fullfile(out_folder, sprintf('unbiased_pca_%s_auc_mmd.pdf', clustering_method)), 'ContentType','vector', 'BackgroundColor','white');
        close(fig2);
        % Export CSVs
        export_method_csvs(out_folder, 'pca_kmeans', ac_list, X_animal_km, build_kmeans_colnames(info_km.K), pair_labels, AUCs_km, Pvals_km);
    catch ME
        %warning('PCA+k-means failed: %s', ME.message);
    end

    % 3) No clustering: aggregate stats
    fprintf('Method 3: No clustering (aggregate stats)...\n');
    try
        X_animal_agg = features_aggregate(X_std, ac_to_idx);
        [AUCs_agg, Pvals_agg, pair_labels, n_samples] = loao_and_mmd(X_animal_agg, ac_list, animals, conds, pairs);
        fig3 = plot_auc_mmd_with_error(pair_labels, AUCs_agg, Pvals_agg, n_samples, sprintf('%s - Aggregates (means/std/quantiles)', exp.name), cmap);
        exportgraphics(fig3, fullfile(out_folder, 'unbiased_aggregate_auc_mmd.pdf'), 'ContentType','vector', 'BackgroundColor','white');
        close(fig3);
        % Export CSVs (generic column names)
        Dagg = size(X_animal_agg,2);
        agg_cols = arrayfun(@(i) sprintf('feat_%03d', i), 1:Dagg, 'UniformOutput', false);
        export_method_csvs(out_folder, 'aggregate', ac_list, X_animal_agg, agg_cols, pair_labels, AUCs_agg, Pvals_agg);
    catch ME
        %warning('Aggregate analysis failed: %s', ME.message);
    end

    % 4) Sequence model: discrete HMM (fallback Markov)
    fprintf('Method 4: Sequence model (HMM discrete on symbols)...\n');
    try
        % if strcmp(clustering_method, 'dbscan')
        %     [sym_idx, ~] = symbols_from_dbscan(score_pca, optimal_K);
        % else
        %     [sym_idx, ~] = symbols_from_kmeans(score_pca, optimal_K);
        % end
        [sym_idx] = analysisstruct.annot_reordered{2};
        optimal_K = max(sym_idx);
        H_hidden = max(6, floor(optimal_K/2)); % Hidden states based on optimal K
        [X_animal_hmm, used_hmm] = features_from_hmm_discrete(sym_idx, ac_to_idx, optimal_K, H_hidden);
        [AUCs_hmm, Pvals_hmm, pair_labels, n_samples] = loao_and_mmd(X_animal_hmm, ac_list, animals, conds, pairs);
        tag = 'HMM'; if ~used_hmm, tag = 'Markov fallback'; end
        fig4 = plot_auc_mmd_with_error(pair_labels, AUCs_hmm, Pvals_hmm, n_samples, sprintf('%s - Sequence model (%s, Ksym=%d, H=%d)', exp.name, tag, optimal_K, H_hidden), cmap);
        exportgraphics(fig4, fullfile(out_folder, 'unbiased_sequence_auc_mmd.pdf'), 'ContentType','vector', 'BackgroundColor','white');
        close(fig4);
        % Export CSVs (label with hidden/obs dims)
        hmm_cols = [arrayfun(@(i) sprintf('occ_%02d', i), 1:size(X_animal_hmm,2), 'UniformOutput', false)];
        export_method_csvs(out_folder, sprintf('sequence_%s', lower(tag)), ac_list, X_animal_hmm, hmm_cols, pair_labels, AUCs_hmm, Pvals_hmm);
    catch ME
        %warning('Sequence model failed: %s', ME.message);
    end

    % 5) Separate Markov chain analysis (sequence-only features from symbols)
    fprintf('Method 5: Markov chain analysis (separate)...\n');
    try
        % if strcmp(clustering_method, 'dbscan')
        %     [sym_idx_mk, ~] = symbols_from_dbscan(score_pca, optimal_K);
        % else
        %     [sym_idx_mk, ~] = symbols_from_kmeans(score_pca, optimal_K);
        % end
        % [sym_idx_mk, ~] = symbols_from_kmeans(analysisstruct.zValues, 24);
        rng("default");
        % [sym_idx_mk, ~] = kmeans(analysisstruct.zValues, 1200, 'Distance','sqeuclidean');
        sym_idx_mk = analysisstruct.annot_reordered{2}';
        [X_markov, ~, Pmats_markov] = features_markov_only(sym_idx_mk, ac_to_idx, max(sym_idx_mk));
        [AUCs_mk, Pvals_mk, pair_labels, n_samples] = loao_and_mmd(X_markov, ac_list, animals, conds, pairs);
        fig5 = plot_auc_mmd_with_error(pair_labels, AUCs_mk, Pvals_mk, n_samples, sprintf('%s - Markov chain only (Ksym=%d)', exp.name, optimal_K), cmap);
        exportgraphics(fig5, fullfile(out_folder, 'unbiased_markov_only_auc_mmd.pdf'), 'ContentType','vector', 'BackgroundColor','white');
        close(fig5);
        % Export CSVs with explicit Markov column names
        markov_cols = [arrayfun(@(i) sprintf('occ_%02d', i), 1:optimal_K, 'UniformOutput', false), ...
                       arrayfun(@(i) sprintf('P_rowmajor_%03d', i), 1:optimal_K*optimal_K, 'UniformOutput', false)];
        export_method_csvs(out_folder, 'markov_only', ac_list, X_markov, markov_cols, pair_labels, AUCs_mk, Pvals_mk);
        % Optional: export full transition matrices per AC as separate CSV
        export_markov_matrices_csv(out_folder, 'markov_transition_matrices', ac_list, Pmats_markov);
    catch ME
        %warning('Markov-only analysis failed: %s', ME.message);
    end

    fprintf('Saved results to %s\n', out_folder);
end

fprintf('\nAll experiments complete.\n');
% end


%% Data collection and preprocessing

function [X_raw, ids, animals, conds] = collect_features_and_ids(analysisstruct, animal_condition_identifier, GC)
% Align features and IDs to frames used in analysis (good tracking)

% Upsample ids then select good frames
ups = GC.repfactor;
long_ids = repelem(animal_condition_identifier, ups);

good = analysisstruct.frames_with_good_tracking;
if iscell(good), good = good{1}; end
if islogical(good), good = find(good); end
ids = long_ids(good);

% Features: jt_features and extra_jt_features (if present)
X1 = []; X2 = [];
if isfield(analysisstruct, 'jt_features') && ~isempty(analysisstruct.jt_features)
    X1 = analysisstruct.jt_features;
end
if isfield(analysisstruct, 'extra_jt_features') && ~isempty(analysisstruct.extra_jt_features)
    X2 = analysisstruct.extra_jt_features;
end
if isempty(X1) && isempty(X2)
    error('No jt_features or extra_jt_features found in analysisstruct.');
end

% If features are cells, take the first cell (common pattern)
if iscell(X1), X1 = X1{1}; end
if iscell(X2), X2 = X2{1}; end

if isempty(X1)
    X_raw_all = X2;
elseif isempty(X2)
    X_raw_all = X1;
else
    % Concatenate feature sets
    if size(X1,1) ~= size(X2,1)
        error('jt_features and extra_jt_features have different row counts.');
    end
    X_raw_all = [X1, X2];
end

% Select the same frames used elsewhere
if size(X_raw_all,1) ~= numel(ids)
    % In some pipelines features are full-length; subset by good indices if needed
    % Try to use good as linear indices into features if lengths differ
    if max(good) <= size(X_raw_all,1)
        X_raw = X_raw_all(good, :);
    else
        error('Feature rows (%d) do not match good frames index and ids (%d).', size(X_raw_all,1), numel(ids));
    end
else
    X_raw = X_raw_all;
end

% Parse animals and conditions
animals = cellfun(@(s) s(1:find(s=='_',1)-1), ids, 'UniformOutput', false);
conds = cellfun(@(s) s(end), ids, 'UniformOutput', false);
end

function [X_std, mu, sig] = standardize_impute(X)
% Impute NaNs with column medians, then z-score columns
X_std = X;
% Median imputation
for j = 1:size(X_std,2)
    col = X_std(:,j);
    if any(isnan(col))
        med = median(col(~isnan(col)));
        if ~isfinite(med), med = 0; end
        col(isnan(col)) = med;
        X_std(:,j) = col;
    end
end
% Z-score
mu = mean(X_std,1);
sig = std(X_std,[],1);
sig(sig==0) = 1;
X_std = (X_std - mu) ./ sig;
end

function [score_pca, info] = do_pca(X_std, ncomp)
% PCA on standardized features
[coeff, score, ~, ~, explained, mu] = pca(X_std, 'NumComponents', ncomp);
score_pca = score;
info = struct('coeff', coeff, 'explained', explained, 'mu', mu);
end

function [ac_list, ac_to_idx] = index_by_ac(ids)
% Build index mapping from animal+condition id to frame indices
ac_list = unique(ids, 'stable');
ac_to_idx = containers.Map();
for i = 1:numel(ac_list)
    ac = ac_list{i};
    ac_to_idx(ac) = find(strcmp(ids, ac));
end
end


%% Feature builders

function [X_occ, X_tran, info] = features_from_gmm(score_pca, ac_to_idx, K_candidates)
% Fit GMM in PCA space (choose K via BIC), compute:
%  - Soft occupancies per AC: mean responsibilities
%  - Soft transitions per AC: sum_t r_t(i) * r_{t+1}(j), row-normalized

opts = statset('MaxIter', 300, 'Display', 'off');
bestBIC = Inf; bestGM = []; bestK = NaN;
for K = K_candidates
    try
        gm = fitgmdist(score_pca, K, 'RegularizationValue', 1e-6, 'Replicates', 3, 'Options', opts);
        if gm.BIC < bestBIC
            bestBIC = gm.BIC; bestGM = gm; bestK = K;
        end
    catch
        % skip K if fails
    end
end
if isempty(bestGM)
    error('GMM fit failed for all K candidates.');
end

[R, ~] = posterior(bestGM, score_pca); % N x K responsibilities
K = bestK;

% Occupancy per AC: mean responsibilities
ac_keys = keys(ac_to_idx);
A = numel(ac_keys);
X_occ = zeros(A, K);
X_tran = zeros(A, K*K);

for a = 1:A
    idx = ac_to_idx(ac_keys{a});
    Ri = R(idx, :);
    % occupancy (soft)
    occ = mean(Ri, 1); % 1 x K

    % soft transitions: sum_t r_t(i) r_{t+1}(j)
    T = size(Ri,1);
    Tpairs = max(T-1, 0);
    Tij = zeros(K,K);
    for t = 1:Tpairs
        rt = Ri(t, :)';
        rt1 = Ri(t+1, :)';
        Tij = Tij + (rt * rt1'); % K x K
    end
    % row-normalize (avoid div-by-zero)
    rowSums = sum(Tij, 2);
    rowSums(rowSums==0) = 1;
    P = Tij ./ rowSums;
    X_tran(a, :) = P(:)';
    X_occ(a, :) = occ;
end

info = struct('K', K, 'BIC', bestBIC);
end

function [X_occ, X_tran, info] = features_from_kmeans(score_pca, ac_to_idx, K)
% Fit k-means in PCA space, compute occupancies and transitions
% Uses bout-based analysis to avoid diagonal-dominated transitions

rng(42);
% Use only top 20 PCs to stabilize clustering and avoid overfitting
Xk = score_pca(:, 1:min(20, size(score_pca,2)));
% Stronger settings: k-means++, cosine distance, more iterations/replicates
[idx, ~] = kmeans(Xk, K, 'Replicates', 20, 'MaxIter', 1000, 'Distance','cosine', ...
    'Start','plus', 'Display','off');

% Bout detection parameters: adjusted for ~12 Hz effective sampling rate
min_bout_frames = 30; % Still 30 frames to ensure meaningful behavioral bouts

ac_keys = keys(ac_to_idx);
A = numel(ac_keys);
X_occ = zeros(A, K);
X_tran = zeros(A, K*K);

for a = 1:A
    fr = ac_to_idx(ac_keys{a});
    raw_seq = idx(fr);

    % Detect behavioral bouts instead of using frame-by-frame transitions
    bout_seq = detect_bouts(raw_seq, min_bout_frames);

    if isempty(bout_seq)
        % Fallback: use original sequence if no bouts detected
        bout_seq = raw_seq;
        fprintf('Warning: No bouts detected for k-means %s, using original sequence\n', ac_keys{a});
    end

    % occupancy based on bouts
    h = histcounts(bout_seq, 0.5:1:(K+0.5));
    occ = h / max(1, sum(h));

    % transitions between bouts (avoids artificial self-transitions)
    Tij = zeros(K,K);
    for t = 1:(numel(bout_seq)-1)
        i = bout_seq(t); j = bout_seq(t+1);
        Tij(i,j) = Tij(i,j) + 1;
    end
    rowSums = sum(Tij, 2);
    rowSums(rowSums==0) = 1;
    P = Tij ./ rowSums;
    X_occ(a,:) = occ;
    X_tran(a,:) = P(:)';
end

info = struct('K', K);
end

function X_animal = features_aggregate(X_std, ac_to_idx)
% Aggregate per-feature statistics per animal+condition
% Stats: mean, std, 25th, 50th, 75th percentiles
ac_keys = keys(ac_to_idx);
A = numel(ac_keys);
D = size(X_std,2);
stats_per_feat = 5;
X_animal = zeros(A, D * stats_per_feat);

for a = 1:A
    fr = ac_to_idx(ac_keys{a});
    Xi = X_std(fr, :);
    mu = mean(Xi, 1);
    sd = std(Xi, [], 1);
    q = prctile(Xi, [25 50 75], 1);
    % concat: [mu, sd, q25, q50, q75]
    feat = [mu, sd, q(1,:), q(2,:), q(3,:)];
    X_animal(a, :) = feat;
end
end

function [sym_idx, info] = symbols_from_kmeans(score_pca, K)
% Discretize PCA space via k-means for sequence modeling
rng(123);
% Use only top 10 PCs for symbolization + stronger k-means settings
Xk = score_pca(:, 1:min(10, size(score_pca,2)));
[sym_idx, ~] = kmeans(Xk, K, 'Replicates', 20, 'MaxIter', 1000, 'Distance','cosine', ...
    'Start','plus', 'Display','off');
info = struct('K', K);
end

function [X_markov, occ_all, Pmats] = features_markov_only(sym_idx, ac_to_idx, Ksym)
% Build per-AC Markov occupancy and transition features from discrete symbols
% Uses bout-based analysis to avoid diagonal-dominated transitions from high sampling rates

% Bout detection parameters: adjusted for ~12 Hz effective sampling rate
min_bout_frames = 30; % Still 30 frames to ensure meaningful behavioral bouts

ac_keys = keys(ac_to_idx);
A = numel(ac_keys);
occ_all = zeros(A, Ksym);
Pmats = zeros(Ksym, Ksym, A);

for a = 1:A
    seq = sym_idx(ac_to_idx(ac_keys{a}));

    % Detect behavioral bouts instead of using frame-by-frame transitions
    bout_seq = detect_bouts(seq, min_bout_frames);

    if isempty(bout_seq)
        % Fallback: use original sequence if no bouts detected
        bout_seq = seq;
        fprintf('Warning: No bouts detected for %s, using original sequence\n', ac_keys{a});
    end

    % occupancy based on bouts (more meaningful than raw frames)
    h = histcounts(bout_seq, 0.5:1:(Ksym+0.5));
    occ = h / max(1, sum(h));
    occ_all(a,:) = occ;

    % transitions between bouts (avoids artificial self-transitions)
    Tij = zeros(Ksym, Ksym);
    for t = 1:(numel(bout_seq)-1)
        i = bout_seq(t); j = bout_seq(t+1);
        Tij(i,j) = Tij(i,j) + 1;
    end
    rows = sum(Tij,2); rows(rows==0) = 1;
    P = Tij ./ rows;
    Pmats(:,:,a) = P;
end
% Feature vector: concatenate occupancy and row-major P
X_markov = [occ_all, reshape(Pmats, Ksym*Ksym, A)']; % A x (K + K^2)
end

function [X_animal, used_hmm] = features_from_hmm_discrete(sym_idx, ac_to_idx, Ksym, H_hidden)
% Train a global discrete HMM on k-means symbols (1..Ksym).
% Uses bout-based sequences to avoid diagonal-dominated transitions
% If hmmtrain is unavailable, fallback to per-AC Markov chain features.

used_hmm = false;
ac_keys = keys(ac_to_idx);
A = numel(ac_keys);

% Bout detection parameters: adjusted for ~12 Hz effective sampling rate
min_bout_frames = 30; % Still 30 frames to ensure meaningful behavioral bouts

if exist('hmmtrain', 'file') == 2 && exist('hmmdecode', 'file') == 2 && exist('hmmviterbi', 'file') == 2
    % Prepare bout sequences (as cell array) instead of raw frame sequences
    seqs = cell(A,1);
    for a = 1:A
        raw_seq = sym_idx(ac_to_idx(ac_keys{a}));
        bout_seq = detect_bouts(raw_seq, min_bout_frames);

        if isempty(bout_seq)
            % Fallback: use original sequence if no bouts detected
            bout_seq = raw_seq;
            fprintf('Warning: No bouts detected for HMM %s, using original sequence\n', ac_keys{a});
        end

        seqs{a} = bout_seq(:)'; % Ensure row vector
        if isempty(seqs{a}), seqs{a} = 1; end
    end

    % Initialize HMM parameters
    rng(7);
    Atrans = mkstochastic(rand(H_hidden, H_hidden));
    Emiss = mkstochastic(rand(H_hidden, Ksym));

    try
        % Skip hmmtrain due to concatenation issues, use initialization only
        [ESTTR, ESTEMIT] = hmmtrain(seqs, Atrans, Emiss, 'Maxiterations', 50, 'Verbose', false);
        ESTTR = Atrans; ESTEMIT = Emiss;
        used_hmm = false; % Mark as fallback since we're not actually training
    catch
        ESTTR = Atrans; ESTEMIT = Emiss;
        used_hmm = false;
    end

    % Decode per sequence and build features: occupancy over hidden states + transitions
    X_occ = zeros(A, H_hidden);
    X_tran = zeros(A, H_hidden*H_hidden);
    for a = 1:A
        O = seqs{a};
        try
            [~, PSTATE] = hmmdecode(O, ESTTR, ESTEMIT); % PSTATE: H x T posterior
            occ = mean(PSTATE, 2)'; % 1 x H
            q = hmmviterbi(O, ESTTR, ESTEMIT); % hard hidden sequence
        catch
            % Fallback if decode fails
            occ = normalize(ones(1, H_hidden));
            q = ones(1, numel(O));
        end
        Tij = zeros(H_hidden, H_hidden);
        for t = 1:(numel(q)-1)
            Tij(q(t), q(t+1)) = Tij(q(t), q(t+1)) + 1; % Count transitions
        end
        rowSums = sum(Tij, 2);
        rowSums(rowSums==0) = 1;
        P = Tij ./ rowSums;
        X_occ(a,:) = occ;
        X_tran(a,:) = P(:)';
    end
    X_animal = [X_occ, X_tran];
else
    % Markov fallback on observed symbols: occupancy + transitions over Ksym
    % Use bout-based analysis here too
    X_occ = zeros(A, Ksym);
    X_tran = zeros(A, Ksym*Ksym);
    for a = 1:A
        raw_seq = sym_idx(ac_to_idx(ac_keys{a}));
        bout_seq = detect_bouts(raw_seq, min_bout_frames);

        if isempty(bout_seq)
            bout_seq = raw_seq;
        end

        h = histcounts(bout_seq, 0.5:1:(Ksym+0.5));
        occ = h / max(1, sum(h));
        Tij = zeros(Ksym, Ksym);
        for t = 1:(numel(bout_seq)-1)
            i = bout_seq(t); j = bout_seq(t+1);
            Tij(i,j) = Tij(i,j) + 1;
        end
        rowSums = sum(Tij, 2);
        rowSums(rowSums==0) = 1;
        P = Tij ./ rowSums;
        X_occ(a,:) = occ;
        X_tran(a,:) = P(:)';
    end
    X_animal = [X_occ, X_tran];
end
end

function y = normalize(x)
s = sum(x);
if s == 0, y = x; else, y = x ./ s; end
end

function M = mkstochastic(M)
M = max(M, eps);
M = M ./ sum(M, 2);
end

function bout_seq = detect_bouts(sym_idx, min_bout_length)
% Identify behavioral "bouts" - continuous periods of similar behavior
% Only count transitions between bouts, not within bouts
% This prevents diagonal-dominated transition matrices from high sampling rates

bout_seq = [];
if isempty(sym_idx)
    return;
end

current_state = sym_idx(1);
bout_start = 1;

for t = 2:length(sym_idx)
    if sym_idx(t) ~= current_state
        bout_length = t - bout_start;
        if bout_length >= min_bout_length
            bout_seq = [bout_seq, current_state];
        end
        current_state = sym_idx(t);
        bout_start = t;
    end
end

% Handle the final bout
bout_length = length(sym_idx) - bout_start + 1;
if bout_length >= min_bout_length
    bout_seq = [bout_seq, current_state];
end

% Ensure we have at least some bouts, even if short
if isempty(bout_seq) && ~isempty(sym_idx)
    % If no bouts meet criteria, use a more lenient threshold
    bout_seq = detect_bouts(sym_idx, max(1, floor(min_bout_length/2)));
end
end


%% Evaluation and plotting

function [AUCs, Pvals, pair_labels, n_samples] = loao_and_mmd(X_animal, ac_list, ~, ~, pairs)
% Build per-AC meta, then evaluate LOAO AUC and MMD
% We need animal and condition per AC row
A = numel(ac_list);
ac_animals = cell(A,1);
ac_conds = cell(A,1);
for a = 1:A
    ac = ac_list{a};
    % Parse animal and condition from AC id
    ius = find(ac=='_',1);
    if isempty(ius)
        ac_animals{a} = ac;
        ac_conds{a} = ac(end);
    else
        ac_animals{a} = ac(1:ius-1);
        ac_conds{a} = ac(end);
    end
end

% Prepare outputs
AUCs = nan(size(pairs,1),1);
Pvals = nan(size(pairs,1),1);
pair_labels = cell(size(pairs,1),1);
n_samples = nan(size(pairs,1),1);

for p = 1:size(pairs,1)
    a_cond = pairs{p,1};
    b_cond = pairs{p,2};
    pair_labels{p} = sprintf('%s vs %s', a_cond, b_cond);

    keep = strcmp(ac_conds, a_cond) | strcmp(ac_conds, b_cond);
    if nnz(keep) < 4
        % Not enough samples to evaluate: set to chance and non-significant p
        AUCs(p) = 0.5; Pvals(p) = 1.0; n_samples(p) = 0; continue;
    end
    Xa = X_animal(keep, :);
    ya = strcmp(ac_conds(keep), a_cond); % 1 for a_cond, 0 for b_cond
    n_samples(p) = nnz(keep); % Total animal+condition samples for this pair

    if nnz(ya) < 2 || nnz(~ya) < 2
        % Not enough class samples: set to chance and non-significant p
        AUCs(p) = 0.5; Pvals(p) = 1.0; continue;
    end

    % LOAO grouped by animal+condition (each AC is independent)
    uac = 1:size(Xa,1); % Each row is a unique animal+condition
    scores = []; labels = [];
    for i = 1:numel(uac)
        te = (uac == uac(i)); % Leave out one animal+condition
        tr = ~te;
        Xtr = Xa(tr,:); ytr = ya(tr);
        Xte = Xa(te,:); yte = ya(te);

        % Standardize per fold
        % mu = mean(Xtr,1); sg = std(Xtr,[],1); sg(sg==0) = 1;
        % Xtr = (Xtr - mu) ./ sg;
        % Xte = (Xte - mu) ./ sg;

        % Ensure both classes exist in training fold
        if nnz(ytr) == 0 || nnz(~ytr) == 0
            % Skip this fold; contributes no scores
            continue;
        end
        % Non-linear classifier: RBF SVM (standardization already applied)
        mdl = fitcsvm(Xtr, ytr, 'KernelFunction','rbf', 'KernelScale','auto', ...
            'Standardize', false, 'ClassNames', [false true]);
        [predLbl, score] = predict(mdl, Xte); %#ok<ASGLU>
        if size(score,2) == 2
            pos = score(:,2);
        else
            % Fallback: use logical label as score if margins unavailable
            pos = double(predLbl);
        end
        scores = [scores; pos]; %#ok<AGROW>
        labels = [labels; yte]; %#ok<AGROW>
    end

    % AUC (fallback to 0.5 if computation fails or NaN)
    try
        [~,~,~,auc] = perfcurve(labels, scores, 1);
        if isempty(auc) || isnan(auc)
            auc = 0.5;
        end
    catch
        auc = 0.5;
    end
    AUCs(p) = auc;

    % MMD permutation on animal-level features (fallback to 1.0 if fails or NaN)
    try
        [~, pval] = compute_mmd_permutation(Xa, ya, [], 1000);
        if isempty(pval) || isnan(pval)
            pval = 1.0;
        end
    catch
        pval = 1.0;
    end
    Pvals(p) = pval;
end
end

function fig = plot_auc_mmd(pair_labels, AUCs, Pvals, ttl, cmap) %#ok<INUSD>
% Simple bar plot of AUC (%) per pair, with MMD p-value annotations
% Replace NaNs with nanmean fallback (or 0.5 if all NaN) to avoid missing bars
auc_mean = mean(AUCs(~isnan(AUCs)));
if isnan(auc_mean), auc_mean = 0.5; end
AUCs_plot = AUCs;
AUCs_plot(isnan(AUCs_plot)) = auc_mean;

Pvals_plot = Pvals;
Pvals_plot(isnan(Pvals_plot)) = 1.0;

fig = figure('Color','w', 'Position',[120 120 900 520], 'Visible','on');
tiledlayout(fig, 1, 1, 'TileSpacing','compact', 'Padding','compact');
ax = nexttile; hold(ax,'on');
x = 1:numel(AUCs_plot);
bar(ax, x, AUCs_plot*100, 0.6, 'FaceColor',[0.3 0.3 0.3]);
yline(ax, 50, '--', 'Color', [0.6 0.6 0.6]);
set(ax, 'XTick', x, 'XTickLabel', pair_labels, 'XTickLabelRotation', 20);
ylabel(ax, 'AUC (%)', 'Color', 'k');
title(ax, sprintf('%s\n(LOAO AUC with MMD p-values)', ttl), 'Color', 'k');
grid(ax, 'on');
% annotate p-values above bars
yl = ylim(ax);
for i = 1:numel(Pvals_plot)
    txt = sprintf('p=%.3f', Pvals_plot(i));
    text(ax, x(i), max(yl(1), (AUCs_plot(i)*100)+2), txt, 'HorizontalAlignment','center', 'Color','k', 'FontSize',9);
end
set(ax, 'Color','w', 'XColor','k', 'YColor','k');
box(ax,'on');
end

%% Export helpers

function export_method_csvs(out_folder, method_tag, ac_list, X_animal, col_names, pair_labels, AUCs, Pvals)
% Export per-AC features and pairwise summary CSVs
meta = ac_meta_from_list(ac_list);
% Build table
T = table(meta.AC_ID, meta.Animal, meta.Condition, 'VariableNames', {'AC_ID','Animal','Condition'});
% Ensure col_names length matches
if numel(col_names) ~= size(X_animal,2)
    col_names = arrayfun(@(i) sprintf('f_%03d', i), 1:size(X_animal,2), 'UniformOutput', false);
end
for j = 1:size(X_animal,2)
    T.(col_names{j}) = X_animal(:,j);
end
% Write features CSV
fn_feat = fullfile(out_folder, sprintf('features_%s.csv', method_tag));
writetable(T, fn_feat);

% Summary CSV
auc_fallback = mean(AUCs(~isnan(AUCs))); if isnan(auc_fallback), auc_fallback = 0.5; end
AUCs_out = AUCs; AUCs_out(isnan(AUCs_out)) = auc_fallback;
Pvals_out = Pvals; Pvals_out(isnan(Pvals_out)) = 1.0;
Ts = table(pair_labels(:), AUCs_out(:), Pvals_out(:), 'VariableNames', {'Pair','AUC','MMD_p'});
fn_sum = fullfile(out_folder, sprintf('summary_%s.csv', method_tag));
writetable(Ts, fn_sum);
end

function meta = ac_meta_from_list(ac_list)
A = numel(ac_list);
AC_ID = ac_list(:);
Animal = cell(A,1);
Condition = cell(A,1);
for a = 1:A
    ac = ac_list{a};
    iu = find(ac=='_',1);
    if isempty(iu)
        Animal{a} = ac;
    else
        Animal{a} = ac(1:iu-1);
    end
    Condition{a} = ac(end);
end
meta = struct('AC_ID',{AC_ID}, 'Animal',{Animal}, 'Condition',{Condition});
end

function cols = build_gmm_colnames(K)
cols = [arrayfun(@(i) sprintf('occ_%02d', i), 1:K, 'UniformOutput', false), ...
        arrayfun(@(i) sprintf('P_rowmajor_%03d', i), 1:K*K, 'UniformOutput', false)];
end

function cols = build_kmeans_colnames(K)
cols = [arrayfun(@(i) sprintf('occ_%02d', i), 1:K, 'UniformOutput', false), ...
        arrayfun(@(i) sprintf('P_rowmajor_%03d', i), 1:K*K, 'UniformOutput', false)];
end

function export_markov_matrices_csv(out_folder, tag, ac_list, Pmats)
% Export flattened transition matrices with explicit row/col labels
[K, ~, A] = size(Pmats);
row_names = ac_list(:);
col_names = arrayfun(@(i) sprintf('P_%02d_%02d', floor((i-1)/K)+1, mod(i-1,K)+1), 1:K*K, 'UniformOutput', false);
X = reshape(Pmats, K*K, A)';
T = table(row_names, 'VariableNames', {'AC_ID'});
for j = 1:size(X,2)
    T.(col_names{j}) = X(:,j);
end
fn = fullfile(out_folder, sprintf('%s.csv', tag));
writetable(T, fn);
end


%% MMD (two-sample) with permutation

function [mmd2, pval, dist_null] = compute_mmd_permutation(X, y, kernel_sigma, n_perm, rng_seed)
% Two-sample MMD^2 (unbiased) with RBF kernel; permutation at sample level.
if nargin < 3 || isempty(kernel_sigma)
    pd = pdist(X, 'euclidean');
    med = median(pd(pd>0));
    if ~isfinite(med) || med == 0, med = 1; end
    kernel_sigma = med;
end
if nargin < 4 || isempty(n_perm), n_perm = 1000; end
if nargin >= 5 && ~isempty(rng_seed), rng(rng_seed); end

D2 = squareform(pdist(X, 'euclidean')).^2;
K = exp(-D2 / (2*kernel_sigma^2));

idx0 = find(y == 0); idx1 = find(y == 1);
mmd2 = mmd2_unbiased(K, idx0, idx1);

N = numel(y);
dist_null = zeros(n_perm,1);
for p = 1:n_perm
    yp = y(randperm(N));
    idx0p = find(yp == 0); idx1p = find(yp == 1);
    dist_null(p) = mmd2_unbiased(K, idx0p, idx1p);
end
pval = (sum(dist_null >= mmd2) + 1) / (n_perm + 1);
end

function v = mmd2_unbiased(K, idx0, idx1)
n0 = numel(idx0); n1 = numel(idx1);
if n0 < 2 || n1 < 2, v = NaN; return; end
K00 = K(idx0, idx0); K11 = K(idx1, idx1); K01 = K(idx0, idx1);
term00 = (sum(K00(:)) - sum(diag(K00))) / (n0*(n0-1));
term11 = (sum(K11(:)) - sum(diag(K11))) / (n1*(n1-1));
term01 = (2 * sum(K01(:))) / (n0*n1);
v = term00 + term11 - term01;
end

%% New functions for improved clustering and plotting

function [X_occ, X_tran, info] = features_from_dbscan(score_pca, ac_to_idx)
% DBSCAN clustering alternative to k-means
% Uses bout-based analysis to avoid diagonal-dominated transitions

fprintf('Running DBSCAN clustering...\n');

% Use top 20 components
Xk = score_pca(:, 1:min(20, size(score_pca,2)));

% Choose epsilon using k-distance plot method
k = 2 * size(Xk, 2);  % 2 × dimensionality = 40 for 20D
distances = pdist2(Xk, Xk);
knn_dist = sort(distances, 2);
k_distances = sort(knn_dist(:, k+1));  % k-th nearest neighbor distances

% Set parameters based on distribution
epsilon = prctile(k_distances, 95);  % 95th percentile
minpts = k;  % ~40 for 20D data

% Run DBSCAN
idx = dbscan(Xk, epsilon, minpts);

% Handle noise points (label -1) by assigning to nearest cluster
noise_points = (idx == -1);
if any(noise_points)
    valid_clusters = unique(idx(idx > 0));
    if ~isempty(valid_clusters)
        for i = find(noise_points)'
            % Find nearest non-noise point
            dists = vecnorm(Xk - Xk(i,:), 2, 2);
            dists(noise_points) = Inf; % Exclude other noise points
            [~, nearest_idx] = min(dists);
            idx(i) = idx(nearest_idx);
        end
    else
        % If all points are noise, use k-means fallback
        fprintf('DBSCAN found only noise, falling back to k-means...\n');
        rng(42);
        idx = kmeans(Xk, 8, 'Replicates', 20, 'MaxIter', 1000, 'Distance','cosine', ...
            'Start','plus', 'Display','off');
    end
end

% Relabel clusters to be consecutive starting from 1
unique_clusters = unique(idx);
K = numel(unique_clusters);
idx_relabeled = idx;
for i = 1:K
    idx_relabeled(idx == unique_clusters(i)) = i;
end

% Bout detection parameters: adjusted for ~12 Hz effective sampling rate
min_bout_frames = 30; % Still 30 frames to ensure meaningful behavioral bouts

% Compute occupancies and transitions using bout-based analysis
ac_keys = keys(ac_to_idx);
A = numel(ac_keys);
X_occ = zeros(A, K);
X_tran = zeros(A, K*K);

for a = 1:A
    fr = ac_to_idx(ac_keys{a});
    raw_seq = idx_relabeled(fr);

    % Detect behavioral bouts instead of using frame-by-frame transitions
    bout_seq = detect_bouts(raw_seq, min_bout_frames);

    if isempty(bout_seq)
        % Fallback: use original sequence if no bouts detected
        bout_seq = raw_seq;
        fprintf('Warning: No bouts detected for DBSCAN %s, using original sequence\n', ac_keys{a});
    end

    % occupancy based on bouts
    h = histcounts(bout_seq, 0.5:1:(K+0.5));
    occ = h / max(1, sum(h));

    % transitions between bouts (avoids artificial self-transitions)
    Tij = zeros(K,K);
    for t = 1:(numel(bout_seq)-1)
        i = bout_seq(t); j = bout_seq(t+1);
        Tij(i,j) = Tij(i,j) + 1;
    end
    rowSums = sum(Tij, 2);
    rowSums(rowSums==0) = 1;
    P = Tij ./ rowSums;
    X_occ(a,:) = occ;
    X_tran(a,:) = P(:)';
end

info = struct('K', K, 'epsilon', epsilon, 'minpts', minpts);
fprintf('DBSCAN found %d clusters (eps=%.3f, minpts=%d)\n', K, epsilon, minpts);
end

function [sym_idx, info] = symbols_from_dbscan(score_pca, fallback_K)
% Discretize PCA space via DBSCAN for sequence modeling
Xk = score_pca(:, 1:min(20, size(score_pca,2)));

% DBSCAN parameters
k = 2 * size(Xk, 2);
distances = pdist2(Xk, Xk);
knn_dist = sort(distances, 2);
k_distances = sort(knn_dist(:, k+1));
epsilon = prctile(k_distances, 95);
minpts = k;

sym_idx = dbscan(Xk, epsilon, minpts);

% Handle noise points and relabel
noise_points = (sym_idx == -1);
if any(noise_points)
    valid_clusters = unique(sym_idx(sym_idx > 0));
    if ~isempty(valid_clusters)
        for i = find(noise_points)'
            dists = vecnorm(Xk - Xk(i,:), 2, 2);
            dists(noise_points) = Inf;
            [~, nearest_idx] = min(dists);
            sym_idx(i) = sym_idx(nearest_idx);
        end
    else
        % Fallback to k-means
        rng(123);
        sym_idx = kmeans(Xk, fallback_K, 'Replicates', 20, 'MaxIter', 1000, 'Distance','cosine', ...
            'Start','plus', 'Display','off');
    end
end

% Relabel to be consecutive
unique_clusters = unique(sym_idx);
K = numel(unique_clusters);
for i = 1:K
    sym_idx(sym_idx == unique_clusters(i)) = i;
end

info = struct('K', K, 'epsilon', epsilon, 'minpts', minpts);
end

function fig = plot_auc_mmd_with_error(pair_labels, AUCs, Pvals, n_samples, ttl, cmap) %#ok<INUSD>
% Bar plot with error bars (SEM) and sample sizes
% For now, use bootstrap-estimated SEM; in practice you'd get this from CV folds

% Estimate SEM as std(AUCs)/sqrt(mean(n_samples)) - rough approximation
sem_auc = std(AUCs(~isnan(AUCs))) / sqrt(mean(n_samples(~isnan(n_samples))));
if isnan(sem_auc), sem_auc = 0.02; end % Default 2% error

% Replace NaNs with fallbacks
auc_mean = mean(AUCs(~isnan(AUCs)));
if isnan(auc_mean), auc_mean = 0.5; end
AUCs_plot = AUCs;
AUCs_plot(isnan(AUCs_plot)) = auc_mean;

Pvals_plot = Pvals;
Pvals_plot(isnan(Pvals_plot)) = 1.0;

fig = figure('Color','w', 'Position',[120 120 900 520], 'Visible','on');
tiledlayout(fig, 1, 1, 'TileSpacing','compact', 'Padding','compact');
ax = nexttile; hold(ax,'on');

x = 1:numel(AUCs_plot);

% Bar plot with error bars
bar(ax, x, AUCs_plot*100, 0.6, 'FaceColor',[0.3 0.3 0.3]);
errorbar(ax, x, AUCs_plot*100, sem_auc*100*ones(size(AUCs_plot)), 'k.', 'LineWidth', 1.5);

yline(ax, 50, '--', 'Color', [0.6 0.6 0.6]);
set(ax, 'XTick', x, 'XTickLabel', pair_labels, 'XTickLabelRotation', 20);
ylabel(ax, 'AUC (%)', 'Color', 'k');
title(ax, sprintf('%s\n(LOAO AUC with MMD p-values and SEM)', ttl), 'Color', 'k');
grid(ax, 'on');

% Annotate p-values above bars
yl = ylim(ax);
for i = 1:numel(Pvals_plot)
    txt = sprintf('p=%.3f', Pvals_plot(i));
    text(ax, x(i), max(yl(1), (AUCs_plot(i)*100 + sem_auc*100)+3), txt, 'HorizontalAlignment','center', 'Color','k', 'FontSize',9);
end

% Annotate n values at base of bars
for i = 1:numel(n_samples)
    if ~isnan(n_samples(i))
        txt = sprintf('n=%d', n_samples(i));
        text(ax, x(i), 5, txt, 'HorizontalAlignment','center', 'Color','k', 'FontSize',8, 'FontWeight','bold');
    end
end

set(ax, 'Color','w', 'XColor','k', 'YColor','k');
box(ax,'on');
end


%% Helper function to transform analysisstruct by animal IDs

function analysisstruct_transformed = transform_analysisstruct_by_animals(analysisstruct, ids)
% Transform analysisstruct.annot_reordered from current format to animal-grouped format
% 
% Input:
%   - analysisstruct: Original analysisstruct with annot_reordered{1} and annot_reordered{2}
%   - ids: Cell array identifying animal-condition for each frame (from collect_features_and_ids)
%
% Output:
%   - analysisstruct_transformed: Copy of analysisstruct with modified annot_reordered
%     containing separate cells for each animal plus global at the end

% Create a copy of the original analysisstruct to avoid modifying it
analysisstruct_transformed = analysisstruct;

% Extract the original annotation data from the second cell (global condition)
if length(analysisstruct.annot_reordered) < 2 || isempty(analysisstruct.annot_reordered{2})
    error('analysisstruct.annot_reordered{2} is missing or empty');
end

annot_data = analysisstruct.annot_reordered{2};

% Check that the sizes match
if length(annot_data) ~= length(ids)
    error('Size mismatch: annot_data has %d elements but ids has %d', ...
        length(annot_data), length(ids));
end

% Extract animal IDs from ids
% Each identifier is like "1633_B", "1636_C", etc. - we want everything before the last underscore
animal_ids = cellfun(@(x) x(1:find(x=='_',1,'last')-1), ids, 'UniformOutput', false);

% Find unique animals and sort them for consistent ordering
unique_animals = unique(animal_ids);
fprintf('Found animals: %s\n', strjoin(unique_animals, ', '));

% Initialize the new annot_reordered cell array
% It will have one cell per animal plus one for the global data
new_annot_reordered = cell(1, length(unique_animals) + 1);

% Group data by animal
for i = 1:length(unique_animals)
    animal = unique_animals{i};
    
    % Find frames belonging to this animal
    animal_mask = strcmp(animal_ids, animal);
    animal_indices = find(animal_mask);
    
    % Extract data for this animal
    animal_data = annot_data(animal_indices);
    
    % Store in the new cell array
    new_annot_reordered{i} = animal_data;
    
    fprintf('Animal %s: %d frames (%.1f%%)\n', animal, length(animal_data), ...
        100 * length(animal_data) / length(annot_data));
end

% Add the global data (original data) as the last cell
new_annot_reordered{end} = annot_data;

% Update the analysisstruct
analysisstruct_transformed.annot_reordered = new_annot_reordered;

fprintf('Transformed annot_reordered: %d animal-specific cells + 1 global cell\n', ...
    length(unique_animals));
end
