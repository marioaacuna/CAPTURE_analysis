function run_unbiased_state_space()
% Unbiased animal-level behavioral discriminability without t-SNE.
% Methods:
%   1) PCA + GMM (states = soft posteriors) -> occupancy + soft transitions
%   2) PCA + k-means (hard states)         -> occupancy + transitions
%   3) No clustering (aggregate stats per feature) -> mean/std/quantiles
%   4) Sequence model (discrete HMM on k-means symbols; fallback Markov)
%
% Evaluation:
%   - Leave-one-animal-out (grouped by animal ID) logistic AUC
%   - MMD permutation test (animal-level)
%
% Outputs:
%   - Figures per method with AUC and MMD for pairs: F vs B, F vs S/H,
%     N vs B, N vs H, F vs C, N vs G.
%   - Saved under GC.temp_root/figs_unbiased_state_space/<exp.name>/

clc;
GC = general_configs();

% Match your existing experiment list
experiments = struct();
experiments(1).name = 'BSFC_300hz';
experiments(1).folder = '0_preprocessing_BSFC_300hz';
experiments(2).name = 'BHNG_300hz';
experiments(2).folder = '0_preprocessing_BHNG_300hz';

for exp_idx = 1:numel(experiments)
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

    % Preprocess: impute NaNs, standardize, PCA
    [X_std, col_mu, col_sig] = standardize_impute(X_raw); %#ok<NASGU>
    ncomp = min(50, size(X_std,2));
    [score_pca, pca_info] = do_pca(X_std, ncomp); %#ok<NASGU>

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
            'F','B';
            'F',ctrlPain;
            'F','C';
        };
    elseif contains(exp.name, 'BHNG')
        candidate_pairs = {
            'N','B';
            'N',ctrlPain;
            'N','G';
        };
    else
        % Fallback: anchor to present main condition
        if any(strcmp(uconds,'F'))
            candidate_pairs = {
                'F','B';
                'F',ctrlPain;
                'F','C';
            };
        else
            candidate_pairs = {
                'N','B';
                'N',ctrlPain;
                'N','G';
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
        K_candidates = [8 12 16 24 32];
        [X_occ_tr, X_tran_tr, info_gmm] = features_from_gmm(score_pca, ac_to_idx, K_candidates);
        X_animal_gmm = [X_occ_tr, X_tran_tr];
        [AUCs_gmm, Pvals_gmm, pair_labels] = loao_and_mmd(X_animal_gmm, ac_list, animals, conds, pairs);
        fig1 = plot_auc_mmd(pair_labels, AUCs_gmm, Pvals_gmm, sprintf('%s - PCA+GMM (K=%d, %dD PCA)', exp.name, info_gmm.K, size(score_pca,2)), cmap);
        exportgraphics(fig1, fullfile(out_folder, 'unbiased_pca_gmm_auc_mmd.pdf'), 'ContentType','vector', 'BackgroundColor','white');
        close(fig1);
        % Export CSVs
        export_method_csvs(out_folder, 'pca_gmm', ac_list, X_animal_gmm, build_gmm_colnames(info_gmm.K), pair_labels, AUCs_gmm, Pvals_gmm);
    catch ME
        %warning('PCA+GMM failed: %s', ME.message);
    end

    % 2) PCA + k-means
    fprintf('Method 2: PCA + k-means (hard states)...\n');
    try
        K_km = 24; % reasonable default, can be tuned
        [X_occ_km, X_tran_km, info_km] = features_from_kmeans(score_pca, ac_to_idx, K_km);
        X_animal_km = [X_occ_km, X_tran_km];
        [AUCs_km, Pvals_km, pair_labels] = loao_and_mmd(X_animal_km, ac_list, animals, conds, pairs);
        fig2 = plot_auc_mmd(pair_labels, AUCs_km, Pvals_km, sprintf('%s - PCA+k-means (K=%d, %dD PCA)', exp.name, info_km.K, size(score_pca,2)), cmap);
        exportgraphics(fig2, fullfile(out_folder, 'unbiased_pca_kmeans_auc_mmd.pdf'), 'ContentType','vector', 'BackgroundColor','white');
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
        [AUCs_agg, Pvals_agg, pair_labels] = loao_and_mmd(X_animal_agg, ac_list, animals, conds, pairs);
        fig3 = plot_auc_mmd(pair_labels, AUCs_agg, Pvals_agg, sprintf('%s - Aggregates (means/std/quantiles)', exp.name), cmap);
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
    fprintf('Method 4: Sequence model (HMM discrete on k-means symbols)...\n');
    try
    K_sym = 24;
    [sym_idx, ~] = symbols_from_kmeans(score_pca, K_sym);
        H_hidden = 12; % hidden states for HMM
        [X_animal_hmm, used_hmm] = features_from_hmm_discrete(sym_idx, ac_to_idx, K_sym, H_hidden);
        [AUCs_hmm, Pvals_hmm, pair_labels] = loao_and_mmd(X_animal_hmm, ac_list, animals, conds, pairs);
        tag = 'HMM'; if ~used_hmm, tag = 'Markov fallback'; end
        fig4 = plot_auc_mmd(pair_labels, AUCs_hmm, Pvals_hmm, sprintf('%s - Sequence model (%s, Ksym=%d, H=%d)', exp.name, tag, K_sym, H_hidden), cmap);
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
        K_sym_mk = 24;
    [sym_idx_mk, ~] = symbols_from_kmeans(score_pca, K_sym_mk);
    [X_markov, ~, Pmats_markov] = features_markov_only(sym_idx_mk, ac_to_idx, K_sym_mk);
        [AUCs_mk, Pvals_mk, pair_labels] = loao_and_mmd(X_markov, ac_list, animals, conds, pairs);
        fig5 = plot_auc_mmd(pair_labels, AUCs_mk, Pvals_mk, sprintf('%s - Markov chain only (Ksym=%d)', exp.name, K_sym_mk), cmap);
        exportgraphics(fig5, fullfile(out_folder, 'unbiased_markov_only_auc_mmd.pdf'), 'ContentType','vector', 'BackgroundColor','white');
        close(fig5);
        % Export CSVs with explicit Markov column names
        markov_cols = [arrayfun(@(i) sprintf('occ_%02d', i), 1:K_sym_mk, 'UniformOutput', false), ...
                       arrayfun(@(i) sprintf('P_rowmajor_%03d', i), 1:K_sym_mk*K_sym_mk, 'UniformOutput', false)];
        export_method_csvs(out_folder, 'markov_only', ac_list, X_markov, markov_cols, pair_labels, AUCs_mk, Pvals_mk);
        % Optional: export full transition matrices per AC as separate CSV
        export_markov_matrices_csv(out_folder, 'markov_transition_matrices', ac_list, Pmats_markov);
    catch ME
        %warning('Markov-only analysis failed: %s', ME.message);
    end

    fprintf('Saved results to %s\n', out_folder);
end

fprintf('\nAll experiments complete.\n');
end


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
rng(42);
[idx, ~] = kmeans(score_pca, K, 'Replicates', 10, 'MaxIter', 300, 'Display', 'off');
ac_keys = keys(ac_to_idx);
A = numel(ac_keys);
X_occ = zeros(A, K);
X_tran = zeros(A, K*K);

for a = 1:A
    fr = ac_to_idx(ac_keys{a});
    seq = idx(fr);
    % occupancy
    h = histcounts(seq, 0.5:1:(K+0.5));
    occ = h / max(1, sum(h));
    % transitions (hard)
    Tij = zeros(K,K);
    for t = 1:(numel(seq)-1)
        i = seq(t); j = seq(t+1);
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
[sym_idx, ~] = kmeans(score_pca, K, 'Replicates', 10, 'MaxIter', 300, 'Display', 'off');
info = struct('K', K);
end

function [X_markov, occ_all, Pmats] = features_markov_only(sym_idx, ac_to_idx, Ksym)
% Build per-AC Markov occupancy and transition features from discrete symbols
ac_keys = keys(ac_to_idx);
A = numel(ac_keys);
occ_all = zeros(A, Ksym);
Pmats = zeros(Ksym, Ksym, A);
for a = 1:A
    seq = sym_idx(ac_to_idx(ac_keys{a}));
    % occupancy
    h = histcounts(seq, 0.5:1:(Ksym+0.5));
    occ = h / max(1, sum(h));
    occ_all(a,:) = occ;
    % transitions
    Tij = zeros(Ksym, Ksym);
    for t = 1:(numel(seq)-1)
        i = seq(t); j = seq(t+1);
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
% If hmmtrain is unavailable, fallback to per-AC Markov chain features.
used_hmm = false;
ac_keys = keys(ac_to_idx);
A = numel(ac_keys);

if exist('hmmtrain', 'file') == 2 && exist('hmmdecode', 'file') == 2 && exist('hmmviterbi', 'file') == 2
    % Prepare sequences (as cell array)
    seqs = cell(A,1);
    for a = 1:A
        seqs{a} = sym_idx(ac_to_idx(ac_keys{a}))';
        if isempty(seqs{a}), seqs{a} = 1; end
    end

    % Initialize HMM parameters
    rng(7);
    Atrans = mkstochastic(rand(H_hidden, H_hidden));
    Emiss = mkstochastic(rand(H_hidden, Ksym));

    try
        [ESTTR, ESTEMIT] = hmmtrain(seqs, Atrans, Emiss, 'Maxiterations', 50, 'Verbose', false);
        used_hmm = true;
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
            Tij(q(t), q(t+1)) = Tij(q(t), q(t+1)) + 1;
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
    X_occ = zeros(A, Ksym);
    X_tran = zeros(A, Ksym*Ksym);
    for a = 1:A
        seq = sym_idx(ac_to_idx(ac_keys{a}));
        h = histcounts(seq, 0.5:1:(Ksym+0.5));
        occ = h / max(1, sum(h));
        Tij = zeros(Ksym, Ksym);
        for t = 1:(numel(seq)-1)
            i = seq(t); j = seq(t+1);
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


%% Evaluation and plotting

function [AUCs, Pvals, pair_labels] = loao_and_mmd(X_animal, ac_list, ~, ~, pairs)
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

for p = 1:size(pairs,1)
    a_cond = pairs{p,1};
    b_cond = pairs{p,2};
    pair_labels{p} = sprintf('%s vs %s', a_cond, b_cond);

    keep = strcmp(ac_conds, a_cond) | strcmp(ac_conds, b_cond);
    if nnz(keep) < 4
        % Not enough samples to evaluate: set to chance and non-significant p
        AUCs(p) = 0.5; Pvals(p) = 1.0; continue;
    end
    Xa = X_animal(keep, :);
    ya = strcmp(ac_conds(keep), a_cond); % 1 for a_cond, 0 for b_cond
    anims = ac_animals(keep);

    if nnz(ya) < 2 || nnz(~ya) < 2
        % Not enough class samples: set to chance and non-significant p
        AUCs(p) = 0.5; Pvals(p) = 1.0; continue;
    end

    % LOAO grouped by animal
    uanim = unique(anims);
    scores = []; labels = [];
    for i = 1:numel(uanim)
        te = strcmp(anims, uanim{i});
        tr = ~te;
        Xtr = Xa(tr,:); ytr = ya(tr);
        Xte = Xa(te,:); yte = ya(te);

        % Standardize per fold
        mu = mean(Xtr,1); sg = std(Xtr,[],1); sg(sg==0) = 1;
        Xtr = (Xtr - mu) ./ sg;
        Xte = (Xte - mu) ./ sg;

        mdl = fitclinear(Xtr, ytr, 'Learner','logistic', 'Regularization','lasso', 'Solver','sparsa', ...
            'GradientTolerance',1e-6, 'BetaTolerance',1e-6, 'Lambda', 'auto');
        [~,score] = predict(mdl, Xte);
        if size(score,2) == 2
            pos = score(:,2);
        else
            % If predict gives only labels, fall back to distance approximation
            pos = double(predict(mdl, Xte));
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