%% Run individuality analysis summary across multiple zValues fields per experiment
% Loads previously saved analysis_<zfield>.mat files (no re-clustering),
% computes individuality metrics for each zfield, saves metrics, and
% generates a summary plot across zfields emphasizing F–C (BSFC) and N–G (BHNG).

clear; clc;
GC = general_configs();

% Experiments and folders (must match where data were saved)
experiments = struct();
experiments(1).name = 'BSFC_300hz';
experiments(1).folder = '0_preprocessing_BSFC_300hz';
experiments(2).name = 'BHNG_300hz';
experiments(2).folder = '0_preprocessing_BHNG_300hz';

% zValues variants to process (must match saved analysis_<zfield>.mat files)
zvals_fields =  {'zValues_all_features', ...
                 'zValues_important_features', ...
                 'zValues_jt_features_only', ...
                 'zValues_hand_made', ...
                 'zValues_only_important'};

for exp_idx = 1:numel(experiments)
    exp = experiments(exp_idx);
    fprintf('\n=== Experiment: %s ===\n', exp.name);

    data_folder = fullfile(GC.project_path, 'data', exp.folder);
    if ~exist(data_folder, 'dir')
        warning('Data folder not found: %s', data_folder);
        continue;
    end

    % Load animal identifiers (required for individuality)
    fn_predictions = fullfile(data_folder, 'agg_predictions.mat');
    if ~exist(fn_predictions, 'file')
        warning('Predictions file not found: %s', fn_predictions);
        continue;
    end
    S_pred = load(fn_predictions, 'animal_condition_identifier');
    if ~isfield(S_pred, 'animal_condition_identifier')
        warning('animal_condition_identifier var missing in %s', fn_predictions);
        continue;
    end

    % Output folder (same base as clustering outputs)
    exp_out_folder = fullfile(GC.temp_root, 'clusters_many_zvals', exp.name);
    if ~exist(exp_out_folder, 'dir'); mkdir(exp_out_folder); end

    indiv_results = struct('zfield', {}, 'metrics', {});

    % Process each saved analysis_<zfield>.mat
    for zi = 1:numel(zvals_fields)
        zf = zvals_fields{zi};
        fn_analysis = fullfile(data_folder, sprintf('analysis_%s.mat', zf));
        if ~exist(fn_analysis, 'file')
            fprintf('Skipping %s (missing %s)\n', zf, fn_analysis);
            continue;
        end

        S = load(fn_analysis, 'analysisstruct');
        if ~isfield(S, 'analysisstruct')
            warning('analysisstruct missing in %s', fn_analysis);
            continue;
        end
        analysisstruct = S.analysisstruct; %#ok<NASGU>

        try
            metrics = compute_individuality_analysis(S.analysisstruct, S_pred.animal_condition_identifier, GC);
        catch ME
            warning('Individuality computation failed for %s: %s', zf, ME.message);
            continue;
        end

        % Save metrics under per-zfield folder
        zf_folder = fullfile(exp_out_folder, zf);
        if ~exist(zf_folder, 'dir'); mkdir(zf_folder); end
        save(fullfile(zf_folder, sprintf('individuality_%s.mat', zf)), 'metrics');

        % Export CSVs for downstream statistics for this zval category
        try
            export_individuality_csvs(zf_folder, zf, metrics);
        catch ME
            warning('CSV export failed for %s: %s', zf, ME.message);
        end

        % Collect for summary
        indiv_results(end+1).zfield = zf; %#ok<SAGROW>
        indiv_results(end).metrics = metrics;

        fprintf('Done %s\n', zf);
    end

    % Summary across zfields
    if ~isempty(indiv_results)
        try
            plot_individuality_summary_across_zvals(exp.name, indiv_results, exp_out_folder);
            fprintf('Saved summary plot for %s\n', exp.name);
        catch ME
            warning('Summary plotting failed for %s: %s', exp.name, ME.message);
        end
    else
        fprintf('No zfields processed for %s.\n', exp.name);
    end
end

fprintf('\nAll experiments complete.\n');


function metrics = compute_individuality_analysis(analysisstruct, animal_condition_identifier, GC)
% Compute individuality metrics using GLOBAL data (all conditions, all animals)
% Returns a struct with cluster- and frame-level metrics.

individuality_threshold = 0.8; % 80%

% Build global animal+condition ids aligned to frames used in analysis
upsampling_factor = GC.repfactor;
long_ids = repelem(animal_condition_identifier, upsampling_factor);
if iscell(analysisstruct.frames_with_good_tracking)
    good_idx = analysisstruct.frames_with_good_tracking{1};
else
    good_idx = analysisstruct.frames_with_good_tracking;
end

function export_individuality_csvs(out_folder, zf, metrics)
% Create per-AC, per-cluster, and per-condition CSVs capturing key metrics

% 1) Per animal+condition metrics
keys_total = keys(metrics.animal_total_frames);
A = numel(keys_total);
AC_ID = cell(A,1); Animal = cell(A,1); Condition = cell(A,1);
total_frames = zeros(A,1); individual_frames = zeros(A,1); frame_indiv_pct = nan(A,1);
total_clusters_present = zeros(A,1); total_dominant_clusters = zeros(A,1); individual_clusters = zeros(A,1);
for ii = 1:A
    ac = keys_total{ii};
    AC_ID{ii} = ac;
    us = find(ac=='_',1);
    if isempty(us), Animal{ii} = ac; else, Animal{ii} = ac(1:us-1); end
    Condition{ii} = ac(end);
    total_frames(ii) = metrics.animal_total_frames(ac);
    if isKey(metrics.animal_individual_frames, ac)
        individual_frames(ii) = metrics.animal_individual_frames(ac);
    else
        individual_frames(ii) = 0;
    end
    if total_frames(ii) > 0
        frame_indiv_pct(ii) = individual_frames(ii) / total_frames(ii);
    else
        frame_indiv_pct(ii) = NaN;
    end
    if isKey(metrics.animal_total_clusters_present, ac)
        total_clusters_present(ii) = metrics.animal_total_clusters_present(ac);
    else
        total_clusters_present(ii) = 0;
    end
    if isKey(metrics.animal_total_dominant_clusters, ac)
        total_dominant_clusters(ii) = metrics.animal_total_dominant_clusters(ac);
    else
        total_dominant_clusters(ii) = 0;
    end
    if isKey(metrics.animal_individual_clusters, ac)
        individual_clusters(ii) = metrics.animal_individual_clusters(ac);
    else
        individual_clusters(ii) = 0;
    end
end
T_ac = table(AC_ID, Animal, Condition, total_frames, individual_frames, frame_indiv_pct, ...
             total_clusters_present, total_dominant_clusters, individual_clusters);
fn1 = fullfile(out_folder, sprintf('per_ac_%s.csv', zf));
writetable(T_ac, fn1);

% 2) Per-cluster composition and dominance
clus_ids = metrics.cluster_ids(:);
Nclus = numel(clus_ids);
dom_ac = metrics.cluster_dominant_animal(:);
dom_pct = metrics.cluster_dominance_percentage(:);
is_indiv = metrics.cluster_is_individual(:);

% Collect all AC ids appearing in any cluster composition
comp_all_AC = {};
for cc = 1:Nclus
    comp = metrics.cluster_animal_composition{cc};
    if isempty(comp), continue; end
    fns = fieldnames(comp);
    for jj = 1:numel(fns)
        ac_clean = regexprep(fns{jj}(4:end), '__', '_');
        if ~any(strcmp(comp_all_AC, ac_clean))
            comp_all_AC{end+1} = ac_clean; %#ok<AGROW>
        end
    end
end

% Build table with dynamic columns per AC
Cluster_ID = clus_ids;
Dominant_AC = cell(Nclus,1); Dominant_Animal = cell(Nclus,1); Dominant_Condition = cell(Nclus,1);
Dominance_Pct = dom_pct; Is_Individual = is_indiv;
Comp = zeros(Nclus, numel(comp_all_AC));
for cc = 1:Nclus
    Dominant_AC{cc} = dom_ac{cc};
    if strcmp(Dominant_AC{cc}, 'None')
        Dominant_Animal{cc} = '';
        Dominant_Condition{cc} = '';
    else
        us = find(Dominant_AC{cc}=='_',1);
        if isempty(us), Dominant_Animal{cc} = Dominant_AC{cc}; else, Dominant_Animal{cc} = Dominant_AC{cc}(1:us-1); end
        Dominant_Condition{cc} = Dominant_AC{cc}(end);
    end
    comp = metrics.cluster_animal_composition{cc};
    if ~isempty(comp)
        fns = fieldnames(comp);
        for jj = 1:numel(fns)
            ac_clean = regexprep(fns{jj}(4:end), '__', '_');
            kk = find(strcmp(comp_all_AC, ac_clean));
            if ~isempty(kk)
                Comp(cc,kk) = comp.(fns{jj});
            end
        end
    end
end
T_clus = table(Cluster_ID, Dominant_AC, Dominant_Animal, Dominant_Condition, Dominance_Pct, Is_Individual);
% Append composition columns
for jj = 1:numel(comp_all_AC)
    colname = ['Comp_' regexprep(comp_all_AC{jj}, '[^a-zA-Z0-9_]', '_')];
    T_clus.(colname) = Comp(:,jj);
end
fn2 = fullfile(out_folder, sprintf('clusters_%s.csv', zf));
writetable(T_clus, fn2);

% 3) Per-condition aggregates for this zval (cluster share and frame-level)
% Infer conditions present from AC IDs
conds = unique(cellfun(@(ac) ac(end), keys_total, 'UniformOutput', false));
conds = conds(:)';
cluster_share = nan(numel(conds),1);
cond_total_frames = zeros(numel(conds),1);
cond_ind_frames = zeros(numel(conds),1);
for ci = 1:numel(conds)
    cc = conds{ci};
    % frame aggregates
    for ii = 1:A
        if AC_ID{ii}(end) == cc
            cond_total_frames(ci) = cond_total_frames(ci) + total_frames(ii);
            cond_ind_frames(ci) = cond_ind_frames(ci) + individual_frames(ii);
        end
    end
    % cluster share among individual clusters
    indiv_mask = is_indiv;
    if any(indiv_mask)
        dom_cond = cellfun(@(s) s(end), dom_ac(indiv_mask), 'UniformOutput', false);
        cluster_share(ci) = sum(strcmp(dom_cond, cc)) / sum(indiv_mask);
    else
        cluster_share(ci) = NaN;
    end
end
frame_pct = zeros(numel(conds),1);
for ci = 1:numel(conds)
    if cond_total_frames(ci) > 0
        frame_pct(ci) = cond_ind_frames(ci) / cond_total_frames(ci);
    else
        frame_pct(ci) = NaN;
    end
end
Condition = conds(:);
Cluster_Share = cluster_share(:);
Frame_Indiv_Pct = frame_pct(:);
Total_Frames = cond_total_frames(:);
Individual_Frames = cond_ind_frames(:);
Overall_Frame_Indiv_Pct = repmat(metrics.overall_frame_individuality_pct, numel(Condition), 1);
T_cond = table(Condition, Cluster_Share, Frame_Indiv_Pct, Total_Frames, Individual_Frames, Overall_Frame_Indiv_Pct);
fn3 = fullfile(out_folder, sprintf('by_condition_%s.csv', zf));
writetable(T_cond, fn3);
end
if islogical(good_idx)
    good_idx = find(good_idx);
end
global_animal_condition_ids = long_ids(good_idx);

% Cluster assignments (global, final annotation)
cluster_assignments = analysisstruct.annot_reordered{end,end};
cluster_assignments = cluster_assignments(:);
if numel(cluster_assignments) ~= numel(global_animal_condition_ids)
    error('Mismatch between cluster_assignments (%d) and global_animal_condition_ids (%d).', ...
        numel(cluster_assignments), numel(global_animal_condition_ids));
end

unique_clusters = unique(cluster_assignments);
unique_clusters = unique_clusters(unique_clusters > 0); % drop background/noise if 0

global_animal_names = cellfun(@(x) x(1:find(x=='_',1)-1), global_animal_condition_ids, 'UniformOutput', false);

% Cluster-level composition and dominance
cluster_animal_composition = cell(length(unique_clusters), 1);
cluster_dominant_animal = cell(length(unique_clusters), 1);
cluster_dominance_percentage = zeros(length(unique_clusters), 1);
cluster_is_individual = false(length(unique_clusters), 1);

for c_idx = 1:length(unique_clusters)
    cid = unique_clusters(c_idx);
    c_frames = (cluster_assignments == cid);
    total_frames_in_cluster = sum(c_frames);

    if total_frames_in_cluster == 0
        cluster_dominant_animal{c_idx} = 'None';
        cluster_dominance_percentage(c_idx) = 0;
        cluster_is_individual(c_idx) = false;
        cluster_animal_composition{c_idx} = struct();
        continue;
    end

    counts_map = containers.Map();
    for i = 1:length(global_animal_condition_ids)
        if c_frames(i)
            key = global_animal_condition_ids{i};
            if isKey(counts_map, key)
                counts_map(key) = counts_map(key) + 1;
            else
                counts_map(key) = 1;
            end
        end
    end

    if isempty(keys(counts_map))
        cluster_dominant_animal{c_idx} = 'None';
        cluster_dominance_percentage(c_idx) = 0;
        cluster_is_individual(c_idx) = false;
        cluster_animal_composition{c_idx} = struct();
        continue;
    end

    ids = keys(counts_map);
    vals = cell2mat(values(counts_map));
    [max_count, mxi] = max(vals);
    dominant_id = ids{mxi};
    dominance = max_count / total_frames_in_cluster;

    cluster_dominant_animal{c_idx} = dominant_id;
    cluster_dominance_percentage(c_idx) = dominance;
    cluster_is_individual(c_idx) = dominance >= individuality_threshold;

    comp = struct();
    for j = 1:numel(ids)
        valid_field_name = ['ID_' regexprep(ids{j}, '[^a-zA-Z0-9_]', '_')];
        comp.(valid_field_name) = vals(j) / total_frames_in_cluster;
    end
    cluster_animal_composition{c_idx} = comp;
end

% Per animal+condition cluster stats
animal_condition_individual_clusters = containers.Map();
animal_condition_total_clusters_present = containers.Map();
animal_condition_total_dominant_clusters = containers.Map();

for c_idx = 1:length(unique_clusters)
    comp = cluster_animal_composition{c_idx};
    if ~isempty(comp)
        fn = fieldnames(comp);
        for f = 1:numel(fn)
            ac_id = fn{f}(4:end); % strip 'ID_'
            ac_id = regexprep(ac_id, '__', '_');
            if isKey(animal_condition_total_clusters_present, ac_id)
                animal_condition_total_clusters_present(ac_id) = animal_condition_total_clusters_present(ac_id) + 1;
            else
                animal_condition_total_clusters_present(ac_id) = 1;
            end
        end
    end

    dom = cluster_dominant_animal{c_idx};
    if ~strcmp(dom, 'None')
        if isKey(animal_condition_total_dominant_clusters, dom)
            animal_condition_total_dominant_clusters(dom) = animal_condition_total_dominant_clusters(dom) + 1;
        else
            animal_condition_total_dominant_clusters(dom) = 1;
        end
        if cluster_is_individual(c_idx)
            if isKey(animal_condition_individual_clusters, dom)
                animal_condition_individual_clusters(dom) = animal_condition_individual_clusters(dom) + 1;
            else
                animal_condition_individual_clusters(dom) = 1;
            end
        end
    end
end

% Frame-level stats
animal_condition_individual_frames = containers.Map();
animal_condition_total_frames = containers.Map();

for i = 1:length(global_animal_condition_ids)
    ac_id = global_animal_condition_ids{i};
    if isKey(animal_condition_total_frames, ac_id)
        animal_condition_total_frames(ac_id) = animal_condition_total_frames(ac_id) + 1;
    else
        animal_condition_total_frames(ac_id) = 1;
    end
end

for c_idx = 1:length(unique_clusters)
    if cluster_is_individual(c_idx) && ~strcmp(cluster_dominant_animal{c_idx}, 'None')
        cid = unique_clusters(c_idx);
        dom_ac = cluster_dominant_animal{c_idx};
        c_frames = (cluster_assignments == cid);
        frames_for_dom = 0;
        for i = 1:length(global_animal_condition_ids)
            if c_frames(i) && strcmp(global_animal_condition_ids{i}, dom_ac)
                frames_for_dom = frames_for_dom + 1;
            end
        end
        if isKey(animal_condition_individual_frames, dom_ac)
            animal_condition_individual_frames(dom_ac) = animal_condition_individual_frames(dom_ac) + frames_for_dom;
        else
            animal_condition_individual_frames(dom_ac) = frames_for_dom;
        end
    end
end

% Compute per-animal+condition frame individuality percentage and overall summary
frame_indiv_pct_map = containers.Map();
keys_total = keys(animal_condition_total_frames);
for k = 1:numel(keys_total)
    key = keys_total{k};
    tot = animal_condition_total_frames(key);
    ind = 0;
    if isKey(animal_condition_individual_frames, key)
        ind = animal_condition_individual_frames(key);
    end
    if tot > 0
        frame_indiv_pct_map(key) = ind / tot;
    else
        frame_indiv_pct_map(key) = NaN;
    end
end

vals_total = values(animal_condition_total_frames);
total_global_frames = sum(cell2mat(vals_total));
if isempty(keys(animal_condition_individual_frames))
    total_individual_frames = 0;
else
    vals_individual = values(animal_condition_individual_frames);
    total_individual_frames = sum(cell2mat(vals_individual));
end
overall_frame_indiv_pct = total_individual_frames / max(total_global_frames, 1);

% Pack results
metrics = struct();
metrics.threshold = individuality_threshold;
metrics.cluster_ids = unique_clusters;
metrics.cluster_animal_composition = cluster_animal_composition;
metrics.cluster_dominant_animal = cluster_dominant_animal;
metrics.cluster_dominance_percentage = cluster_dominance_percentage;
metrics.cluster_is_individual = cluster_is_individual;
metrics.animal_individual_clusters = animal_condition_individual_clusters;
metrics.animal_total_clusters_present = animal_condition_total_clusters_present;
metrics.animal_total_dominant_clusters = animal_condition_total_dominant_clusters;
metrics.animal_individual_frames = animal_condition_individual_frames;
metrics.animal_total_frames = animal_condition_total_frames;
metrics.animal_frame_individuality_pct = frame_indiv_pct_map;
metrics.total_global_frames = total_global_frames;
metrics.total_individual_frames = total_individual_frames;
metrics.overall_frame_individuality_pct = overall_frame_indiv_pct;
metrics.global_animal_condition_ids = global_animal_condition_ids; %#ok<STRNU>
metrics.global_animal_names = global_animal_names; %#ok<STRNU>

end

function plot_individuality_summary_across_zvals(exp_name, indiv_results, out_folder)
% Comprehensive summary plot across z-values variants per experiment.

% Determine experiment condition set and highlight pair
if contains(exp_name, 'BSFC')
    conds = {'B','S','F','C'};
    highlight = {'F','C'};
    title_suffix = 'BSFC (F vs C emphasis)';
elseif contains(exp_name, 'BHNG')
    conds = {'B','H','N','G'};
    highlight = {'N','G'};
    title_suffix = 'BHNG (N vs G emphasis)';
else
    conds = {'B','S','F','C'};
    highlight = {};
    title_suffix = exp_name;
end

nz = numel(indiv_results);
zcats = cell(1, nz);
cluster_share = nan(nz, numel(conds));
frame_indiv_pct = nan(nz, numel(conds));

for i = 1:nz
    zcats{i} = indiv_results(i).zfield;
    M = indiv_results(i).metrics;

    % Cluster-level: share of individual clusters dominated by each condition
    ci = M.cluster_is_individual(:)';
    dom = M.cluster_dominant_animal(:)';
    if any(ci)
        total_indiv = sum(ci);
        dom_cond = cellfun(@(s) s(end), dom(ci), 'UniformOutput', false);
        for c = 1:numel(conds)
            cluster_share(i,c) = sum(strcmp(dom_cond, conds{c})) / total_indiv;
        end
    end

    % Frame-level: percent frames in individual behaviors per condition
    indiv_frames_map = M.animal_individual_frames;
    total_frames_map = M.animal_total_frames;
    keys_total = keys(total_frames_map);
    for c = 1:numel(conds)
        cond_char = conds{c};
        tot_frames_c = 0; ind_frames_c = 0;
        for k = 1:numel(keys_total)
            ac_id = keys_total{k};
            if ac_id(end) == cond_char
                tot_frames_c = tot_frames_c + total_frames_map(ac_id);
                if isKey(indiv_frames_map, ac_id)
                    ind_frames_c = ind_frames_c + indiv_frames_map(ac_id);
                end
            end
        end
        if tot_frames_c > 0
            frame_indiv_pct(i,c) = ind_frames_c / tot_frames_c;
        else
            frame_indiv_pct(i,c) = NaN;
        end
    end
end

% Colors per condition (requested mapping):
% Baseline B: gray; S/H: green; F/N: red; C/G: blue
cmap = containers.Map();
cmap('B') = [0.5, 0.5, 0.5];
cmap('S') = [0.4660, 0.6740, 0.1880];
cmap('F') = [0.8500, 0.3250, 0.0980];
cmap('C') = [0, 0.4470, 0.7410];
cmap('H') = [0.4660, 0.6740, 0.1880];
cmap('N') = [0.8500, 0.3250, 0.0980];
cmap('G') = [0, 0.4470, 0.7410];

fig = figure('Color','w','Position',[100 100 1100 900],'Visible','off');

subplot(3,1,1);
hold on;
for c = 1:numel(conds)
    lw = 2; if any(strcmp(highlight, conds{c})), lw = 3.5; end
    plot(1:nz, cluster_share(:,c)*100, '-o', 'LineWidth', lw, 'Color', cmap(conds{c}), 'DisplayName', conds{c});
end
hold off;
set(gca, 'XTick', 1:nz, 'XTickLabel', zcats, 'XTickLabelRotation', 30);
ylabel('% of individual clusters');
title(sprintf('Cluster-level individuality share across z-values - %s', title_suffix));
grid on; legend('Location','eastoutside'); ylim([0 100]);
set(gca, 'Color', 'w', 'XColor','k', 'YColor','k');

subplot(3,1,2);
hold on;
for c = 1:numel(conds)
    lw = 2; if any(strcmp(highlight, conds{c})), lw = 3.5; end
    plot(1:nz, frame_indiv_pct(:,c)*100, '-o', 'LineWidth', lw, 'Color', cmap(conds{c}), 'DisplayName', conds{c});
end
hold off;
set(gca, 'XTick', 1:nz, 'XTickLabel', zcats, 'XTickLabelRotation', 30);
ylabel('% frames in individual behaviors');
xlabel('z-values variant');
title(sprintf('Frame-level individuality across z-values - %s', title_suffix));
grid on; legend('Location','eastoutside'); ylim([0 30]);
set(gca, 'Color', 'w', 'XColor','k', 'YColor','k');

% Pairwise difference subplot for highlighted conditions
if numel(highlight) == 2
    idx1 = find(strcmp(conds, highlight{1}));
    idx2 = find(strcmp(conds, highlight{2}));
    if ~isempty(idx1) && ~isempty(idx2)
        subplot(3,1,3);
        hold on;
        d_cluster = (cluster_share(:,idx1) - cluster_share(:,idx2))*100; % percentage points
        d_frame = (frame_indiv_pct(:,idx1) - frame_indiv_pct(:,idx2))*100;
        plot(1:nz, d_cluster, '-s', 'LineWidth', 2.5, 'Color', [0.2 0.2 0.2], 'DisplayName', sprintf('Individual clusters: %s - %s', highlight{1}, highlight{2}));
        plot(1:nz, d_frame, '--d', 'LineWidth', 2.5, 'Color', [0.3 0.6 0.9], 'DisplayName', sprintf('Frames: %s - %s', highlight{1}, highlight{2}));
        yline(0, ':', 'Color', [0.5 0.5 0.5]);
        hold off;
        set(gca, 'XTick', 1:nz, 'XTickLabel', zcats, 'XTickLabelRotation', 30);
        ylabel('Δ (pp)');
        xlabel('z-values variant');
        title(sprintf('Pairwise difference (highlight): %s - %s', highlight{1}, highlight{2}));
        % Explain sign of Δ: positive means first condition > second
        text(0.01, 0.95, sprintf('Positive Δ means %s > %s', highlight{1}, highlight{2}), ...
            'Units','normalized', 'Color','k', 'FontSize', 9, 'VerticalAlignment','top');
        grid on; legend('Location','eastoutside');
        set(gca, 'Color', 'w', 'XColor','k', 'YColor','k');
    end
end

% Save
out_fn = fullfile(out_folder, sprintf('individuality_summary_across_zvals_%s.pdf', exp_name));
exportgraphics(fig, out_fn, 'ContentType','vector', 'BackgroundColor','white');

close(fig);

% Also generate a separate figure comparing frame individuality differences across ALL condition pairs
% Build pairwise differences matrix: rows=pairs, cols=z-values variants
pairs = nchoosek(1:numel(conds), 2);
num_pairs = size(pairs,1);
if num_pairs >= 1 && nz >= 1
    diff_matrix = zeros(num_pairs, nz);
    pair_labels = cell(num_pairs,1);
    for p = 1:num_pairs
        i1 = pairs(p,1); i2 = pairs(p,2);
        % Δ in percentage points (cond_i1 - cond_i2)
        diff_matrix(p,:) = (frame_indiv_pct(:,i1) - frame_indiv_pct(:,i2))' * 100;
        pair_labels{p} = sprintf('%s-%s', conds{i1}, conds{i2});
    end

    fig2 = figure('Color','w','Position',[120 120 1100 max(450, 120 + 30*num_pairs)], 'Visible','off');
    imagesc(diff_matrix);
    colormap(parula);
    maxabs = max(abs(diff_matrix(:)));
    if isfinite(maxabs) && maxabs > 0
        caxis([-maxabs, maxabs]);
    end
    cb = colorbar; grid off;
    set(cb, 'Color','k');
    set(gca, 'YTick', 1:num_pairs, 'YTickLabel', pair_labels);
    set(gca, 'XTick', 1:nz, 'XTickLabel', zcats, 'XTickLabelRotation', 30);
    ylabel('Condition pairs');
    xlabel('z-values variant');
    title(sprintf('Pairwise differences in frame individuality (pp) across z-values - %s', title_suffix));
    set(gca, 'Color', 'w', 'XColor','k', 'YColor','k');

    out_fn2 = fullfile(out_folder, sprintf('individuality_frame_diffs_allpairs_%s.pdf', exp_name));
    exportgraphics(fig2, out_fn2, 'ContentType','vector', 'BackgroundColor','white');
    close(fig2);
end
end
