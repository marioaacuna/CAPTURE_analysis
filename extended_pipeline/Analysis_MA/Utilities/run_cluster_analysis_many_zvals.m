%% Run clustering and plotting for multiple zValues fields per experiment
% Loops through configured experiments and selected zValues variants in
% analysisstruct, performs clustering, plots, and saves outputs and structs.
% Data comes from compute_zvals_diff_features.m 

clear; clc;
global GC
GC = general_configs();

% Experiments and folders (same as in plot_tsnemaps_per_condition_drugs_many_zvals)
experiments = struct();
experiments(1).name = 'BSFC_300hz';
experiments(1).folder = '0_preprocessing_BSFC_300hz';
experiments(2).name = 'BHNG_300hz';
experiments(2).folder = '0_preprocessing_BHNG_300hz';

% zValues variants to process
zvals_fields =  {'zValues_all_features', ...
                 'zValues_important_features', ...
                 'zValues_jt_features_only', ...
                 'zValues_hand_made', ...
                 'zValues_only_important'};

% Plot options
plot_poses = true;

for exp_idx = 1:numel(experiments)
    exp = experiments(exp_idx);
    fprintf('--- Experiment: %s ---\n', exp.name);

    data_folder = fullfile(GC.project_path, 'data', exp.folder);
    if ~exist(data_folder, 'dir')
        warning('Data folder not found: %s', data_folder);
        continue;
    end

    fn_analysis = fullfile(data_folder, 'raw_concat_analysis.mat');
    if ~exist(fn_analysis, 'file')
        warning('Analysis file not found: %s', fn_analysis);
        continue;
    end

    S = load(fn_analysis, 'analysisstruct');
    if ~isfield(S, 'analysisstruct')
        warning('analysisstruct var not found in %s', fn_analysis);
        continue;
    end
    %analysisstruct = S.analysisstruct; %#ok<NASGU>

    % Export and save folders per experiment
    exp_out_folder = fullfile(GC.temp_root, 'clusters_many_zvals', exp.name);
    if ~exist(exp_out_folder, 'dir'); mkdir(exp_out_folder); end

    for zi = 1:numel(zvals_fields)
        zf = zvals_fields{zi};
        analysisstruct = S.analysisstruct; 

        if ~isfield(analysisstruct, zf)
            fprintf('Skipping %s (field not found)\n', zf);
            continue;
        end
        fprintf('Processing zValues field: %s\n', zf);

        % Cluster on this zValues
        analysisstruct_temp = cluster_analysis_for_zval(analysisstruct, zf, GC);

        % Plot and export
        zf_folder = fullfile(exp_out_folder, zf);
        if ~exist(zf_folder, 'dir'); mkdir(zf_folder); end
        try
            plot_and_export_clusters(analysisstruct_temp, zf_folder, exp.name, plot_poses);
        catch ME
            warning('Plotting failed for %s: %s', zf, ME.message);
        end

        % Save updated analysisstruct (per zval) alongside original data
        analysis_filename = fullfile(data_folder, sprintf('analysis_%s.mat', zf));
        try
            disp('saving analysis struct ...');
            analysisstruct_out = analysisstruct_temp; 
            analysisstruct = analysisstruct_temp; % % back to analysisstruct
            save(analysis_filename, 'analysisstruct', '-v7.3');
        catch ME
            warning('Failed saving analysisstruct for %s: %s', zf, ME.message);
        end
    end

    fprintf('Done experiment: %s\n', exp.name);
end

fprintf('All experiments complete.\n');
