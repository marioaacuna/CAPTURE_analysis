% 1. Initialization
% clear;
close all;
clc;
GC = general_configs;
rootpath = GC.preprocessing_rootpath;

% Load Data
load(GC.filename_analysis, 'analysisstruct')
load(GC.filename_ratception, 'ratception_struct');
load(GC.filename_predictions, 'animal_condition_identifier');
input_params.repfactor = GC.repfactor;

% Preprocess data
markers_aligned_ds = load_aligned_markers(ratception_struct.markers_aligned_preproc, input_params.repfactor, 10);
markers_not_aligned_ds = load_aligned_markers(ratception_struct.markers_preproc, input_params.repfactor, 15);
% Extract conditions
frame_identifiers = animal_condition_identifier;
conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);
unique_conditions = unique(conditions);

%% 2.1 calculate walking bouts based on 2D positions (x,y)
% Basic usage
% walking_bouts = detectWalkingFrom2D(markers_not_aligned_ds.SpineM);

% With custom parameters
params = struct();
params.velocity_percentile = 87; % More stringent threshold
params.min_bout_duration = 0.1; % Longer minimum bout
params.sampling_rate = 100; % Hz
params.smoothing_window = 5; % frames
params.z_smoothing_window = 5;
params.z_threshold_percentile = 99.3;  % Threshold for Z displacement
params.direction_threshold = 75;  % Max angle deviation from heading (degrees)
params.cycle_window = [-0.5 0.5]; % Time window for cycle normalization (seconds)
params.n_normalized_points = 200; % Number of points after interpolation
params.interp_method = 'pchip'; % Interpolation method (pchip/linear/spline)
params.markers_to_study = {'WristL', 'WristR', 'KneeL', 'KneeR','AnkleL', 'AnkleR', 'HindpawL', 'HindpawR'}; % Markers to study

free_SpineM = markers_not_aligned_ds.SpineM;
free_Snout = markers_not_aligned_ds.Snout;
[walking_bouts, metrics] = detectWalkingFrom2D(free_SpineM,free_Snout, params);

%% 3. Analysis of angles at walking in egocentric reference
markers = markers_aligned_ds;

% Calculate angles and perform analysis
analyze_angles(markers, walking_bouts, conditions, unique_conditions, frame_identifiers, animal_condition_identifier);

disp('done')

% Output a small video with a walking example

s = struct();
s.markers_aligned_preproc = markers_aligned_ds;
s.markernames = ratception_struct.markernames;
s.markercolor =ratception_struct.markercolor; 
s.markers_preproc = markers_not_aligned_ds;
s.links = ratception_struct.links;
h = figure;
bout_starts = find(diff([0; walking_bouts]) == 1);
bout_ends = find(diff([walking_bouts; 0]) == -1);

frames_to_plot = [bout_starts(1):bout_ends(1),bout_starts(3): bout_ends(3)];
g = animate_markers_aligned_fullmovie_demo(s,frames_to_plot, h , '');

% Create GIF
gif_filename = 'walking_analysis.gif';
h = figure;
for frame = 1:length(frames_to_plot)
    % Plot 
    animate_markers_aligned_fullmovie_demo(s,frames_to_plot(frame), h , '');

    
    % Capture the plot as an image
    frame_img = getframe(h);
    img = frame2im(frame_img);
    [imind, cm] = rgb2ind(img, 256);
    
    % Write to the GIF File
    if frame == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

% % Create MP4
% mp4_filename = 'walking_analysis.mp4';
% v = VideoWriter(mp4_filename, 'MPEG-4');
% v.FrameRate = 10; % Adjust frame rate as needed
% open(v);
% for frame = 1:length(frames_to_plot)
%     % Plot 
       % animate_markers_aligned_fullmovie_demo(s,frames_to_plot(frame), h , '');
% 
%     % Capture the plot as an image
%     frame_img = getframe(h);
% 
%     % Write to the video file
%     writeVideo(v, frame_img);
% end
% close(v);

%% 4. Kinematic analysis of gait cycles
cycles_all = analyze_gait_cycles_3d(markers, walking_bouts, params);

%% 4.1 Kinematics of gait cycles per condition
CYCLES = struct();
% Separate data by conditions
for c = 1:length(unique_conditions)
    condition = unique_conditions{c};
    condition_mask = strcmp(conditions, condition) ;
    
    % Perform gait analysis for each condition
    fprintf('Gait analysis for condition: %s\n', condition);
    
    % Restrict markers to current condition
    condition_markers = structfun(@(x) x(condition_mask, :), markers, 'UniformOutput', false);
    
    % Analyze gait cycles for the current condition
    cycles = analyze_gait_cycles_3d(condition_markers, walking_bouts(strcmp(conditions, condition)), params);
    
    CYCLES.(condition) = cycles;
end

disp('done')

%% 5. Plotting cycles per marker for each condition
marker_names = params.markers_to_study;
n_markers = length(marker_names);
n_conditions = length(unique_conditions);

% Determine subplot grid size
n_cols = ceil(sqrt(n_markers));
n_rows = ceil(n_markers / n_cols);

% Define colors for each condition
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], ...
          [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
legendHandles = [];  % Store handles for legend

figure;
for m = 1:n_markers
    subplot(n_rows, n_cols, m);
    hold on;
    for c = 1:n_conditions
        condition = unique_conditions{c};
        cycles = CYCLES.(condition);
        
        % Extract data for the current marker
        marker_data = arrayfun(@(x) x.(marker_names{m}), cycles, 'UniformOutput', false);
        marker_data = cat(1, marker_data{:});
        
        % Calculate mean and SEM
        marker_mean = mean(marker_data, 1);
        marker_sem = std(marker_data, 0, 1) / sqrt(size(marker_data, 1));
        
        % Plot mean and SEM as shaded area
        t = linspace(params.cycle_window(1), params.cycle_window(2), params.n_normalized_points);
        fill([t, fliplr(t)], [marker_mean + marker_sem, fliplr(marker_mean - marker_sem)], ...
            colors{c}, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        h_line = plot(t, marker_mean, 'Color', colors{c}, 'LineWidth', 1.5);
        
        if m == 1  % Only store handles from first subplot for legend
            legendHandles = [legendHandles h_line];
        end
    end
    title(marker_names{m});
    xlabel('Time (s)');
    ylabel('Position (AU)');
    grid on;
end
legend(legendHandles, unique_conditions, 'Location', 'best');
sgtitle('Mean Gait Cycle per Marker and Condition');
hold off;

%% 6. Statistical Analysis and Plotting for Each Marker
toggle_toolbox('spm1d', 'on')
for m = 1:n_markers
    marker_name = marker_names{m};
    
    % Basic time-point analysis
    [p_values, significant_timepoints] = analyze_timepoints(CYCLES, marker_name, unique_conditions, params);
    
    % Extract features across all conditions
    features = extract_cycle_features(CYCLES, marker_name, unique_conditions);
    
    % Create visualization of statistical results
    plot_statistical_results(CYCLES, marker_name, unique_conditions, params);
    
    % Add super title for the marker
    sgtitle(['Statistical Analysis for Marker: ', marker_name]);
    
    % Compute pattern similarity
    similarity_results = analyze_pattern_similarity(CYCLES, marker_name, unique_conditions);
    
    % Visualize results
    plot_similarity_results(similarity_results, unique_conditions, marker_name);
    
    % Print similarity results
    fprintf('\nDTW Distances between conditions for %s:\n', marker_name);
    for i = 1:length(similarity_results.condition_pairs)
        pair = similarity_results.condition_pairs{i};
        fprintf('%s vs %s: %.3f\n', pair{1}, pair{2}, pair{3});
    end
    
    % Run SPM analysis
    spm_results = analyze_waveform_spm(CYCLES, marker_name, unique_conditions);
    
    % Visualize results
    plot_spm_results(spm_results, unique_conditions, params, marker_name);
    
    % Plot pairwise comparisons
    plot_spm_pairwise_results(spm_results, unique_conditions, params, marker_name);
    
    % Print summary of significant clusters
    fprintf('\nSPM Analysis Results for %s:\n', marker_name);
    fprintf('Number of significant clusters: %d\n', length(spm_results.inference.clusters));
    
    for i = 1:length(spm_results.inference.clusters)
        cluster = spm_results.inference.clusters{i};
        fprintf('\nCluster %d:\n', i);
        fprintf('P-value: %.4f\n', cluster.P);
        fprintf('Extent: %.2f%% to %.2f%% of gait cycle\n', ...
            cluster.endpoints(1)/params.n_normalized_points * 100, ...
            cluster.endpoints(2)/params.n_normalized_points * 100);
    end
    
    % Print summary of pairwise comparisons
    fprintf('\nPairwise Comparison Results for %s:\n', marker_name);
    pair_names = fieldnames(spm_results.pairwise);
    for p = 1:length(pair_names)
        pair = spm_results.pairwise.(pair_names{p});
        fprintf('\n%s:\n', strrep(pair_names{p}, '_', ' '));
        if pair.inference.nClusters > 0
            for i = 1:length(pair.inference.clusters)
                fprintf('Significant difference (p = %.4f) from %.2f to %.2f\n', ...
                    pair.inference.clusters{i}.P, ...
                    pair.inference.clusters{i}.endpoints(1), ...
                    pair.inference.clusters{i}.endpoints(2));
            end
        else
            fprintf('No significant differences\n');
        end
    end
end
toggle_toolbox('spm1d', 'off')

%% --- END OF SCRIPT ---




% 1. Time-Point Analysis
function [p_values, significant_timepoints] = analyze_timepoints(CYCLES, marker_name, unique_conditions, params)
    % Initialize arrays for timepoint analysis
    n_timepoints = params.n_normalized_points;
    p_values = zeros(1, n_timepoints);
    significant_timepoints = false(1, n_timepoints);
    
    % For each time point
    for t = 1:n_timepoints
        % Extract data at this timepoint for each condition
        condition_data = cell(length(unique_conditions), 1);
        group_labels = [];
        
        for c = 1:length(unique_conditions)
            condition = unique_conditions{c};
            timepoint_values = [];
            
            % Loop through all cycle groups in this condition
            for cycle_group = 1:length(CYCLES.(condition))
                cycle_data = CYCLES.(condition)(cycle_group).(marker_name);
                timepoint_values = [timepoint_values; cycle_data(:,t)];
            end
            
            condition_data{c} = timepoint_values;
            group_labels = [group_labels; repmat(c, length(timepoint_values), 1)];
        end
        
        % Combine all data for this timepoint
        all_data = vertcat(condition_data{:});
        
        % Perform ANOVA
        [p_values(t), ~, ~] = anova1(all_data, group_labels, 'off');
        significant_timepoints(t) = p_values(t) < 0.05;
    end
end

% 2. Feature-Based Analysis
function features = extract_cycle_features(CYCLES, marker_name, unique_conditions)
    features = struct();
    
    for c = 1:length(unique_conditions)
        condition = unique_conditions{c};
        peak_vals = [];
        rom_vals = [];
        
        % Loop through all cycle groups
        for cycle_group = 1:length(CYCLES.(condition))
            cycles = CYCLES.(condition)(cycle_group).(marker_name);
            
            % For each cycle in this group
            for cycle = 1:size(cycles,1)
                cycle_data = cycles(cycle,:);
                
                % Calculate features
                peak_vals = [peak_vals; max(cycle_data)];
                rom_vals = [rom_vals; range(cycle_data)];
            end
        end
        
        % Store features for this condition
        features.(condition).peaks = peak_vals;
        features.(condition).rom = rom_vals;
        features.(condition).mean_peak = mean(peak_vals);
        features.(condition).std_peak = std(peak_vals);
        features.(condition).mean_rom = mean(rom_vals);
        features.(condition).std_rom = std(rom_vals);
    end
end

% Example usage:
% marker_name = 'WristL';
% [p_vals, sig_points] = analyze_timepoints(CYCLES, marker_name, unique_conditions, params);
% features = extract_cycle_features(CYCLES, marker_name, unique_conditions);

% Optional: Visualization of statistical results
function plot_statistical_results(CYCLES, marker_name, unique_conditions, params)
    [p_values, sig_timepoints] = analyze_timepoints(CYCLES, marker_name, unique_conditions, params);
    
    figure;
    
    % Plot p-values
    subplot(2,1,1);
    plot(linspace(params.cycle_window(1), params.cycle_window(2), params.n_normalized_points), ...
         -log10(p_values), 'b-', 'LineWidth', 1.5);
    hold on;
    yline(-log10(0.05), 'r--', 'Significance threshold');
    xlabel('Time (s)');
    ylabel('-log10(p-value)');
    title(['Statistical Significance Across Gait Cycle', ' ', marker_name]);
    grid on;
    
    % Plot feature comparisons
    features = extract_cycle_features(CYCLES, marker_name, unique_conditions);
    
    subplot(2,1,2);
    peak_means = [];
    peak_sems = [];
    for c = 1:length(unique_conditions)
        condition = unique_conditions{c};
        peak_means(c) = features.(condition).mean_peak;
        peak_sems(c) = features.(condition).std_peak / sqrt(length(features.(condition).peaks));
    end
    
    bar(peak_means);
    hold on;
    errorbar(1:length(unique_conditions), peak_means, peak_sems, 'k.', 'LineWidth', 1.5);
    xlabel('Condition');
    ylabel('Peak Value');
    title(['Peak Comparison Across Conditions', ' ', marker_name]);
    set(gca, 'XTickLabel', unique_conditions);
    
    sgtitle(['Statistical Analysis for Marker: ', marker_name]);
end




function similarity_results = analyze_pattern_similarity(CYCLES, marker_name, unique_conditions)
    % Initialize storage for results
    similarity_results = struct();
    similarity_results.dtw_distances = zeros(length(unique_conditions));
    similarity_results.condition_pairs = {};
    similarity_results.mean_cycles = struct();
    
    % First compute mean cycle for each condition
    for c = 1:length(unique_conditions)
        condition = unique_conditions{c};
        all_cycles = [];
        
        % Gather all cycles from this condition
        for cycle_group = 1:length(CYCLES.(condition))
            cycles = CYCLES.(condition)(cycle_group).(marker_name);
            all_cycles = [all_cycles; cycles];
        end
        
        % Store mean cycle
        similarity_results.mean_cycles.(condition) = mean(all_cycles, 1);
    end
    
    % Compute DTW distances between all condition pairs
    pair_idx = 1;
    for c1 = 1:length(unique_conditions)
        for c2 = c1+1:length(unique_conditions)
            cond1 = unique_conditions{c1};
            cond2 = unique_conditions{c2};
            
            % Get mean cycles
            cycle1 = similarity_results.mean_cycles.(cond1);
            cycle2 = similarity_results.mean_cycles.(cond2);
            
            % Calculate DTW distance
            dtw_dist = dtw(cycle1', cycle2');
            
            % Store results
            similarity_results.dtw_distances(c1,c2) = dtw_dist;
            similarity_results.dtw_distances(c2,c1) = dtw_dist;
            similarity_results.condition_pairs{pair_idx} = {cond1, cond2, dtw_dist};
            pair_idx = pair_idx + 1;
        end
    end
end

% Function to visualize the similarity results
function plot_similarity_results(similarity_results, unique_conditions,marker_name)
    figure('Position', [100 100 800 800]);
    
    % Plot distance matrix
    subplot(2,1,1);
    imagesc(similarity_results.dtw_distances);
    colorbar;
    clim([0 80])
    title('DTW Distances Between Conditions');
    set(gca, 'XTick', 1:length(unique_conditions), 'XTickLabel', unique_conditions);
    set(gca, 'YTick', 1:length(unique_conditions), 'YTickLabel', unique_conditions);
    colormap('parula');
    
    % Plot mean cycles overlay
    subplot(2,1,2);
    hold on;
    colors = lines(length(unique_conditions));
    
    for c = 1:length(unique_conditions)
        condition = unique_conditions{c};
        plot(similarity_results.mean_cycles.(condition), 'Color', colors(c,:), ...
            'LineWidth', 2, 'DisplayName', condition);
    end
    
    title('Mean Cycles Comparison');
    xlabel('Normalized Time Points');
    ylabel('Position');
    legend('Location', 'best');
    grid on;
    hold off;
    
    sgtitle(['Pattern Similarity Results for marker ', marker_name]);
end


function spm_results = analyze_waveform_spm(CYCLES, marker_name, unique_conditions)
    % Initialize storage
    spm_results = struct();
    
    % Prepare data arrays for SPM
    Y = [];  % All cycles data
    A = [];  % Condition labels
    group_count = 1;
    
    % Gather data from all conditions
    for c = 1:length(unique_conditions)
        condition = unique_conditions{c};
        
        % Get all cycles from this condition
        for cycle_group = 1:length(CYCLES.(condition))
            cycles = CYCLES.(condition)(cycle_group).(marker_name);
            
            % Add cycles to data array
            Y = [Y; cycles];
            A = [A; repmat(c, size(cycles,1), 1)];
            
            % Store cycle count for this condition
            if c == 1
                group_count = group_count + size(cycles,1);
            end
        end
    end
    
    % Perform SPM analysis
    spm_results.anova = spm1d.stats.anova1(Y, A);
    spm_results.inference = spm_results.anova.inference(0.05);
      
    % Add pairwise comparisons
    spm_results.pairwise = struct();
    
    % Perform all pairwise comparisons
    for i = 1:length(unique_conditions)-1
        for j = i+1:length(unique_conditions)
            cond1 = unique_conditions{i};
            cond2 = unique_conditions{j};
            
            % Get data for these two conditions
            idx = (A == i) | (A == j);
            Y_pair = Y(idx, :);
            A_pair = A(idx);
            
            % Perform t-test
            t_test = spm1d.stats.ttest2(Y_pair(A_pair == i, :), ...
                                      Y_pair(A_pair == j, :));
            inference = t_test.inference(0.05);
            
            % Store results
            pair_name = sprintf('%s_vs_%s', cond1, cond2);
            spm_results.pairwise.(pair_name).ttest = t_test;
            spm_results.pairwise.(pair_name).inference = inference;
        end
    end
    
    % Store additional information
    spm_results.Y = Y;
    spm_results.A = A;
    spm_results.group_count = group_count;
end

function plot_spm_results(spm_results, unique_conditions, params, marker_name)
    figure('Position', [100 100 1200 800]);
    
    % Plot 1: Mean cycles by condition
    subplot(2,1,1);
    hold on;
    colors = lines(length(unique_conditions));
    
    for c = 1:length(unique_conditions)
        condition_cycles = spm_results.Y(spm_results.A == c, :);
        mean_cycle = mean(condition_cycles, 1);
        std_cycle = std(condition_cycles, 0, 1) / sqrt(size(condition_cycles,1));
        
        t = linspace(params.cycle_window(1), params.cycle_window(2), size(condition_cycles,2));
        
        % Plot mean and standard deviation
        fill([t fliplr(t)], [mean_cycle+std_cycle fliplr(mean_cycle-std_cycle)], ...
            colors(c,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(t, mean_cycle, 'Color', colors(c,:), 'LineWidth', 2, ...
            'DisplayName', unique_conditions{c});
    end
    
    title('Mean Cycles by Condition');
    xlabel('Time');
    ylabel('Position');
    legend('Location', 'best');
    grid on;
    
    % Plot 2: SPM results
    subplot(2,1,2);
    
    % Plot SPM{F} statistic
    plot(t, spm_results.anova.z, 'k-', 'LineWidth', 1.5);
    hold on;
    
    % Add threshold line
    hline = refline([0 spm_results.inference.zstar]);
    hline.Color = 'r';
    hline.LineStyle = '--';
    
    % Normalize cluster time points to cycle window
    n_points = length(spm_results.anova.z);
    time_to_window = @(x) params.cycle_window(1) + (x/n_points) * (params.cycle_window(2) - params.cycle_window(1));
    
    % Highlight significant clusters if they exist
    if spm_results.inference.nClusters > 0
        clusters = spm_results.inference.clusters;
        for i = 1:length(clusters)
            % Convert cluster endpoints to normalized time
            start_point = time_to_window(clusters{i}.endpoints(1));
            end_point = time_to_window(clusters{i}.endpoints(2));
            
            % Create normalized time vector for this cluster
            t_cluster = linspace(start_point, end_point, 100);
            
            % Get corresponding F-values through interpolation
            z_cluster = interp1(t, spm_results.anova.z, t_cluster);
            
            % Highlight significant region
            fill([t_cluster fliplr(t_cluster)], ...
                [z_cluster fliplr(zeros(size(z_cluster)))], ...
                'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            
            % Add text annotation for cluster p-value
            text(mean([start_point, end_point]), max(z_cluster), ...
                sprintf('p = %.4f', clusters{i}.P), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom');
        end
    end
    
    title('SPM Analysis Results');
    xlabel('Time');
    ylabel('F-statistic');
    grid on
   
    sgtitle(['SPM Analysis Results for ', marker_name]);
end

function plot_spm_pairwise_results(spm_results, unique_conditions, params, marker_name)
    % Get all pairwise comparison names
    pair_names = fieldnames(spm_results.pairwise);
    n_pairs = length(pair_names);
    
    % Create figure with subplots for each comparison
    figure('Position', [100 100 1200 800]);
    
    for p = 1:n_pairs
        subplot(ceil(n_pairs/5),ceil(n_pairs/2), p);
        pair = spm_results.pairwise.(pair_names{p});
        
        % Time vector
        t = linspace(params.cycle_window(1), params.cycle_window(2), size(spm_results.Y, 2));
        
        % Plot t-statistic
        plot(t, pair.ttest.z, 'k-', 'LineWidth', 1.5);
        hold on;
        
        % Add threshold lines
        yline(pair.inference.zstar, 'r--');
        yline(-pair.inference.zstar, 'r--');
        
        % Highlight significant clusters
        if pair.inference.nClusters > 0
            for i = 1:length(pair.inference.clusters)
                cluster = pair.inference.clusters{i};
                % Convert cluster endpoints to normalized time
                start_point = params.cycle_window(1) + ...
                    (cluster.endpoints(1)/size(spm_results.Y,2)) * ...
                    (params.cycle_window(2) - params.cycle_window(1));
                end_point = params.cycle_window(1) + ...
                    (cluster.endpoints(2)/size(spm_results.Y,2)) * ...
                    (params.cycle_window(2) - params.cycle_window(1));
                
                t_cluster = linspace(start_point, end_point, 100);
                z_cluster = interp1(t, pair.ttest.z, t_cluster);
                
                fill([t_cluster fliplr(t_cluster)], ...
                    [z_cluster fliplr(zeros(size(z_cluster)))], ...
                    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                
                % Add p-value
                text(mean([start_point, end_point]), max(z_cluster), ...
                    sprintf('p = %.4f', cluster.P), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom');
            end
        end
        
        title(strrep(pair_names{p}, '_', ' '));
        xlabel('Time');
        ylabel('t-statistic');
        grid on;
    end
    
    sgtitle(['Pairwise Comparison Results for ', marker_name]);
end