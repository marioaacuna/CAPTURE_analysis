function analyze_gait_cycles_3d(markers_aligned, walking_bouts, params)
% The 3 dimensions in aligned_cycles(1).wrist_l are:
% 
% Heading direction (forward/backward movement)
% Mediolateral direction (left/right movement)
% Vertical direction (up/down movement)   
% 

% % Initialize markers and parameters
   params = initialize_params(params);
   markers = initialize_markers(markers_aligned);
   
   % Process walking bouts
   all_cycles = process_walking_bouts(markers, walking_bouts, params);
   
   % Align and analyze cycles
   aligned_cycles = align_cycles_3d(all_cycles, params);
   
   % Visualize results
   if params.plot_results
       plot_3d_cycles(aligned_cycles, params);
   end
end


function params = initialize_params(params)
   default_params = struct(...
       'sampling_rate', 100,...       % Data acquisition frequency (Hz)
       'min_segment_duration', 1.2,...    % Minimum duration of valid walking bout (seconds)
       'min_cycle_length', 0.2,...     % Minimum duration of gait cycle (seconds)
       'max_cycle_length', 1,...     % Maximum duration of gait cycle (seconds)
       'window_size', 5,...            % Smoothing window size (samples)
       'cycle_window', [-0.2 0.2],...  % Time window for cycle normalization (seconds)
       'n_normalized_points', 200,... % Number of points after interpolation
       'interp_method', 'pchip',...    % Interpolation method (pchip/linear/spline)
       'cycle_detection_method', 'peak_detection',... % Method for detecting cycles
       'alignment_method', 'fixed_window',... % Time normalization approach
       'plot_results', true,...        % Enable/disable plotting
       'plot_dimensions', {{'X (Forward)', 'Y (Lateral)', 'Z (Vertical)'}},... % Axis labels
       'marker_colors', {{[0.5 0 0.5], [1 0 0], [0 1 1], [1 1 0]}}); % Colors for each marker

   params = merge_structs(default_params, params);
end

function markers = initialize_markers(markers_aligned)
    markers = struct('wrist_l', markers_aligned.WristL, ...
                    'wrist_r', markers_aligned.WristR, ...
                    'ankle_l', markers_aligned.AnkleL, ...
                    'ankle_r', markers_aligned.AnkleR, ...
                    'spine_m', markers_aligned.SpineM);

end

function all_cycles = process_walking_bouts(markers, walking_bouts, params)
   bout_starts = find(diff([0; walking_bouts]) == 1);
   bout_ends = find(diff([walking_bouts; 0]) == -1);
   
   all_cycles = struct();
   cycle_count = 0;
   
   for bout_idx = 1:length(bout_starts)
       segment_idx = bout_starts(bout_idx):bout_ends(bout_idx);
       
       if ~is_valid_segment(segment_idx, params)
           continue
       end
       
       norm_markers = normalize_markers_3d(markers, segment_idx);
       [cycles, valid] = detect_cycles(norm_markers, params);
       
       all_cycles = store_valid_cycles(cycles, valid, all_cycles, cycle_count);
       cycle_count = cycle_count + length(valid);
   end
end

function [cycles, valid] = detect_cycles(norm_markers, params)
   marker_names = fieldnames(norm_markers);
   cycles = struct();
   
   % Get cycle indices based on combined marker motions
   switch params.cycle_detection_method
       case 'zero_crossing'
           [starts, ends] = detect_by_zero_crossing(norm_markers, marker_names, params);
       case 'peak_detection'
           [starts, ends] = detect_by_peaks(norm_markers, marker_names, params);
       case 'phase_space'
           [starts, ends] = detect_by_phase(norm_markers, marker_names, params);
   end
   
   [cycles, valid] = extract_valid_cycles(norm_markers, starts, ends, params);
end

function [starts, ends] = detect_by_zero_crossing(markers, marker_names, params)
   combined_signals = [];
   
   % Combine signals from all markers
   for i = 1:length(marker_names)
       if ~strcmp(marker_names{i}, 'spine_m')
           marker = markers.(marker_names{i});
           sagittal = compute_sagittal_position(marker);
           combined_signals = [combined_signals sagittal];
       end
   end
   
   % Use mean of all markers for robust cycle detection
   mean_signal = mean(combined_signals, 2);
   zero_crosses = find(diff(sign(mean_signal)) > 0);
   
   starts = zero_crosses(1:end-1);
   ends = zero_crosses(2:end);
end
function [starts, ends] = detect_by_peaks(markers, marker_names, params)
   combined_signals = [];
   
   for i = 1:length(marker_names)
       if ~strcmp(marker_names{i}, 'spine_m')
           marker = markers.(marker_names{i});
           sagittal = compute_sagittal_position(marker);
           combined_signals = [combined_signals sagittal];
       end
   end
   
   mean_signal = mean(combined_signals, 2);
   [~, peak_locs] = findpeaks(mean_signal, 'MinPeakDistance', round(params.min_cycle_length * params.sampling_rate));
   
   starts = peak_locs(1:end-1);
   ends = peak_locs(2:end);
end

function [starts, ends] = detect_by_phase(markers, marker_names, params)
   combined_signals = [];
   velocities = [];
   
   for i = 1:length(marker_names)
       if ~strcmp(marker_names{i}, 'spine_m')
           marker = markers.(marker_names{i});
           sagittal = compute_sagittal_position(marker);
           velocity = diff(sagittal) * params.sampling_rate;
           velocity = [velocity; velocity(end)];
           
           combined_signals = [combined_signals sagittal];
           velocities = [velocities velocity];
       end
   end
   
   mean_pos = mean(combined_signals, 2);
   mean_vel = mean(velocities, 2);
   
   % Find zero crossings in phase space
   phase_angle = atan2(mean_vel, mean_pos);
   crossings = find(diff(phase_angle) < -pi); % Full cycle completion
   
   starts = crossings(1:end-1);
   ends = crossings(2:end);
end


function aligned = align_cycles_3d(cycles, params)
   t_norm = create_normalized_time(params);
   aligned = initialize_aligned_structure(cycles, params);
   
   for i = 1:length(cycles)
       t_orig = create_cycle_time(cycles(i), params);
       aligned = interpolate_cycle(cycles(i), aligned, i, t_orig, t_norm, params);
   end
end

function norm_markers = normalize_markers_3d(markers, segment_idx)
    marker_names = fieldnames(markers);
    norm_markers = struct();
    
    for i = 1:length(marker_names)
        if ~strcmp(marker_names{i}, 'spine_m')
            % Just extract the segment data as it's already normalized
            norm_markers.(marker_names{i}) = markers.(marker_names{i})(segment_idx,:);
        end
    end
end

% function marker_norm = normalize_trajectories_3d(marker, spine)
%    marker_cent = marker - spine;
%    spine_diff = compute_spine_direction(spine);
% 
%    marker_norm = zeros(size(marker));
%    for i = 1:size(spine_diff,1)
%        R = compute_rotation_matrix(spine_diff(i,:));
%        marker_norm(i,:) = (R * marker_cent(i,:)')';
%    end
% end

function spine_diff = compute_spine_direction(spine)
   spine_diff = diff(spine);
   spine_diff = [spine_diff; spine_diff(end, :)]; % Match dimensions
end

function R = compute_rotation_matrix(direction)
   x_dir = normalize_vector(direction);
   approx_up = [0 0 1];
   y_dir = normalize_vector(cross(approx_up, x_dir));
   z_dir = cross(x_dir, y_dir);
   R = [x_dir; y_dir; z_dir]';
end

function [cycles, valid] = extract_valid_cycles(norm_markers, starts, ends, params)
   cycles = struct();
   valid = false(length(starts), 1);
   marker_names = fieldnames(norm_markers);
   
   for i = 1:length(starts)
       cycle_length = ends(i) - starts(i);
       
       if is_valid_cycle_length(cycle_length, params)
           valid(i) = true;
           for m = 1:length(marker_names)
               if ~strcmp(marker_names{m}, 'spine_m')
                   cycles(i).(marker_names{m}) = ...
                       norm_markers.(marker_names{m})(starts(i):ends(i), :);
               end
           end
       end
   end
end

function valid = is_valid_cycle_length(cycle_length, params)
   min_samples = round(params.min_cycle_length * params.sampling_rate);
   max_samples = round(params.max_cycle_length * params.sampling_rate);
   valid = cycle_length >= min_samples && cycle_length <= max_samples;
end

function all_cycles = store_valid_cycles(cycles, valid, all_cycles, cycle_count)
   % Get valid cycles only
   valid_cycles = cycles(valid);
   
   % Initialize all_cycles if empty
   if isempty(fieldnames(all_cycles))
       all_cycles = valid_cycles;
       return;
   end
   
   % Ensure structures have same fields
   if ~isempty(valid_cycles)
       fields1 = fieldnames(all_cycles);
       fields2 = fieldnames(valid_cycles);
       if ~isequal(fields1, fields2)
           valid_cycles = harmonize_fields(all_cycles, valid_cycles);
       end
       
       % Append valid cycles
       for i = 1:length(valid_cycles)
           all_cycles(cycle_count + i) = valid_cycles(i);
       end
   end
end

function harmonized = harmonize_fields(template, target)
   fields = fieldnames(template);
   harmonized = target;
   
   for i = 1:length(fields)
       if ~isfield(harmonized, fields{i})
           [harmonized.(fields{i})] = deal([]);
       end
   end
end

function t_norm = create_normalized_time(params)
   t_norm = linspace(params.cycle_window(1), params.cycle_window(2), params.n_normalized_points);
end

function t_orig = create_cycle_time(cycle, params)
   marker_names = fieldnames(cycle);
   first_marker = marker_names{1};
   t_orig = linspace(0, size(cycle.(first_marker), 1) / params.sampling_rate, size(cycle.(first_marker), 1));
end
function aligned = initialize_aligned_structure(cycles, params)
   marker_names = fieldnames(cycles(1));
   aligned = struct();
   
   for i = 1:length(cycles)
       for m = 1:length(marker_names)
           if ~strcmp(marker_names{m}, 'duration')
               aligned(i).(marker_names{m}) = zeros(params.n_normalized_points, 3);
           end
       end
   end
end

function aligned = interpolate_cycle(cycle, aligned, idx, t_orig, t_norm, params)
    marker_names = fieldnames(cycle);
    
    for m = 1:length(marker_names)
        if ~strcmp(marker_names{m}, 'duration')
            data = cycle.(marker_names{m});
            if size(data, 1) >= 2  % Check minimum points requirement
                aligned(idx).(marker_names{m}) = interp1(t_orig, data, t_norm, params.interp_method);
            else
                warning('Cycle %d, marker %s has insufficient points for interpolation', idx, marker_names{m});
                aligned(idx).(marker_names{m}) = zeros(length(t_norm), size(data, 2));
            end
        end
    end
end


function plot_3d_cycles(aligned_cycles, params)
   marker_names = fieldnames(aligned_cycles(1));
   marker_names = marker_names(~strcmp(marker_names, 'duration'));
   
   % Initialize arrays for means and SEMs
   means = struct();
   sems = struct();

   % For each marker, we need to combine all cycles
   n_cycles = length(aligned_cycles);
   n_timepoints = size(aligned_cycles(1).(marker_names{1}), 1);

   % Calculate means and SEMs for each marker
   for m = 1:length(marker_names)
       % Initialize 3D matrix for each marker (timepoints × dimensions × cycles)
       all_data = zeros(n_timepoints, 3, n_cycles);
       % Stack all cycles for this marker
       for cycle = 1:n_cycles
           this_d = aligned_cycles(cycle).(marker_names{m});
           if isempty(this_d), continue, end
           all_data(:,:,cycle) = aligned_cycles(cycle).(marker_names{m});
       end

       % Calculate mean and SEM
       marker_mean = mean(all_data, 3);
       marker_sem = std(all_data, 0, 3) / sqrt(n_cycles);

       % Calculate mean and SEM
       means.(marker_names{m}) = marker_mean;
       sems.(marker_names{m}) = marker_sem;
   end
   
   % Plot
   figure;
   hold on;
   for m = 1:length(marker_names)
       % Plot mean points

       scatter(means.(marker_names{m})(:,1), ...
               means.(marker_names{m})(:,3), ...
               50, params.marker_colors{m}, 'filled');
       
       % % Optional: Add error bars
       % for i = 1:size(means.(marker_names{m}), 1)
       %     scatter3([means.(marker_names{m})(i,1) - sems.(marker_names{m})(i,1), ...
       %           means.(marker_names{m})(i,1) + sems.(marker_names{m})(i,1)], ...
       %          [means.(marker_names{m})(i,2), means.(marker_names{m})(i,2)], ...
       %          [means.(marker_names{m})(i,3), means.(marker_names{m})(i,3)], ...
       %          'Color', params.marker_colors{m}, 'LineWidth', 1);
       % end
   end
   
   xlabel(params.plot_dimensions{1});
   ylabel(params.plot_dimensions{2});
   grid on;
   axis equal;
   view(45, 30);
   legend(marker_names, 'Location', 'best');
   title('Mean Gait Cycle');
   hold off;
end


% Helper functions
function merged = merge_structs(default_struct, input_struct)
   merged = default_struct;
   if ~isempty(input_struct)
       fields = fieldnames(input_struct);
       for i = 1:length(fields)
           merged.(fields{i}) = input_struct.(fields{i});
       end
   end
end

function valid = is_valid_segment(segment_idx, params)
   valid = length(segment_idx) >= round(params.min_segment_duration * params.sampling_rate);
end

function v_norm = normalize_vector(v)
   v_norm = v / (norm(v) + eps);
end

function sagittal_pos = compute_sagittal_position(marker)
   sagittal_pos = sqrt(marker(:,1).^2 + marker(:,3).^2);
end