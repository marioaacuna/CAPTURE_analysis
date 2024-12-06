function all_cycles = analyze_gait_cycles_3d(markers_aligned, walking_bouts, params)
% The 3 dimensions in aligned_cycles(1).wrist_l are:
% 
% Heading direction (forward/backward movement)
% Mediolateral direction (left/right movement)
% Vertical direction (up/down movement)   
% 

% % Initialize markers and parameters
   params = initialize_params(params);
   markers = initialize_markers(markers_aligned, params);
   
   % Process walking bouts
   all_cycles = process_walking_bouts(markers, walking_bouts, params);
   
   % Align and analyze cycles
   % aligned_cycles = align_cycles_3d(all_cycles, params);
   
   % Visualize results
   if params.plot_results
       plot_3d_cycles(aligned_cycles, params);
   end
end


function params = initialize_params(params)
   default_params = struct(...
       'sampling_rate', 100,...       % Data acquisition frequency (Hz)
       'min_segment_duration', 0.9,...    % Minimum duration of valid walking bout (seconds)
       'min_cycle_length', 0.25,...     % Minimum duration of gait cycle (seconds)
       'max_cycle_length', 2,...     % Maximum duration of gait cycle (seconds)
       'window_size', 5,...            % Smoothing window size (samples)
       'cycle_window', [-0.5 0.5],...  % Time window for cycle normalization (seconds)
       'n_normalized_points', 200,... % Number of points after interpolation
       'interp_method', 'pchip',...    % Interpolation method (pchip/linear/spline)
       'cycle_detection_method', 'phase_space',... % Method for detecting cycles, phase_space, peak_detection, zero_crossing
       'alignment_method', 'fixed_window',... % Time normalization approach
       'plot_results', false,...        % Enable/disable plotting
       'plot_dimensions', {{'Time (AU)', 'Sagittal'}},... % Axis labels
       'marker_colors', {{[0.5 0 0.5], [1 0 0], [0 1 1], [1 1 0]}}); % Colors for each marker

   params = merge_structs(default_params, params);
end

function markers = initialize_markers(markers_aligned, params)
    markers = struct();
    for i = 1:length(params.markers_to_study)
        marker_name = params.markers_to_study{i};
        if isfield(markers_aligned, marker_name)
            markers.(marker_name) = markers_aligned.(marker_name);
        end
    end
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
       [cycles] = detect_cycles(norm_markers, params);
       
       all_cycles = store_valid_cycles(cycles, all_cycles, cycle_count);
       cycle_count = cycle_count + 1;
   end
end

function marker_cycles = detect_cycles(markers, params)
    marker_cycles = struct();
    marker_names = fieldnames(markers);
    for i = 1:length(marker_names)
        marker = markers.(marker_names{i});
        sagittal = compute_sagittal_position(marker);
        norm_sagittal = sagittal - mean(sagittal);
        
        % Find both positive and negative crossings
        pos_crossings = find(diff(sign(norm_sagittal)) > 0);
        % neg_crossings = find(diff(sign(norm_sagittal)) < 0);
        
        % Initialize storage for cycles
        cycles = {};
        
        % Extract cycles using positive crossings (could also use negative)
        for icr = 1:length(pos_crossings)-1
            start_idx = pos_crossings(icr);
            end_idx = pos_crossings(icr+1);
            
            % Extract cycle
            cycle = norm_sagittal(start_idx:end_idx);
            
            % Check if cycle length is valid
            if length(cycle) >= params.min_cycle_length*params.sampling_rate &&...
                    length(cycle) <= params.max_cycle_length*params.sampling_rate
                cycles{end+1} = cycle;
            end
        end
        
    
        % Align cycles to same length using interpolation
        n_points = params.n_normalized_points;  % e.g., 100 points
        aligned_cycles = zeros(length(cycles), n_points);
        
        for ic = 1:length(cycles)
            % Create time vectors for interpolation
            t_orig = linspace(0, 1, length(cycles{ic}));
            t_normalized = linspace(0, 1, n_points);
            
            % Interpolate to normalized length
            aligned_cycles(ic,:) = interp1(t_orig, cycles{ic}, t_normalized, 'pchip');
        end

        marker_cycles.(marker_names{i})= aligned_cycles;
    end
end

function [starts, ends] = detect_by_zero_crossing(markers, marker_names, params)
    % Initialize storage for each marker's crossings
    marker_crossings = struct();
    
    % Process each marker independently
    for i = 1:length(marker_names)
        if ~strcmp(marker_names{i}, 'spine_m')
            marker = markers.(marker_names{i});
            sagittal = compute_sagittal_position(marker);
            
            % Center the signal around its mean to ensure meaningful crossings
            centered_signal = sagittal - mean(sagittal);
            
            % Find zero crossings (positive slope)
            crossings = find(diff(sign(centered_signal)) > 0);
            
            % Store if we have at least one complete cycle
            if length(crossings) >= 2
                marker_crossings.(marker_names{i}).starts = crossings(1:end-1);
                marker_crossings.(marker_names{i}).ends = crossings(2:end);
            end
        end
    end
    
    % Return consolidated crossings
    [starts, ends] = consolidate_marker_cycles(marker_crossings);
end

function [marker_cycles] = detect_by_phase(markers, marker_names, params)
    marker_cycles = struct();
    
    for i = 1:length(marker_names)
        if ~strcmp(marker_names{i}, 'spine_m')
            marker = markers.(marker_names{i});
            sagittal = compute_sagittal_position(marker);
            
            % Compute velocity with improved smoothing
            velocity = smooth_velocity(sagittal, params.sampling_rate);
            
            % Normalize signals preserving zero-crossings
            norm_pos = sagittal - mean(sagittal);
            norm_vel = velocity - mean(velocity);
            
            % Compute phase angle
            phase = atan2(norm_vel, norm_pos);
            
            % Instead of looking for negative transitions in unwrapped phase,
            % look for complete rotations in the phase plane
            transitions = find(abs(diff(phase)) > pi);
            
            % Filter transitions to ensure they represent true cycles
            valid_transitions = filter_transitions(transitions, phase, params);
            
            if length(valid_transitions) >= 2
                marker_cycles.(marker_names{i}).starts = valid_transitions(1:end-1);
                marker_cycles.(marker_names{i}).ends = valid_transitions(2:end);


            end
        end
    end

    
    [starts, ends] = consolidate_marker_cycles(marker_cycles);
end

function valid = filter_transitions(transitions, phase, params)
    if isempty(transitions)
        valid = [];
        return
    end
    
    % Minimum samples between transitions
    min_samples = round(params.min_cycle_length * params.sampling_rate);
    
    % Find transitions that are sufficiently separated
    transition_gaps = diff(transitions);
    valid_gaps = transition_gaps >= min_samples;
    
    % Include first transition and those with valid gaps
    valid = transitions([true; valid_gaps]);
    
    % Verify direction of rotation at transitions
    for i = 1:length(valid)-1
        segment = phase(valid(i):valid(i+1));
        if ~is_valid_rotation(segment)
            valid(i+1) = [];
        end
    end
end

function valid = is_valid_rotation(phase_segment)
    % Unwrap just this segment
    unwrapped = unwrap(phase_segment);
    % Check if total phase change is approximately 2π
    phase_change = abs(unwrapped(end) - unwrapped(1));
    valid = abs(phase_change - 2*pi) < pi/2;
end
function v_smooth = smooth_velocity(signal, fs)
    % Compute velocity with Savitzky-Golay filtering
    window = min(21, length(signal)-1);
    if mod(window, 2) == 0
        window = window - 1;
    end
    v_smooth = sgolayfilt(gradient(signal) * fs, 3, window);
end

function [starts, ends] = consolidate_marker_cycles(marker_cycles)
    % Combine cycles from all markers while preserving individual timing
    all_starts = [];
    all_ends = [];
    
    marker_names = fieldnames(marker_cycles);
    for i = 1:length(marker_names)
        if isfield(marker_cycles.(marker_names{i}), 'starts')
            all_starts = [all_starts; marker_cycles.(marker_names{i}).starts];
            all_ends = [all_ends; marker_cycles.(marker_names{i}).ends];
        end
    end
    
    % Sort and pair corresponding starts/ends
    [starts, sort_idx] = sort(all_starts);
    ends = all_ends(sort_idx);
end

function [starts, ends] = detect_by_peaks(markers, marker_names, params)
     all_peaks = [];
    
    % Detect peaks for each marker separately
    for i = 1:length(marker_names)
        if ~strcmp(marker_names{i}, 'spine_m')
            marker = markers.(marker_names{i});
            sagittal = compute_sagittal_position(marker);
            
            % Find peaks for this marker
            [~, marker_peaks] = findpeaks(sagittal, ...
                'MinPeakDistance', round(params.min_cycle_length * params.sampling_rate));
            
            all_peaks = [all_peaks; marker_peaks];
        end
    end
    
    % Sort all detected peaks chronologically
    all_peaks = sort(all_peaks);
    
    % Remove peaks that are too close together (optional)
    min_separation = round(params.min_cycle_length * params.sampling_rate / 2);
    valid_peaks = [true; diff(all_peaks) >= min_separation];
    all_peaks = all_peaks(valid_peaks);
    
    % Create starts and ends
    starts = all_peaks(1:end-1);
    ends = all_peaks(2:end);
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

function all_cycles = store_valid_cycles(cycles, all_cycles, cycle_count)
   % Get valid cycles only
   valid_cycles = cycles;
   
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
       % all_data = zeros(n_timepoints, 3, n_cycles);
       % % Stack all cycles for this marker
       % for cycle = 1:n_cycles
       %     this_d = aligned_cycles(cycle).(marker_names{m});
       %     if isempty(this_d), continue, end
       %     all_data(:,:,cycle) = aligned_cycles(cycle).(marker_names{m});
       % end
       all_data = [];
       for cycle = 1:n_cycles
           this_d = aligned_cycles(cycle).(marker_names{m});
           all_data = [all_data;this_d];

       end

       % Calculate mean and SEM
       marker_mean = mean(all_data, 1);
       marker_sem = std(all_data, 0, 1) / sqrt(size(all_data,1));

       % Calculate mean and SEM
       means.(marker_names{m}) = marker_mean;
       sems.(marker_names{m}) = marker_sem;
   end
   
   % Plot
   figure;
   hold on;
   for m = 1:length(marker_names)
       % Plot mean points
       subplot(2,2,m)

       plot(means.(marker_names{m}));
       title(marker_names{m})
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