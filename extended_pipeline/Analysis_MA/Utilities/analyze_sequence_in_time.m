function [stats] = analyze_sequence_in_time(clusters_seq)
% Constants
MINUTES_RECORDED = 30;  % Total recording time in minutes
ORIGINAL_FPS = 100;     % Original recording rate

% Calculate temporal scaling
total_recording_frames = MINUTES_RECORDED * 60 * ORIGINAL_FPS;
scaling_factor = total_recording_frames / (length(clusters_seq));
seconds_per_index = (MINUTES_RECORDED * 60) / length(clusters_seq);

% Find state changes
change_points = find(diff(clusters_seq) ~= 0);

% Calculate dwell times without problematic concatenation
dwell_frames = zeros(length(change_points) + 1, 1);
dwell_frames(1) = change_points(1);  % First dwell
for i = 2:length(change_points)
    dwell_frames(i) = change_points(i) - change_points(i-1);
end
dwell_frames(end) = length(clusters_seq) - change_points(end);  % Last dwell

% Convert to actual time
dwell_times = dwell_frames * seconds_per_index;

% Calculate meaningful statistics
stats = struct();
stats.mean_dwell_time_sec = mean(dwell_times);
stats.median_dwell_time_sec = median(dwell_times);
stats.transitions_per_minute = length(change_points) / MINUTES_RECORDED;
stats.dwell_times = dwell_times;

% State transition analysis
transitions = zeros(length(change_points), 2);
for i = 1:length(change_points)
    transitions(i,:) = [clusters_seq(change_points(i)), ...
        clusters_seq(change_points(i)+1)];
end

% Analyze transition patterns
[unique_trans, ~, ic] = unique(transitions, 'rows');
trans_counts = accumarray(ic, 1);
[sorted_counts, idx] = sort(trans_counts, 'descend');

stats.transition_patterns = struct(...
    'transitions', unique_trans(idx,:), ...
    'counts', sorted_counts, ...
    'frequency_per_minute', sorted_counts / MINUTES_RECORDED);

% Add temporal context
stats.temporal_context = struct(...
    'total_frames', length(clusters_seq), ...
    'transition_count', length(change_points), ...
    'transition_rate', length(change_points)/length(clusters_seq), ...
    'seconds_per_index', seconds_per_index);
end