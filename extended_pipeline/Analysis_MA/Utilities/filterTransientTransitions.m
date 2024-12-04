function [filtered_sequence, transition_info] = filterTransientTransitions(cluster_sequence)
    % Filter out A→B→A patterns while preserving transition information
    %
    % Input:
    %   cluster_sequence: Original sequence of cluster IDs
    % 
    % Outputs:
    %   filtered_sequence: Cleaned sequence with transient states removed
    %   transition_info: Structure containing transition analysis
    
    % Initialize filtered sequence
    filtered_sequence = cluster_sequence;
    sequence_length = length(cluster_sequence);
    
    % Track modifications for analysis
    modification_points = [];
    original_patterns = {};
    
    % First pass: Identify transient transitions
    i = 1;
    while i < sequence_length - 1
        current_state = filtered_sequence(i);
        next_state = filtered_sequence(i + 1);
        after_next = filtered_sequence(i + 2);
        
        if (current_state ~= next_state) && (current_state == after_next)
            % Store information about the modification
            modification_points = [modification_points; i];
            original_patterns{end+1} = [current_state, next_state, after_next];
            
            % Remove transient state
            filtered_sequence(i+1) = current_state;
            
            % Skip next position as we've already processed it
            i = i + 2;
        else
            i = i + 1;
        end
    end
    
    % Calculate transition statistics
    transition_info = struct();
    transition_info.total_modifications = length(modification_points);
    transition_info.modification_points = modification_points;
    transition_info.original_patterns = original_patterns;
    
    if ~isempty(modification_points)
        % Analyze temporal distribution of modifications
        transition_info.temporal_distribution = diff(modification_points);
        transition_info.modification_rate = length(modification_points) / sequence_length;
        
        % Analyze pattern frequencies
        [unique_patterns, ~, pattern_idx] = unique(cell2mat(original_patterns), 'rows');
        pattern_counts = accumarray(pattern_idx, 1);
        transition_info.pattern_frequencies = [unique_patterns, pattern_counts];
    end
end