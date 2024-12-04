function movement_features = analyze_movement_patterns(data)
    markers = data;
    movement_features = struct();
    
    % Primary analysis pairs focusing on affected limb
    pain_related_pairs = {
        {'KneeL', 'KneeR'},...      % Direct comparison of affected vs unaffected knee
        {'AnkleL', 'AnkleR'},...    % Compensatory ankle movements
        {'HindpawL', 'HindpawR'},...% Weight bearing differences
        {'ForepawL', 'ForepawR'}    % Forepaw comparison
    };
    
    % Compensatory movement pairs
    compensatory_pairs = {
        {'SpineM', 'KneeL'},...     % Spine compensation for left leg
        {'ShoulderL', 'KneeL'},...  % Front-hind coordination on affected side
        {'ShoulderR', 'KneeL'}      % Cross-body compensation
    };
    
    % Analyze pain-related asymmetry
    for i = 1:length(pain_related_pairs)
        affected = markers.(pain_related_pairs{i}{1});
        unaffected = markers.(pain_related_pairs{i}{2});
        
        % Compute asymmetry index (modified to highlight affected side)
        asym_idx = compute_pain_asymmetry(affected, unaffected);
        movement_features.pain_asymmetry.(sprintf('asym_%s_%s', pain_related_pairs{i}{1}, pain_related_pairs{i}{2})) = asym_idx;
        
        % Phase analysis for compensatory patterns
        [phase_diff, coherence] = compute_phase_relationship(affected, unaffected);
        movement_features.compensation.(sprintf('phase_%s_%s', pain_related_pairs{i}{1}, pain_related_pairs{i}{2})) = phase_diff;
        movement_features.compensation.(sprintf('coh_%s_%s', pain_related_pairs{i}{1}, pain_related_pairs{i}{2})) = coherence;
    end
    
    % Analyze compensatory coordination
    for i = 1:length(compensatory_pairs)
        joint1 = markers.(compensatory_pairs{i}{1});
        joint2 = markers.(compensatory_pairs{i}{2});
        
        [phase_diff, coherence] = compute_phase_relationship(joint1, joint2);
        movement_features.compensation.(sprintf('coord_%s_%s', compensatory_pairs{i}{1}, compensatory_pairs{i}{2})) = phase_diff;
    end
end

function asym_idx = compute_pain_asymmetry(affected, unaffected)
    % Positive values indicate greater affected side movement (potential guarding)
    % Negative values indicate compensation by unaffected side
    diff = affected - unaffected;
    % mean_val = (abs(affected) + abs(unaffected)) / 2;
    % asym_idx = diff ./ (mean_val + eps);
    asym_idx = diff;
end

function [phase_diff, coherence] = compute_phase_relationship(signal1, signal2)
    hilbert1 = hilbert(detrend(signal1));
    hilbert2 = hilbert(detrend(signal2));
    phase1 = unwrap(angle(hilbert1));
    phase2 = unwrap(angle(hilbert2));
    phase_diff = phase1 - phase2;
    [Cxy, ~] = mscohere(signal1, signal2);
    coherence = mean(Cxy);
end