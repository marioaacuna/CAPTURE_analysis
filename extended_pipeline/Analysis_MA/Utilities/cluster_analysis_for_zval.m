function analysisstruct_out = cluster_analysis_for_zval(analysisstruct_in, zval_field, GC)
%CLUSTER_ANALYSIS_FOR_ZVAL Run clustering pipeline for a specific zValues field.
%   analysisstruct_out = cluster_analysis_for_zval(analysisstruct_in, zval_field, GC)
%   - Copies the input analysisstruct
%   - Removes all fields starting with 'zValues' except zval_field
%   - Renames the selected zValues field to 'zValues'
%   - Sets/refreshes clustering parameters from GC
%   - Runs compute_analysis_clusters_demo and returns the updated struct

% Safety checks
if nargin < 3 || isempty(GC)
    GC = general_configs(); %#ok<NASGU>
end

if ~ischar(zval_field) && ~isstring(zval_field)
    error('zval_field must be a string name of a field inside analysisstruct.');
end
zval_field = char(zval_field);

if ~isfield(analysisstruct_in, zval_field)
    error('Field "%s" not found in analysisstruct.', zval_field);
end

% Create a temp copy and prune zValues* fields except the one we want
analysisstruct_temp = analysisstruct_in; %#ok<NASGU>
fns = fieldnames(analysisstruct_temp);
z_like = startsWith(fns, 'zValues');
fns_to_remove = fns(z_like & ~strcmp(fns, zval_field));
if ~isempty(fns_to_remove)
    analysisstruct_temp = rmfield(analysisstruct_temp, fns_to_remove);
end

% Rename selected field to 'zValues' if needed
if ~strcmp(zval_field, 'zValues')
    analysisstruct_temp.zValues = analysisstruct_temp.(zval_field);
    analysisstruct_temp = rmfield(analysisstruct_temp, zval_field);
else
    % Ensure existence
    analysisstruct_temp.zValues = analysisstruct_temp.zValues;
end

% Ensure clustering parameters exist (refresh from GC where applicable)
try
    analysisstruct_temp.params.density_res       = GC.density_res;      
    analysisstruct_temp.params.density_width     = GC.density_width;    
    analysisstruct_temp.params.expansion_factor  = GC.expansion_factor; 
    analysisstruct_temp.params.density_threshold = GC.density_threshold;
catch
    % If GC missing any field, ignore and keep existing
end

% Matched conditions setup if missing
if isfield(analysisstruct_temp, 'condition_inds') && ~isempty(analysisstruct_temp.condition_inds)
    uniq = unique(analysisstruct_temp.condition_inds);
    analysisstruct_temp.matchedconds     = {uniq};
    analysisstruct_temp.conditions_to_run = uniq;
else
    % Fallback to single condition
    analysisstruct_temp.matchedconds      = {1};
    analysisstruct_temp.conditions_to_run = 1;
end

% Run clustering
fprintf('%% INIT clustering for %s %%\n', zval_field);
params = struct('reorder', 1);
analysisstruct_temp = compute_analysis_clusters_demo(analysisstruct_temp, params); % adds annot_reordered, etc.
fprintf('%% Done clustering for %s %%\n', zval_field);

% Return
analysisstruct_out = analysisstruct_temp;

end
