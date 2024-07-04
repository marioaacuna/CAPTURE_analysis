% plot density map
% conditions is coming from :
%   Extract unique animal IDs and conditions
%   [animal_ids, ~, animal_indices] = unique(cellfun(@(x) x(1:end-2), frame_identifiers, 'UniformOutput', false));
%   conditions = cellfun(@(x) x(end), frame_identifiers, 'UniformOutput', false);

conds = unique(conditions);


for icond = 1:2
    fighand_in =figure(444 + icond);
    set(fighand_in,'Color','w')
    title(this_cond)

    this_cond = conds(icond);
    Zvals =analysisstruct.zValues(contains(conditions, this_cond),:);
    % subplot(1,2,icond)
     plotdensitymaps({Zvals},1,fighand_in,analysisstruct.params.density_width,...
        max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor,analysisstruct.params.density_res);
    % this will normalize it to all values
end