% plot density map
% conditions is coming from :Script 1 line 222 etc
% conditions ie. actua conditions, concatenating animals in same cond
% condition_inds = cond_inds; % sorting per condition
% for iff = 1:length(animal_list_used_after_analysis)
%      animal_ID = animal_list_used_after_analysis{iff};
%     condition_inds(iff) = double(endsWith(animal_ID, '_F') +1);
% 
% end


conds = unique(condition_inds);


for icond = 1:2
    fighand_in =figure(444 + icond);
    set(fighand_in,'Color','w')
    

    this_cond = conds(icond);
    Zvals =analysisstruct.zValues(ismember(condition_inds, this_cond),:);
    % Zvals =analysisstruct.zValues(contains(condition_inds, this_cond),:);
    % subplot(1,2,icond)
     plotdensitymaps({Zvals},1,fighand_in,analysisstruct.params.density_width,...
        max(analysisstruct.zValues(:))*analysisstruct.params.expansion_factor,analysisstruct.params.density_res);
    % this will normalize it to all values
     title(this_cond)
end