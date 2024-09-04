% Figure for plotting density maps per condition
%% Load data
clear
close all
global GC

load(GC.filename_analysis, 'analysisstruct');

% Load predictions
load(GC.filename_predictions, 'predictions', 'animal_condition_identifier');

% Load ratception structure
load(GC.filename_ratception, 'ratception_struct');

%%
% change to animal_condition_identifier
long_animal_frames_identifier = repelem(animal_condition_identifier, GC.repfactor);

animal_list_used_after_analysis =  long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});

% animal_list_used_after_analysis = animal_condition_identifier;


% for inter animal analysis:
animal_list = unique(animal_list_used_after_analysis, 'stable');
cond_inds = zeros(1,length(analysisstruct.condition_inds)); % sorting per animal
for iid = 1:length(animal_list)
    animal_ID = animal_list{iid};
    idx = ismember(animal_list_used_after_analysis, animal_ID);
    cond_inds(idx) = iid;
end

% analysisstruct.condition_inds = cond_inds;

% conditions ie. actua conditions, concatenating animals in same cond
condition_inds = cond_inds; % sorting per condition
for iff = 1:length(animal_list_used_after_analysis)
     animal_ID = animal_list_used_after_analysis{iff};
    condition_inds(iff) = double(endsWith(animal_ID, '_F') +1);
   
end

%%


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