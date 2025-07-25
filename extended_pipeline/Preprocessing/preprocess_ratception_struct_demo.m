function [ratception_struct,analysisstruct_app] = preprocess_ratception_struct_demo(aa,preprocessing_parameters,params)

repfactor = params.repfactor;

%% Ensure input is single precision
markernames = fieldnames(aa.predictions);
for i = 1:length(markernames)
    if isa(aa.predictions.(markernames{i}), 'double')
        aa.predictions.(markernames{i}) = single(aa.predictions.(markernames{i}));
    end
end

ratception_struct = [];
ratception_struct.fps = 300;
ratception_struct.markers_preproc = aa.predictions;
ratception_struct.markernames = markernames;
ratception_struct.markers = ratception_struct.markers_preproc;

ratception_structtemp = ratception_struct;
chunksize = 10^5;
totalsize = size(ratception_struct.markers_preproc.SpineF,1);
ratception_struct_full = [];
filelength = 0;

for rk = 1:ceil(totalsize/chunksize)
    disp(['chunk ', num2str(rk), ' from ', num2str(ceil(totalsize/chunksize))])
    ratception_struct_temp = ratception_structtemp;
    
    for ll = 1:numel(markernames)
        % Force single precision in repelem output
        temp_data = ratception_struct_temp.markers.(markernames{ll})(1+chunksize*(rk-1):min(chunksize*rk,totalsize),:);
        ratception_struct_temp.markers.(markernames{ll}) = single(repelem(temp_data, repfactor, 1));
    end
    ratception_struct_temp.markers_preproc = ratception_struct_temp.markers;

    ratception_struct_temppreproc = compute_preprocessed_mocapstruct_demo(ratception_struct_temp, preprocessing_parameters);

    if (rk == 1)
        ratception_struct_full = ratception_struct_temppreproc;
        filelength = size(ratception_struct_temppreproc.aligned_mean_position, 1);
    else
        for ll = 1:numel(markernames)
            % Ensure concatenation preserves single precision
            ratception_struct_full.markers_preproc.(markernames{ll}) = single(cat(1, ...
                ratception_struct_full.markers_preproc.(markernames{ll}), ...
                ratception_struct_temppreproc.markers_preproc.(markernames{ll})));
            
            ratception_struct_full.markers_aligned_preproc.(markernames{ll}) = single(cat(1, ...
                ratception_struct_full.markers_aligned_preproc.(markernames{ll}), ...
                ratception_struct_temppreproc.markers_aligned_preproc.(markernames{ll})));
        end

        % Force single precision on large arrays
        ratception_struct_full.aligned_rotation_matrix = single(cat(3, ...
            ratception_struct_full.aligned_rotation_matrix, ...
            ratception_struct_temppreproc.aligned_rotation_matrix));
        
        ratception_struct_full.aligned_mean_position = single(cat(1, ...
            ratception_struct_full.aligned_mean_position, ...
            ratception_struct_temppreproc.aligned_mean_position));
    end
    
    % Clear temporary structures
    clear ratception_struct_temppreproc temp_data
end

ratception_struct = ratception_struct_full;
clear ratception_struct_full

end