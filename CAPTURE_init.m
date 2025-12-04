addpath(genpath(pwd))
rmpath(genpath('./original_pipeline'))
global GC logger
GC = general_configs();
maxNumCompThreads(24)
disp('Capture init.ed')
