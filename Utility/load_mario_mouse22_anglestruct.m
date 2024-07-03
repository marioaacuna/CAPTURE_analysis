function anglestruct = load_mario_mouse22_anglestruct()

%% these are the XYZ reference frame angles -- unfortunately you can't use them for appendages, as their orientation will change
anglestruct.saggital_names = {'head_sagg','neck_sagg','spine_sagg'};
anglestruct.saggital_pairs =  {[2,3],[3,4],[4,5]}; %head, neck, spine angles , look in the z-y plane
%saggital_include = [1 1 1];

%transverse/overhead
anglestruct.transverse_names = {'head_trans','neck_trans','spine_trans','KneeL_trans','KneeR_trans',...
    'shouldl_trans','shouldr_trans','offset1_trans','offset2_trans'};
anglestruct.transverse_pairs =  {[2,3],[3,4],[4,5],[5,6],[5,7],[4 8 ], [4 9 ]}; %head, neck, spine angles , look in the z-y plane
%transverse_include = [1 1 1 1 1 0 0 1 1];

%coronal/along spine (front view)
anglestruct.coronal_names = {'head_coronal','KneeL_coronal','KneeR_coronal','shouldl_coronal','shouldr_coronal','offset1_coronal','offset2_coronal'};
anglestruct.coronal_pairs =  {[1,3],[5,6],[5,7],[4 8], [4 9]}; %head, neck, spine angles , look in the z-y plane
%coronal_include = [1 1 1 0 0 1 1];

%% angles to include for mai tsne
anglestruct.include_angles = {'head_sagg','neck_sagg','spine_sagg','head_trans','neck_trans',...
    'spine_trans','offset1_trans','offset2_trans','head_coronal','KneeL_coronal','KneeR_coronal',...
    'offset1_coronal','offset2_coronal','KneeL_pitch','KneeL_yaw','KneeR_pitch','KneeR_yaw'};% TO ADD: HIP YAW
%jointangle_struct.KneeL_pitch,jointangle_struct.KneeL_yaw

%alljt names
% allangles_names = {'lelbow_all','larm_all','relbow_all','rarm_all','lknee_all','lshin_all','rknee_all','rshin_all'};
% allangles_pairs =  {[8,10],[10,11],[9,12],[12,13],[6,14],[14,15],[7,16],[16,17]}; %head, neck, spine angles , look in the z-y plane
% all_include = zeros(1,numel(allangles_names));


%% get the various
anglestruct.segment_pairs = {{'EarR','EarL'},{'Snout','EarR'},{'EarR','SpineF'},{'SpineF','SpineM'} ,...%1-4
    {'Tail_base_','SpineM'},{'Tail_base_','KneeL'},{'Tail_base_','KneeR'},... %5-7
    {'SpineF','ShoulderL'},{'SpineF','ShoulderR'},... %8,9
    {'ElbowL','WristL'},{'ElbowR','WristR'},...%10-13
    {'KneeL','AnkleL'},{'KneeR','AnkleR'}...%14-17
    {'ForepawL','ElbowL'},{'ElbowL','ShoulderL'},...
    {'ForepawR','ElbowR'},{'ElbowR','ShoulderR'},...
    {'HindpawL','AnkleL'},{'HindpawR','AnkleR'},{'Tail_base_','Tail_mid_'}}; %18,19


%hip angles to back of the spine
anglestruct.planar_trios{1}.plane = {{'SpineM','Tail_base_'},{'zvector','zvector'}};
anglestruct.planar_trios{1}.vector = {'Tail_base_','KneeL'};
anglestruct.planar_trios{1}.name1 = 'hipl_pitch';
anglestruct.planar_trios{1}.name2 = 'hipl_yaw';
anglestruct.planar_trios{1}.namesuse = [1 2];

anglestruct.planar_trios{2}.plane = {{'SpineM','Tail_base_'},{'zvector','zvector'}};
anglestruct.planar_trios{2}.vector = {'Tail_base_','KneeR'};
anglestruct.planar_trios{2}.name1 = 'hipr_pitch';
anglestruct.planar_trios{2}.name2 = 'hipr_yaw';
anglestruct.planar_trios{2}.namesuse = [1 2];
%knees in the plane defined by the hips and the back of the
%spine line
anglestruct.planar_trios{3}.plane = {{'Tail_base_','SpineM'},{'KneeR','KneeL'}};
anglestruct.planar_trios{3}.vector = {'KneeL','AnkleL'};
anglestruct.planar_trios{3}.name1 = 'kneel_yaw';
anglestruct.planar_trios{3}.name2 = 'kneel_pitch';
anglestruct.planar_trios{3}.namesuse = [1 2];

anglestruct.planar_trios{4}.plane = {{'Tail_base_','SpineM'},{'KneeL','KneeR'}};
anglestruct.planar_trios{4}.vector = {'KneeR','AnkleR'};
anglestruct.planar_trios{4}.name1 = 'kneer_yaw';
anglestruct.planar_trios{4}.name2 = 'kneer_pitch';
anglestruct.planar_trios{4}.namesuse = [1 2];

%shin line, relative to the plane of knee and hips, with axes
%chosen to get angles >0 to avoid ambiguity/0 crossing
anglestruct.planar_trios{5}.plane = {{'KneeR','KneeL'},{'KneeL','AnkleL'}};
anglestruct.planar_trios{5}.vector = {'HindpawL','AnkleL'};
anglestruct.planar_trios{5}.name1 = 'shinl_yaw';
anglestruct.planar_trios{5}.name2 = '~';
anglestruct.planar_trios{5}.namesuse = [1 ];

anglestruct.planar_trios{6}.plane = {{'KneeL','KneeR'},{'KneeR','AnkleR'}};
anglestruct.planar_trios{6}.vector = {'HindpawR','AnkleR'};
anglestruct.planar_trios{6}.name1 = 'shinr_yaw';
anglestruct.planar_trios{6}.name2 = '~';
anglestruct.planar_trios{6}.namesuse = [1 ];

anglestruct.planar_trios{7}.plane = {{'KneeR','KneeL'},{'KneeL','AnkleL'}};
anglestruct.planar_trios{7}.vector = {'HindpawL','AnkleL'};
anglestruct.planar_trios{7}.name1 = '~';
anglestruct.planar_trios{7}.name2 = 'shinl_pitch';
anglestruct.planar_trios{7}.namesuse = [ 2];

anglestruct.planar_trios{8}.plane = {{'KneeL','KneeR'},{'KneeR','AnkleR'}};
anglestruct.planar_trios{8}.vector = {'HindpawR','AnkleR'};
anglestruct.planar_trios{8}.name1 = '~';
anglestruct.planar_trios{8}.name2 = 'shinr_pitch';
anglestruct.planar_trios{8}.namesuse = [ 2];
% 
% 
% % arm angles
% 
anglestruct.planar_trios{9}.plane = {{'ShoulderR','ShoulderL'},{'ShoulderL','ElbowL'}};
anglestruct.planar_trios{9}.vector = {'ElbowL','ForepawL'};
anglestruct.planar_trios{9}.name1 = 'arml_yaw';
anglestruct.planar_trios{9}.name2 = '~';
anglestruct.planar_trios{9}.namesuse = [1 ];

anglestruct.planar_trios{10}.plane = {{'ShoulderL','ShoulderR'},{'ShoulderR','ElbowR'}};
anglestruct.planar_trios{10}.vector = {'ElbowR','ForepawR'};
anglestruct.planar_trios{10}.name1 = 'armr_yaw';
anglestruct.planar_trios{10}.name2 = '~';
anglestruct.planar_trios{10}.namesuse = [1 ];

anglestruct.planar_trios{11}.plane = {{'ShoulderR','ShoulderL'},{'ShoulderL','ElbowL'}};
anglestruct.planar_trios{11}.vector = {'ElbowL','ForepawL'};
anglestruct.planar_trios{11}.name1 = '~';
anglestruct.planar_trios{11}.name2 = 'arml_pitch';
anglestruct.planar_trios{11}.namesuse = [ 2];

anglestruct.planar_trios{12}.plane = {{'ShoulderL','ShoulderR'},{'ShoulderR','ElbowR'}};
anglestruct.planar_trios{12}.vector = {'ElbowR','ForepawR'};
anglestruct.planar_trios{12}.name1 = '~';
anglestruct.planar_trios{12}.name2 = 'armr_pitch';
anglestruct.planar_trios{12}.namesuse = [ 2];

%elbow angles
% elbow angle
% WHICH ARE YAW AND PITCH?
anglestruct.planar_trios{13}.plane = {{'SpineM','SpineF'},{'ShoulderR','ShoulderL'}};
anglestruct.planar_trios{13}.vector = {'ElbowL','ShoulderL'};
anglestruct.planar_trios{13}.name1 = 'elbowl_yaw';
anglestruct.planar_trios{13}.name2 = 'elbowl_pitch';
anglestruct.planar_trios{13}.namesuse = [1 2];

anglestruct.planar_trios{14}.plane = {{'SpineM','SpineF'},{'ShoulderL','ShoulderR'}};
anglestruct.planar_trios{14}.vector = {'ElbowR','ShoulderR'};
anglestruct.planar_trios{14}.name1 = 'elbowr_yaw';
anglestruct.planar_trios{14}.name2 = 'elbowr_pitch';
anglestruct.planar_trios{14}.namesuse = [1 2];

%% new elbow yaw and pitch -suggested by LLM Claude
% 
% % Left elbow angle
% anglestruct.planar_trios{13}.plane = {{'ShoulderL','ElbowL'},{'ShoulderR','ShoulderL'}};
% anglestruct.planar_trios{13}.vector = {'ElbowL','ForepawL'};
% anglestruct.planar_trios{13}.name1 = 'elbowl_pitch';
% anglestruct.planar_trios{13}.name2 = 'elbowl_yaw';
% anglestruct.planar_trios{13}.namesuse = [1 2];
% 
% % Right elbow angle
% anglestruct.planar_trios{14}.plane = {{'ShoulderR','ElbowR'},{'ShoulderL','ShoulderR'}};
% anglestruct.planar_trios{14}.vector = {'ElbowR','ForepawR'};
% anglestruct.planar_trios{14}.name1 = 'elbowr_pitch';
% anglestruct.planar_trios{14}.name2 = 'elbowr_yaw';
% anglestruct.planar_trios{14}.namesuse = [1 2];

%%
% anglestruct.planar_trios{13}.plane = {{'SpineM','SpineF'},{'ShoulderR','ShoulderL'}};
% anglestruct.planar_trios{13}.vector = {'ShoulderL','ElbowL'};
% anglestruct.planar_trios{13}.name1 = 'elbowl_yaw';
% anglestruct.planar_trios{13}.name2 = 'elbowl_pitch';
% anglestruct.planar_trios{13}.namesuse = [1 2];
% 
% anglestruct.planar_trios{14}.plane = {{'SpineM','SpineF'},{'ShoulderL','ShoulderR'}};
% anglestruct.planar_trios{14}.vector = {'ShoulderR','ElbowR'};
% anglestruct.planar_trios{14}.name1 = 'elbowr_yaw';
% anglestruct.planar_trios{14}.name2 = 'elbowr_pitch';
% anglestruct.planar_trios{14}.namesuse = [1 2];


