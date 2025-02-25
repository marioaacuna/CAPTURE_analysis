
logger('Starting creating of netowrk and test it', 'INFO');
%{
 animal_list = {'326', '327', '328', '330', '332_training', '332', '334', '335', '336'};
 
%}

clc, clear

GC = general_configs(); % load general configurations

%% Load Data
logger('Loading data', 'INFO');
% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');

% Load predictions
load(GC.filename_predictions, 'predictions', 'animal_condition_identifier');

% Load ratception structure
% load(GC.filename_ratception, 'ratception_struct');

x = cat(2, analysisstruct.jt_features, analysisstruct.extra_jt_features);

long_animal_frames_identifier = repelem(animal_condition_identifier,3);
animal_list_used_after_analysis =  long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});

%%
% Find indices of strings ending with 'F' or 'S'
valid_idx = cellfun(@(x) endsWith(x, {'F', 'S'}), animal_list_used_after_analysis);

% Extract corresponding rows from x and get labels
x1 = x(valid_idx, :);
l = cellfun(@(x) x(end), animal_list_used_after_analysis(valid_idx), 'UniformOutput', false);

% downsample
x_sub = x1(1:10:end, :) ;
l_sub = l(1:10:end,:);
% 
% % or not 
% x_sub = x1;
% l_sub = l;

%%
% Assuming l_sub is a cell array of labels
% Convert cell array of labels to categorical or string array if needed
if iscell(l_sub)
    labels = categorical(l_sub);
else
    labels = l_sub;
end

% Create variable names for features
feature_names = arrayfun(@(x) sprintf('Feature_%d', x), 1:66, 'UniformOutput', false);

% Create table in one go
data_table = array2table(x_sub, 'VariableNames', feature_names);
data_table.Label = categorical(labels);
%%
%l_cat = categorical(ismember(l_sub, 'F'))';
%tbl = array2table(x_sub, l_cat);
tbl = data_table;
tbl = splitvars(tbl);
classNames  = categories(tbl{:,'Label'});
numObservations = height(tbl);
numObservationsTrain = floor(0.8*numObservations);
numObservationsValidation = floor(0.1*numObservations);
numObservationsTest = numObservations - numObservationsTrain - numObservationsValidation;

idx = randperm(numObservations);
idxTrain = idx(1:numObservationsTrain);
idxValidation = idx(numObservationsTrain+1:numObservationsTrain+numObservationsValidation);
idxTest = idx(numObservationsTrain+numObservationsValidation+1:end);

numFeatures = size(x_sub,2);
numClasses = numel(classNames);


tblTrain = tbl(idxTrain,:);
tblValidation = tbl(idxValidation,:);
tblTest = tbl(idxTest,:);

%%
numGPUs = gpuDeviceCount("available");
miniBatchSize = 256*numGPUs;
options = trainingOptions("adam", ...
    MiniBatchSize=miniBatchSize, ...
    ValidationData=tblValidation, ...
    Shuffle="every-epoch", ...
    Plots="training-progress", ...
    Metrics="accuracy", ...
    MaxEpochs=100,...
    ValidationFrequency= 50,...
    Verbose=false,...
    ExecutionEnvironment="gpu");
% % 
layers = [
    featureInputLayer(numFeatures,Normalization="none")
    fullyConnectedLayer(50, "WeightsInitializer","glorot", "BiasInitializer","zeros", "WeightLearnRateFactor",1,"WeightL2Factor",1,BiasLearnRateFactor=1,BiasL2Factor=0)
    selfAttentionLayer(4,256,"Name","selfattention")
    reluLayer("Name","relu")
    fullyConnectedLayer(numClasses,"Name","fc_1")
    softmaxLayer("Name","softmax")];

%%
% n = trainnet(tbl,net_7.Layers,"crossentropy",options);
n = trainnet(tbl,layers,"crossentropy",options);

%% test
net = n;
labelName  = 'Label';
scores = minibatchpredict(net,tblTest(:,1:end-1),MiniBatchSize=16);
YPred = scores2label(scores,classNames);

YTest = tblTest{:,labelName};
accuracy = sum(YPred == YTest)/numel(YTest);
fprintf('Accuracy: %s\n', num2str(accuracy))
cm =confusionmat(YTest, YPred);
figure, confusionchart(cm, classNames, 'Normalization','row-normalized','Normalization','row-normalized', 'DiagonalColor', 'b', 'OffDiagonalColor', 'b', 'GridVisible', 'off')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Train a net based on all classes
% Find indices of strings ending with 'F' or 'S'
valid_idx = ones(size(animal_list_used_after_analysis,1),1);

% Extract corresponding rows from x and get labels
x1 = x(:, :);
l = cellfun(@(x) x(end), animal_list_used_after_analysis(:,1), 'UniformOutput', false);

% downsample
x_sub = x1(1:10:end, :) ;
l_sub = l(1:10:end,:);
% 
% % or not 
% x_sub = x1;
% l_sub = l;

%%
% Assuming l_sub is a cell array of labels
% Convert cell array of labels to categorical or string array if needed
if iscell(l_sub)
    labels = categorical(l_sub);
else
    labels = l_sub;
end

% Create variable names for features
feature_names = arrayfun(@(x) sprintf('Feature_%d', x), 1:66, 'UniformOutput', false);

% Create table in one go
data_table = array2table(x_sub, 'VariableNames', feature_names);
data_table.Label = categorical(labels);
%%
%l_cat = categorical(ismember(l_sub, 'F'))';
%tbl = array2table(x_sub, l_cat);
tbl = data_table;
tbl = splitvars(tbl);
classNames_all  = categories(tbl{:,'Label'});
numObservations = height(tbl);
numObservationsTrain_all = floor(0.7*numObservations);
numObservationsValidation_all = floor(0.2*numObservations);
numObservationsTest_all = numObservations - numObservationsTrain_all - numObservationsValidation_all;

idx = randperm(numObservations);
idxTrain = idx(1:numObservationsTrain_all);
idxValidation = idx(numObservationsTrain+1:numObservationsTrain+numObservationsValidation_all);
idxTest = idx(numObservationsTrain+numObservationsValidation_all+1:end);

numFeatures = size(x_sub,2);
numClasses_all = numel(classNames_all);


tblTrain_all = tbl(idxTrain,:);
tblValidation_all = tbl(idxValidation,:);
tblTest_all = tbl(idxTest,:);

%%
numGPUs = gpuDeviceCount("available");
miniBatchSize = 256*numGPUs;
options = trainingOptions("adam", ...
    MiniBatchSize=miniBatchSize, ...
    ValidationData=tblValidation_all, ...
    Shuffle="every-epoch", ...
    Plots="training-progress", ...
    Metrics="accuracy", ...
    MaxEpochs=100,...
    ValidationFrequency= 50,...
    Verbose=false,...
    ExecutionEnvironment="gpu");
% % 
layers = [
    featureInputLayer(numFeatures,Normalization="none")
    fullyConnectedLayer(50, "WeightsInitializer","glorot", "BiasInitializer","zeros", "WeightLearnRateFactor",1,"WeightL2Factor",1,BiasLearnRateFactor=1,BiasL2Factor=0)
    selfAttentionLayer(4,256,"Name","selfattention")
    reluLayer("Name","relu")
    fullyConnectedLayer(numClasses_all,"Name","fc_1")
    softmaxLayer("Name","softmax")];

%%
% n = trainnet(tbl,net_7.Layers,"crossentropy",options);
n_all = trainnet(tbl,layers,"crossentropy",options);

%% test
net = n_all;
labelName  = 'Label';
scores = minibatchpredict(net,tblTest_all(:,1:end-1),MiniBatchSize=16);
YPred = scores2label(scores,classNames_all);

YTest_all = tblTest_all{:,labelName};
accuracy_all = sum(YPred == YTest_all)/numel(YTest_all);
fprintf('Accuracy: %s\n', num2str(accuracy_all))
cm_all =confusionmat(YTest_all, YPred);
figure, confusionchart(cm_all, classNames_all, 'Normalization','row-normalized','Normalization','row-normalized', 'DiagonalColor', 'b', 'OffDiagonalColor', 'b', 'GridVisible', 'off')
%%