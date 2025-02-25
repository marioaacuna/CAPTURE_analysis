logger('Loading data', 'INFO');
% Load analysis structure
load(GC.filename_analysis, 'analysisstruct');

% Load predictions
load(GC.filename_predictions, 'predictions', 'animal_condition_identifier');
%%
long_animal_frames_identifier = repelem(animal_condition_identifier,3);
animal_list_used_after_analysis =  long_animal_frames_identifier(analysisstruct.frames_with_good_tracking{1});
%%
x = analysisstruct.zValues;
l = {};
for il = 1:length(animal_list_used_after_analysis)
    thisl  = animal_list_used_after_analysis{il};
    l{il} = thisl(end);
end

% downsample
x_sub = x(1:10:end, :) ;
l_sub = l(1, 1:10:end);

%% select only F and S

x = analysisstruct.zValues;
x1 = [];
l = {};
counti = 0;
for il = 1:length(animal_list_used_after_analysis)
    thisl  = animal_list_used_after_analysis{il};
    if endsWith(thisl, 'F') || endsWith(thisl, 'S')
        counti = counti+1;
        x1(counti,:) = x(il,:);
        l{counti} = thisl(end);

    else
        continue
    end
end


% downsample
x_sub = x1(1:10:end, :) ;
l_sub = l(1, 1:10:end);

%% figure

%%
clc
conds = {'F', 'S'};
[pred, validationAccuracy, partitionedModel, cm] = get_predictions(analysisstruct.zValues, animal_list_used_after_analysis, conds);

accuracy = sum(diag(cm.NormalizedValues)) / sum(cm.NormalizedValues(:));
fprintf('Val Accuracy: %s\n', num2str(accuracy))

%%

function [validationPredictions, validationAccuracy, partitionedModel, cm] = get_predictions(x,animal_list_used_after_analysis,  conds)
rng default
% x = analysisstruct.zValues;
x1 = [];
l = {};
counti = 0;
for il = 1:length(animal_list_used_after_analysis)
    thisl  = animal_list_used_after_analysis{il};
    if endsWith(thisl, conds(1)) || endsWith(thisl, conds(2))
        counti = counti+1;
        x1(counti,:) = x(il,:);
        l{counti} = thisl(end);

    else
        continue
    end
end

figure
gscatter(x1(:,1), x1(:,2), l', [1,0,0;0,0,1])


% convert names
trainingData = x1;
responseData = l';
% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2'});

predictorNames = {'column_1', 'column_2'};
predictors = inputTable(:, predictorNames);
response = responseData;
isCategoricalPredictor = [false, false];
classNames = unique(response);

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
classificationKNN = fitcknn(...
    predictors, ...
    response, ...
    'Distance', 'Euclidean', ...
    'Exponent', [], ...
    'NumNeighbors', 2000, ...
    'DistanceWeight', 'Equal', ...
    'Standardize', false, ...
    'ClassNames', classNames);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
knnPredictFcn = @(x) predict(classificationKNN, x);
trainedClassifier.predictFcn = @(x) knnPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationKNN = classificationKNN;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2024b.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  [yfit,scores] = c.predictFcn(X) \nreplace ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 2 columns because this model was trained using 2 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');



% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2'});

predictorNames = {'column_1', 'column_2'};
predictors = inputTable(:, predictorNames);
response = responseData;
isCategoricalPredictor = [false, false];

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationKNN, 'KFold', 5);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');


%% Plot confusion
figure,
cm = confusionchart(response',validationPredictions, 'Normalization','row-normalized');


end