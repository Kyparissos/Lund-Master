%% 1. Load mnist data
[train_im,train_classes,train_angles]=digitTrain4DArrayData;
[test_im,test_classes,test_angles]=digitTest4DArrayData;
%% HEP2 data
load(fullfile('databases','hep_proper'));
train_im = X1;
train_classes = Y1;
test_im = X2;
test_classes = Y2;
%% 
augmenter = imageDataAugmenter( ...
    'RandRotation',[0 360], ...
    'RandScale',[0.5 1]);
imageSize = [64 64 1];
augimds = augmentedImageDatastore(imageSize,train_im,train_classes, ...
    'DataAugmentation',augmenter);


%% 2. Select deep learning architecture
% layers = [
%     imageInputLayer([28 28 1]) % Specify input sizes
%     fullyConnectedLayer(10)    % Fully connected is a affine map from 28^2 pixels to 10 numbers
%     softmaxLayer               % Convert to 'probabilities'
%     classificationLayer];      % Specify output layer
%% 
layers = [
    imageInputLayer([64 64 1])  
    convolution2dLayer(3,8,'Padding','same','WeightL2Factor',1)
    batchNormalizationLayer
    reluLayer    
    maxPooling2dLayer(3,'Stride',2) 
    convolution2dLayer(3,16,'Padding','same','WeightL2Factor',1)
    batchNormalizationLayer
    reluLayer    
    maxPooling2dLayer(3,'Stride',2)
    convolution2dLayer(3,32,'Padding','same','WeightL2Factor',1)
    batchNormalizationLayer
    reluLayer    
    fullyConnectedLayer(6)
    softmaxLayer
    classificationLayer];
%% 3. Train deep learning network
miniBatchSize = 100;       
max_epochs = 50;           % Specify how long we should optimize
learning_rate = 0.01;     % Try different learning rates 
options = trainingOptions('adam',...
    'MaxEpochs',max_epochs,...
    'InitialLearnRate',learning_rate, ...
    'Plots', 'training-progress', ...
    L2Regularization=1);   
net = trainNetwork(augimds, layers, options);
%% 4. Test the classifier on the test set
[Y_result2,scores2] = classify(net,test_im);
accuracy2 = sum(Y_result2 == test_classes)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);
%% Test the classifier on the training set
[Y_result1,scores1] = classify(net,train_im);
accuracy1 = sum(Y_result1 == train_classes)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);

