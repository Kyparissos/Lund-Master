%% Load data (training data X1 - labels Y1)
% test data X2 - labels Y2
close all
clear all
load(fullfile('databases','hep_proper_mask'));
X1_masks = Y1;
X2_masks = Y2;
load(fullfile('databases','hep_proper'));
%% 
% c = cvpartition(X1,"KFold",4);

%% Use hand made features

nr_of_training_images = size(X1,4);
for i = 1:nr_of_training_images
   [fv,str]=get_features(X1(:,:,1,i),X1_masks(:,:,1,i));
   X1f(i,:)=fv;
end
% X1f = normalize(X1f,1);
nr_of_test_images = size(X2,4);
for i = 1:nr_of_test_images
   [fv,str]=get_features(X2(:,:,1,i),X2_masks(:,:,1,i));
   X2f(i,:)=fv;
end
% X2f = normalize(X2f,1);
%% Visualize the data
Y = tsne(X1f(:,1:6));
figure(1);
gscatter(Y(:,1), Y(:,2),Y1,'rgbcmk','oooooo');
%% 
figure(1);
gscatter(X1f(:,5), X1f(:,6),Y1,'rgbcmk','oooooo');
% xlabel(str{1});
% ylabel(str{2});
%% Train machine learning models to data

%% Fit a decision tree model
disp(' ');
disp('Decision tree');
rng(1) 
model1 = fitctree(X1f(:,:),Y1,'PredictorNames',str, ...    
    'MaxNumCategories',6, ...
    'OptimizeHyperparameters','auto');
% Test the classifier on the training set
[Y_result1,node_1] = predict(model1,X1f);
accuracy1 = sum(Y_result1 == Y1)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);
% Test the classifier on the test set
[Y_result2,node_2] = predict(model1,X2f);
accuracy2 = sum(Y_result2 == Y2)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);

%% Fit a random forest model
disp(' ');
disp('Random forest');
model2 = TreeBagger(50,X1f,Y1,'OOBPrediction','On',...
    'Method','classification');
% Test the classifier on the training set
[Y_result1,node_1] = predict(model2,X1f);
accuracy1 = sum(Y_result1 == Y1)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);
% Test the classifier on the test set
[Y_result2,node_2] = predict(model2,X2f);
accuracy2 = sum(Y_result2 == Y2)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);

%% Fit a support vector machine model
disp(' ');
disp('Suppport vector machine');
rng(1)
model3 = fitcecoc(X1f,Y1,'OptimizeHyperparameters','auto');
% Test the classifier on the training set
[Y_result1,node_1] = predict(model3,X1f);
accuracy1 = sum(Y_result1 == Y1)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);
% Test the classifier on the test set
[Y_result2,node_2] = predict(model3,X2f);
accuracy2 = sum(Y_result2 == Y2)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);

%% Fit a k-nearest neighbour model
disp(' ');
disp('k-nearest neighbour');
rng(1) 
model4 = fitcknn(X1f,Y1,'OptimizeHyperparameters','auto');
% Test the classifier on the training set
[Y_result1,node_1] = predict(model4,X1f);
accuracy1 = sum(Y_result1 == Y1)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);
% Test the classifier on the test set
[Y_result2,node_2] = predict(model4,X2f);
accuracy2 = sum(Y_result2 == Y2)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);

%% Fit a neural network model
disp(' ');
disp('neural network');
rng(1) 
model5 = fitcnet(X1f,Y1,'OptimizeHyperparameters','auto');
% Test the classifier on the training set
[Y_result1,node_1] = predict(model5,X1f);
accuracy1 = sum(Y_result1 == Y1)/numel(Y_result1);
disp(['The accuracy on the training set: ' num2str(accuracy1)]);
% Test the classifier on the test set
[Y_result2,node_2] = predict(model5,X2f);
accuracy2 = sum(Y_result2 == Y2)/numel(Y_result2);
disp(['The accuracy on the test set: ' num2str(accuracy2)]);






