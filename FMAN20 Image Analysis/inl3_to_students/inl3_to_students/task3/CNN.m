% Clear up
clc;
close all;
clearvars;

% Begin by loading the data
load FaceNonFace
nbr_examples = length(Y);
nbr_trials = 50;
err_rates_test = zeros(nbr_trials, 1);
err_rates_train = zeros(nbr_trials, 1);
for i = 1 : nbr_trials
    
    % First split data into training / testing (80% train, 20% test)
    part = cvpartition(nbr_examples, 'HoldOut', 0.20);
    
    % Extract training and test data given the partition above
    X_train = X(:, part.training);
    X_test = X(:, part.test);
    Y_train = Y(part.training);
    Y_test = Y(part.test);
    nbr_train_examples = length(Y_train);
    nbr_test_examples = length(Y_test);

    % CNN training
    net = trainSimpleCNN(X_train,Y_train);

    % CNN test
    prediction_Y_test = predictSimpleCNN(net,X_test);
    prediction_Y_train = predictSimpleCNN(net,X_train);

    prediction_diff_test = prediction_Y_test - Y_test;
    prediction_diff_train = prediction_Y_train - Y_train;

    err_rate_test = nnz(prediction_diff_test) / nbr_test_examples;
    err_rate_train = nnz(prediction_diff_train) / nbr_train_examples;
    err_rates_test(i, 1) = err_rate_test;
    err_rates_train(i, 1) = err_rate_train;
end

mean_err_test = mean(err_rates_test);
mean_err_train = mean(err_rates_train);