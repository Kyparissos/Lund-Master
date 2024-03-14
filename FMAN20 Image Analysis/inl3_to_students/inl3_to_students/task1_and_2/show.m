clc;
close all;
clearvars;

% Begin by loading the data
load FaceNonFace
nbr_examples = length(Y);
part = cvpartition(nbr_examples, 'HoldOut', 0.20);
    
    % Extract training and test data given the partition above
    X_train = X(:, part.training);
    X_test = X(:, part.test);
    Y_train = Y(:,part.training);
    Y_test = Y(:,part.test);
    nbr_train_examples = length(Y_train);
    nbr_test_examples = length(Y_test);
    classification_data = class_train(X_train, Y_train);

    predictions_test = zeros(1, nbr_test_examples);
    for j = 1 : nbr_test_examples
        % YOU SHOULD IMPLEMENT THE FUNCTION classify!
        predictions_test(j) = classify(X_test(:, j), classification_data);
    end
    i = 3 ;
    im = reshape(X_test(:,i),[19,19]);
    imagesc(im)
    if predictions_test(i)==1
        fprintf('face\n');
    elseif predictions_test(i)==-1
        fprintf('non face\n');
    end
    Y_test(i) ;   
    if predictions_test(i)==Y_test(i)
        fprintf('predicted correctly\n'); 
    else
        fprintf('predicted incorrectly\n');
    end
