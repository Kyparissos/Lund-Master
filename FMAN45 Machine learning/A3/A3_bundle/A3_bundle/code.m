close all
clear all
load('cifar10_baseline.mat')
%% Plot filters
conv_layer = net.layers{1, 2};
conv_weights = conv_layer.params.weights;

batch_size = 16;
for label = 1:batch_size
    subplot(4, 4, label)
    imshow(conv_weights(:,:,label))
end
%% 
%% Find misclassified images
addpath(genpath('./'));
x_test = loadMNISTImages('data/mnist/t10k-images.idx3-ubyte');
z_test = loadMNISTLabels('data/mnist/t10k-labels.idx1-ubyte');
z_test(z_test==0) = 10;
x_test = reshape(x_test, [28, 28, 1, 10000]);

predictions = zeros(numel(z_test), 1);
batch_size = 16;
for label = 1:batch_size:size(z_test)
    data_indices = label : min(label + batch_size - 1, numel(z_test));
    output_by_layer = evaluate(net, x_test(:,:,:,data_indices), z_test(data_indices));
    probabilities = output_by_layer{end-1};
    [~, prediction] = max(probabilities, [], 1);
    predictions(data_indices) = prediction;
end
predictions(predictions == 10) = 0;

%% Plot misclassified images
z_test(z_test==10) = 0;
misclassified = find(predictions ~= z_test, 16);
for label = 1:16
    subplot(4, 4, label)
    imshow(x_test(:, :, :, misclassified(label)));
    text = strcat([...
        'Ground Truth: ',...
        num2str(z_test(misclassified(label))),...
        newline,...
        'Prediction: ',...
        num2str(predictions(misclassified(label)))...
        ]);
    title(text);
end
%% 
confusion_matrix = confusionmat(z_test, predictions);
%% 
figure
cm = confusionchart(z_test, predictions, ...
    'ColumnSummary','column-normalized', ...
    'RowSummary','row-normalized');
%% 
load('cifar10_baseline.mat')
%% Plot filters
conv_layer = net.layers{1, 2};
conv_weights = conv_layer.params.weights;

batch_size = 16;
for label = 1:batch_size
    subplot(4, 4, label)
    imshow(conv_weights(:,:,:,label)*10)
end
%% 
%% Find misclassified images
addpath(genpath('./'));
[X_train, z_train, x_test, z_test, classes] = load_cifar10(1);
data_mean = mean(mean(mean(x_test, 1), 2), 4); 
x_test_norm = bsxfun(@minus, x_test, data_mean);


predictions = zeros(numel(z_test), 1);
batch_size = 16;
for label = 1:batch_size:size(z_test)
    data_indices = label : min(label + batch_size - 1, numel(z_test));
    output_by_layer = evaluate(net, x_test_norm(:,:,:,data_indices), z_test(data_indices));
    probabilities = output_by_layer{end-1};
    [~, prediction] = max(probabilities, [], 1);
    predictions(data_indices) = prediction;
end


%% Plot misclassified images

misclassified = find(predictions ~= z_test, 9);
for label = 1:9
    subplot(3, 3, label)
    imshow(x_test(:, :, :, misclassified(label))/255);
    text = strcat([...
        'Ground Truth: ',...
        classes(z_test(misclassified(label))),...
        newline,...
        'Prediction: ',...
        classes(predictions(misclassified(label)))...
        ]);
    title(text);
end
%% 
confusion_matrix = confusionmat(double(z_test), predictions);
%% 
figure
cm = confusionchart(double(z_test), predictions, ...
    'ColumnSummary','column-normalized', ...
    'RowSummary','row-normalized');
