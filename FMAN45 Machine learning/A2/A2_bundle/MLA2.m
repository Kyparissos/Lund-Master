%% E1
load A2_data.mat

X = train_data_01;

% Map to zero mean images
X = X - mean(X, 2);

[U, ~, ~] = svd(X);
principal_components = U(:, 1:2);
projected = principal_components' * X;

%% Plot
gscatter(projected(1, :), projected(2, :), train_labels_01,'br','ox');
xlabel('First Principal Component');
ylabel('Second Principal Component');
title 'PCA'

%% E2 2 Clusters
K = 2;
[y_2, C_2] = K_means_clustering(X, K);

%% Plot
figure
gscatter(projected(1, :), projected(2, :), y_2, 'br', 'ox');
xlabel('First Principal Component');
ylabel('Second Principal Component');
title 'K Means Clustering (K = 2)'

%% 5 Clusters
K = 5;
[y_5, C_5] = K_means_clustering(X, K);

%% Plot
figure
gscatter(projected(1, :), projected(2, :), y_5, 'brgmc', 'ox*^p');
xlabel('First Principal Component');
ylabel('Second Principal Component');
title 'K Means Clustering (K = 5)'

%% E3 
% 2 Clusters with Centroids
img_C_2a = reshape(C_2(:, 1), [28 28]);
img_C_2b = reshape(C_2(:, 2), [28 28]);

figure
hold on
subplot(1,2,1);
imshow(img_C_2a);
title('Cluster 1')
subplot(1,2,2);
imshow(img_C_2b);
title('Cluster 2')
hold off

%% 5 Clusters with Centroids
img_C_5a = reshape(C_5(:, 1), [28 28]);
img_C_5b = reshape(C_5(:, 2), [28 28]);
img_C_5c = reshape(C_5(:, 3), [28 28]);
img_C_5d = reshape(C_5(:, 4), [28 28]);
img_C_5e = reshape(C_5(:, 5), [28 28]);

figure
hold on
subplot(1,5,1);
imshow(img_C_5a);
title('Cluster 1')
subplot(1,5,2);
imshow(img_C_5b);
title('Cluster 2')
subplot(1,5,3);
imshow(img_C_5c);
title('Cluster 3')
subplot(1,5,4);
imshow(img_C_5d);
title('Cluster 4')
subplot(1,5,5);
imshow(img_C_5e);
title('Cluster 5')
hold off

%% E4
clear;clc
load A2_data

K = 12;

X_train = train_data_01;
X_test = test_data_01;
labels_train = train_labels_01;
labels_test = test_labels_01;

[y_train, C_train] = K_means_clustering(X_train, K);
cluster_labels_train = k_means_classifier(labels_train, y_train, K);

[y_test, C_test] = K_means_clustering(X_test, K);
cluster_labels_test = k_means_classifier(labels_test, y_test, K);

for i=1:K
    n_zeros_train(i) = sum(labels_train(y_train == i) == 0);
    n_ones_train(i) = sum(labels_train(y_train == i) == 1);
    n_zeros_test(i) = sum(labels_test(y_test == i) == 0);
    n_ones_test(i) = sum(labels_test(y_test == i) == 1);

    if cluster_labels_train(i) == 1
        misclassified_train(i) = n_zeros_train(i);
    else
        misclassified_train(i) = n_ones_train(i);
    end

    if cluster_labels_test(i) == 1
        misclassified_test(i) = n_zeros_test(i);
    else
        misclassified_test(i) = n_ones_test(i);
    end

end

misclassified_sum_train = sum(misclassified_train);
misclassification_rate_train = misclassified_sum_train ./ size(X_train, 2);

misclassified_sum_test = sum(misclassified_test);
misclassification_rate_test = misclassified_sum_test ./ size(X_test, 2);

%% 
clear;clc
load A2_data.mat

svm = fitcsvm(train_data_01', train_labels_01);

prediction_train = predict(svm, train_data_01');
prediction_test = predict(svm, test_data_01');

performance_train = classification_svm(prediction_train, train_labels_01');
performance_test = classification_svm(prediction_test, test_labels_01');

%% 

clear;
load A2_data.mat

beta = 6;
svm = fitcsvm(train_data_01', train_labels_01, 'KernelFunction','gaussian', 'KernelScale', beta);

prediction_train = predict(svm, train_data_01');
prediction_test = predict(svm, test_data_01');

performance_train = classification_svm(prediction_train, train_labels_01');
performance_test = classification_svm(prediction_test, test_labels_01');
