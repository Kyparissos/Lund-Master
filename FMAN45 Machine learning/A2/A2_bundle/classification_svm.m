function performance = classification_svm(predicted_labels, true_labels)
    image_0 = true_labels(predicted_labels == 0);
    image_1 = true_labels(predicted_labels == 1);

    performance = zeros(2,2);
    performance(1,1) = sum(image_0 == 0);
    performance(1,2) = sum(image_0 == 1);
    performance(2,1) = sum(image_1 == 0);
    performance(2,2) = sum(image_1 == 1);
end