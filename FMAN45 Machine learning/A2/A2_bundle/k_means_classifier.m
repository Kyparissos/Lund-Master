function cluster_labels = k_means_classifier(X_labels, y, K)
    cluster_labels = zeros(K, 1);
    for i=1:K
        n_zeros = sum(X_labels(y == i) == 0);
        n_ones = sum(X_labels(y == i) == 1);
        if n_zeros < n_ones
            cluster_labels(i) = 1;
        else
            cluster_labels(i) = 0;
        end
    end
end