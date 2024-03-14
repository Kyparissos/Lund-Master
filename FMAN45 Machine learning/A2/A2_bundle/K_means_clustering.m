function [y,C] = K_means_clustering(X,K)

% Calculating cluster centroids and cluster assignments for:
% Input:    X   DxN matrix of input data
%           K   Number of clusters
%
% Output:   y   Nx1 vector of cluster assignments
%           C   DxK matrix of cluster centroids

[D,N] = size(X);

intermax = 50;
conv_tol = 1e-6;
% Initialize
C = repmat(mean(X,2),1,K) + repmat(std(X,[],2),1,K).*randn(D,K);
y = zeros(N,1);
Cold = C;

for kiter = 1:intermax
    % CHANGE
    % Step 1: Assign to clusters
    for i=1:N
        y(i) = step_assign_cluster(X(:, i), C, K);
    end
    
    % Step 2: Assign new clusters
    C = step_compute_mean(X, C, y);
        
    if fcdist(C,Cold) < conv_tol
        return
    end
    Cold = C;
    % DO NOT CHANGE
end

end

function d = fxdist(x,C)
    % CHANGE
    d = norm(x - C);
    % DO NOT CHANGE
end

function d = fcdist(C1,C2)
    % CHANGE
    d = norm(C1 - C2);
    % DO NOT CHANGE
end

function index = step_assign_cluster(x, C, K)
    distances = zeros(K, 1);
    for i=1:K
        distances(i) = fxdist(x, C(:, i));
    end
    index = find(distances == min(distances));
end

function [C_new, C_diff] = step_compute_mean(X, C, y)
    K = size(C, 2);
    C_new = zeros(size(C));
    C_diff = zeros(K, 1);

    for i=1:K
        C_new(:, i) = mean(X(:, i == y), 2);
        C_diff(i) = fcdist(C(:, i), C_new(:, i));
    end
end
