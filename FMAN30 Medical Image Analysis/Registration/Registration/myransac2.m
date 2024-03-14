function [best_R, best_t, best_s,percent] = myransac2(x, y, n)
    % Fit Euclidian transform: n = 3
    % Fit Similarity transform: n = 4
    curr_best = 0;
    best_R = [1,0;0,1];
    best_t = [0;0];
    best_s = 1;
    distance = 5; % Max distance between y_i and the projection of x_i
    rng('shuffle')
    for iter = 1:200
        idxs = randperm(length(x), n);
        [R, t, s] = similarity_registration(x(:, idxs), y(:, idxs));
        new_y = s*R*x + t;
        inlier_idxs = vecnorm(new_y - y) < distance;
        nb_inliers = sum(inlier_idxs, 'all');
        if nb_inliers > curr_best && (s > 0) 
                curr_best = nb_inliers;
                best_s = s;
                best_R = R;
                best_t = t;
        end        
    end
    percent = curr_best/length(x);
%     T = [best_R, best_t; 0, 0, 1];
end
