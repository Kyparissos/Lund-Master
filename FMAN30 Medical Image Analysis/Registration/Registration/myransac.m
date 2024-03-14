% function [best_R, best_t,percent] = myransac(x, y, n)
%     % Fit Euclidian transform: n = 3
%     % Fit Similarity transform: n = 4
%     curr_best = 0;
%     best_R = [1,0;0,1];
%     best_t = [0;0];
%     distance = 10; 
%     rng('shuffle')
%     for iter = 1:1000
%         idxs = randperm(length(x), n);
%         [R, t] = rigid_registration(x(:, idxs), y(:, idxs));
%         new_y = R*x + t;
%         inlier_idxs = vecnorm(new_y - y) < distance;
%         if sum(inlier_idxs, 'all') >= max(0.1*curr_best, n)
%             % Model has potential, fit transformation to all inliers
%             [R, t] = rigid_registration(x(:, inlier_idxs), y(:, inlier_idxs));
%             new_y = R*x + t;
%             nb_inliers = sum(vecnorm(new_y - y) < distance, 'all');
%             if nb_inliers > curr_best
%                 curr_best = nb_inliers;
%                 best_R = R;
%                 best_t = t;
%             end
%         end
%     end
%     percent = curr_best/length(x);
% %     T = [best_R, best_t; 0, 0, 1];
% end
%1
function [best_R, best_t,percent] = myransac(x, y, n)
    curr_best = 0;
    best_R = [1,0;0,1];
    best_t = [0;0];
    distance = 10; % Max distance between y_i and the projection of x_i
    for iter = 1:100
        idxs = randperm(length(x), n);
        [R, t] = rigid_registration(x(:, idxs), y(:, idxs));
        new_y = R*x + t;
        inlier_idxs = vecnorm(new_y - y) < distance;
        nb_inliers = sum(inlier_idxs, 'all');
        if nb_inliers > curr_best
                curr_best = nb_inliers;
                best_R = R;
                best_t = t;
        end        
    end
    percent = curr_best/length(x);
end