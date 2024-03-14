clc
clear
close all
load dmsa_images % manual segmentation 1st 20 imgs, landmarks patient 1
load models  % coordinates of the landmarks 
models(end+1,:) = models(1,:);
%% 1. Align each shape to the first
% real and imag part
% the manual segmentations
aligned_mean_shape = [real(models(:,1)), imag(models(:,1))]'; 

for iter = 1:5
        aligned_shapes{1} = aligned_mean_shape;
        % Align each shape to the first
        for i = 1:length(models)
            shape_to_align = [real(models(:, i)), imag(models(:, i))]';
            [Rs{i}, ts{i}, ss{i}] = find_similarity_transform(shape_to_align, aligned_mean_shape);
            aligned_shapes{i} = ss{i}*Rs{i}*shape_to_align  + ts{i};
        end
        
        % Calculate the mean of the transformed shapes (by calculating the mean valuefor each point)
        % only use the first 14 coordinates
        % linear interpolation
        mean_coords = zeros(2,15);
        for i = 1:length(aligned_shapes)
            mean_coords = mean_coords + aligned_shapes{i};
        end
         mean_transformed_shape = mean_coords/length(aligned_shapes);

        % Align the mean shape to the first (to guarantee convergence)
        [R, t, s] = find_similarity_transform(mean_transformed_shape, aligned_mean_shape);
        aligned_mean_shape = s*R*mean_transformed_shape  + t;
        
end
% After aligning the shapes calculate the covariance function of the points and the
% corresponding eigen-values and eigen-vectors in order to obtain a shape model
model_post_transform = zeros(size(models));
for i=1:length(ss)
    model_post_transform(:, i) = apply_transformation(Rs{i}, ts{i}, ss{i}, models(:, i));
end
 

figure

plot(mean_transformed_shape(1,:), mean_transformed_shape(2,:),'r-', 'LineWidth', 1)
hold on
for i=1:length(ss)
    complex_coords= model_post_transform(:, i);
    x = real(complex_coords);
    y = imag(complex_coords);
    p = plot(x,y);
    p.Color(4) = 0.2;
end
legend({'Mean shape'})
title('Aligned Shapes')

C = covariance(mean_transformed_shape, model_post_transform);
[evecs, evals] = eig(C);


% Make decscending order
evals = flip(diag(evals),1);

P = flip(evecs, 2);

% Plot eigenvalues
figure
scatter(1:30, evals)
title('Eigenvalues descendingly ordered')

% Plot cumulative energy for eigenvalues
% figure
% plot(0:28, [0;cumsum(evals)/sum(evals)], 'LineWidth', 1)
% hold on

% plot(0:28, 0.95*ones(1,29), 'r-.')
% legend({'Cumulative Energy', '95%'})
% s1 = scatter(0:28, [0; cumsum(evals)/sum(evals)], 50, 'r', 'filled');
% set(get(get(s1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% title('Cumulative Energy')

% Plot mode with different b's
for idx = 1:10
    visualise(mean_transformed_shape, idx, P, evals)
end

%% functions

function [R, t, s] = find_similarity_transform(x, y)
    % Algorithm to find the Similarity transformation from x to y
    avg_x = mean(x,2);
    avg_y = mean(y,2);
    x_tilde = x - avg_x;
    y_tilde = y - avg_y;
    [U, ~, V] = svd(y_tilde*x_tilde');
    R = U*diag([1 det(U*V')])*V';
    s = sum(diag(y_tilde'*R*x_tilde))/sum(vecnorm(x_tilde).^2);
    t = avg_y - s*R*avg_x;
end

function new_coords = apply_transformation(R, t, s, coords)
    coords = [real(coords), imag(coords)]';
    new_coords = s*R*coords + t;
    new_coords = complex(new_coords(1, :), new_coords(2, :));
end

function S = covariance(mean_shape, shapes)

shapes_xs = real(shapes);
shapes_ys = imag(shapes);

[N, M] = size(shapes);

S = zeros(2 * N, 2 * N);

X_bar = [mean_shape(1, :)'; mean_shape(2, :)'];

for i=1:M
    X = [shapes_xs(:, i); shapes_ys(:, i)];
    dX = X - X_bar;
    S = S + dX * dX'; 
end

S = S / M;

end

function [] = visualise(mean_shape, mode_idx, P, evals)
    eval = evals(mode_idx);
    mid = length(P(:, 1)) / 2;
    evec = [P(1:mid, mode_idx), P(mid+1:end, mode_idx)]';
        
    p1 = mean_shape - 2 * sqrt(eval) * evec;
    p2 = mean_shape - sqrt(eval) * evec;
    p3 = mean_shape;
    p4 = mean_shape + sqrt(eval) * evec;
    p5 = mean_shape + 2 * sqrt(eval) * evec;
    
    
    figure
    plot(p1(1, :), p1(2, :), 'r')
    hold on
    plot(p2(1, :), p2(2, :), 'g')
    plot(p3(1, :), p3(2, :), 'b--', 'LineWidth', 1)
    plot(p4(1, :), p4(2, :), 'm')
    plot(p5(1, :), p5(2, :), 'k')
    
    hold off
    legend({'Mean - 2*sqrt(lambda)', 'Mean - sqrt(lambda)', 'Mean', 'Mean + sqrt(lambda)', 'Mean + 2*sqrt(lambda)'})
    title(['Principal mode ' num2str(mode_idx)])
end



