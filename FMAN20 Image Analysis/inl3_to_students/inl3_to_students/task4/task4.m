% Task 4: Fit lines to data points, using least squares and RANSAC

% Clear up
clc;
close all;
clearvars;

% Begin by loading data points from linedata.mat
load linedata

N = length(xm); % number of data points

% Plot data
plot(xm, ym, '*'); hold on;
xlabel('x') 
ylabel('y')
title('LS & RANSAC') 
x_fine = [min(xm)-0.05,max(xm)+0.05]; % used when plotting the fitted lines

% Fit a line to these data points with least squares
% Here you should write code to obtain the p_ls coefficients (assuming the
% line has the form y = p_ls(1) * x + p_ls(2)).
%p_ls = [rand(), 6]; % REMOVE AND REPLACE WITH LEAST SQUARES SOLUTION
%plot(x_fine, p_ls(1) * x_fine + p_ls(2))
train_xm = [xm,ones(N,1)];
train_ym = ym;
p_ls = train_xm\train_ym;
plot(x_fine, p_ls(1) * x_fine + p_ls(2))

% Fit a line to these data points using RANSAC.
% sampleSize = 2; % number of points to sample per trial
% maxDistance = 2; % max allowable distance for inliers
% points = [xm,ym];
% fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
% evalLineFcn = ...   % distance evaluation function
%   @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);
% 
% [modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn,sampleSize,maxDistance);
% modelInliers = polyfit(points(inlierIdx,1),points(inlierIdx,2),1);
% inlierPts = points(inlierIdx,:);
% x = [min(inlierPts(:,1)) max(inlierPts(:,1))];
% y = modelInliers(1)*x + modelInliers(2);
% plot(x, y, 'g-')    

bestins = 0;
theshold = 1;
for i=1:100
    idx = randperm(N,2);
    k = (ym(idx(1),1)-ym(idx(2),1))/(xm(idx(1),1)-xm(idx(2),1));
    m = ym(idx(1),1)- k*xm(idx(1),1);
    line = [k -1 m];
    distance = abs(line*[xm';ym';ones(1,N)]);
    ins = sum(distance<theshold);
    if ins > bestins
        bestins = ins;
        trueline = line;
    end
end
pa(1) = trueline(1);
pa(2) = trueline(3);
plot(x_fine,pa(1)*x_fine+pa(2),'k--')


%p_ransac = [rand(), 6]; % REMOVE AND REPLACE WITH RANSAC SOLUTION
%plot(x_fine, p_ransac(1) * x_fine + p_ransac(2), 'k--')

% Legend --> show which line corresponds to what (if you need to
% re-position the legend, you can modify rect below)
h=legend('data points', 'least-squares','RANSAC');
rect = [0.20, 0.65, 0.25, 0.25];
set(h, 'Position', rect)

% After having plotted both lines, it's time to compute errors for the
% respective lines. Specifically, for each line (the least squares and the
% RANSAC line), compute the least square error and the total
% least square error. Note that the error is the sum of the individual
% squared errors for each data point! In total you should get 4 errors. Report these
% in your report, and comment on the results. OBS: Recall the distance formula
% between a point and a line from linear algebra, useful when computing orthogonal
% errors!

% WRITE CODE BELOW TO COMPUTE THE 4 ERRORS
% LS_squared_vertical_errors = 0;
% LS_squared_orthogonal_errors = 0;
% RANSAC_squared_vertical_errors = 0;
% RANSAC_squared_orthogonal_errors = 0;
% for i =1:N
%     
%     LS_squared_vertical_errors = (p_ls(1) * xm(i) + p_ls(2) - ym(i))^2 +  LS_squared_vertical_errors;
%     LS_squared_orthogonal_errors = (abs(p_ls(1)*xm(i)+p_ls(2)-ym(i))/sqrt(p_ls(1)^2+1))^2 + LS_squared_orthogonal_errors;
%     RANSAC_squared_vertical_errors = (modelInliers(1) * xm(i) + modelInliers(2) - ym(i))^2 + RANSAC_squared_vertical_errors;
%     RANSAC_squared_orthogonal_errors = (abs(modelInliers(1)*xm(i) + modelInliers(2)-ym(i))/sqrt(modelInliers(1)^2+1))^2 + RANSAC_squared_orthogonal_errors;
% end
% LS_squared_vertical_errors
% LS_squared_orthogonal_errors
% RANSAC_squared_vertical_errors
% RANSAC_squared_orthogonal_errors