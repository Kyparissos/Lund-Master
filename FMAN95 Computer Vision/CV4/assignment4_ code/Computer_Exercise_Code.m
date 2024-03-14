clc
clear all

%% Computer Exercise01 - Robust Homography Estimation and Stitching
im_a = imread("a.jpg");
im_b = imread("b.jpg");

f1 = detectSIFTFeatures(rgb2gray(im_a));
[d1,valid_points1] = extractFeatures(rgb2gray(im_a), f1);
f2 = detectSIFTFeatures(rgb2gray(im_b));
[d2,valid_points2] = extractFeatures(rgb2gray(im_b), f2);

% Perform the matching
matches = matchFeatures(d1,d2,'MatchThreshold',30,'MaxRatio',0.7);
% Extract the matched keypoints
matchedPoints1 = valid_points1(matches(:,1),:); 
matchedPoints2 = valid_points2(matches(:,2),:);
xA = matchedPoints1.Location'; 
xB = matchedPoints2.Location';
xA = [xA;ones(1,length(matches))];
xB = [xB;ones(1,length(matches))];

% Find a homography describing the transformation
threshold = 5;
best_inliers = [];
corr_points = 4;
best_H = [];
for ii=1:50
    rand = randperm(length(matches), corr_points);
    rand_xA = xA(:, rand);
    rand_xB = xB(:, rand);

    M = [];
    [rows,cols] = size(rand_xA);
%     for i=1:c
%         rand_xA_ = rand_xA(:, i);
%         rand_xB_ = rand_xB(:, i);
%         for j=1:r
%             M(r*(i-1)+j, r*(j-1)+1:r*j) = rand_xA_';
%         end    
%         M(r*(i-1)+1:r*i, 3*r+i) = -rand_xB_;
%     end
    for c=1:cols
        for r=1:rows
            zero__ = zeros(1, 3*2 - length(zeros(1, 3*(r-1))));
            rand_xA_ = rand_xA(:, c);
            rand_xB_ = rand_xB(r, c);
            M = [M;zeros(1, 3*(r-1)) rand_xA_' zero__ zeros(1, c-1) -rand_xB_ zeros(1, cols-c)];
        end
    end 

    [U, S, V] = svd(M);
    v = V(:, end);
    H = [v(1:3,:)'; v(4:6,:)'; v(7:9,:)'];

    trans_xA = pflat(H*xA);
    inliers = sqrt(sum((trans_xA(1:2,:)-xB(1:2,:)).^2))<5;
    n_inliers = sum(inliers == 1);
    if n_inliers > sum(best_inliers)
        best_H = H;
        best_inliers = n_inliers;
    end
end



% Creates a transfomation that matlab can use for images. Note: MATLAB uses transposed transformations.
Htform = projective2d(best_H');
% Sets the size and output bounds of the new image.
Rout = imref2d(size(im_a),[-200 800] ,[ -400 600]);

% Transforms the image
[Atransf] = imwarp(im_a,Htform,'OutputView',Rout); 
% Creates a larger version of the second image
Idtform = projective2d(eye (3));
[Btransf] = imwarp(im_b,Idtform ,'OutputView',Rout);

% Writes both images in the new image. 
% (A somewhat hacky solution is needed since pixels outside the valid image area are not allways zero ...)
AB = Btransf; 
AB(Btransf < Atransf) = Atransf(Btransf < Atransf );
% Plots the new image with the correct axes
imagesc(Rout.XWorldLimits ,Rout.YWorldLimits ,AB);

%% Computer Exercise02 - Robust Essential Matrix Estimation
clear all
A = imread("im1.jpg");
B = imread('im2.jpg');
load compEx2data.mat

x1 = [x{1}; ones(1,size(x{1},2))];
x2 = [x{2}; ones(1,size(x{2},2))];
x1n = K^-1 * x1;
x2n = K^-1 * x2;

bestE = [];
best_in = 0;
best_nbr = 0;
for n = 1:100
    rand = randperm(size(x1n, 2), 5);
    rx1 = x1n(:,rand);
    rx2 = x2n(:,rand);
    E = fivepoint_solver(rx1, rx2);
    for i = 1:size(E,2)
        F = (K^-1)' * E{i} * K^-1;

        l1 = pflat(F' * x2);
        l1 = l1./sqrt(repmat(l1(1,:).^2 + l1(2,:).^2 ,[3 1]));
        l2 = pflat(F * x1);
        l2 = l2./sqrt(repmat(l2(1,:).^2 + l2(2,:).^2 ,[3 1]));

        d_1 = abs(sum(l1.*x1));
        d_2 = abs(sum(l2.*x2));

        inliers = (d_1<5) & (d_2<5);
        nbr_inls = sum(inliers);
        if nbr_inls > best_nbr
            bestE = E{i};
            best_in = inliers;
            best_nbr = nbr_inls;
        end
    end
end

[U2,S2,V2] = svd(bestE);

W = [0 -1 0; 1 0 0; 0 0 1];
P1 = [eye(3,3) zeros(3,1)];
P2 = cell(1,4);
P2{1} = [U2*W*V2' U2(:,3)];
P2{2} = [U2*W*V2' -U2(:,3)];
P2{3} = [U2*W'*V2' U2(:,3)];
P2{4} = [U2*W'*V2' -U2(:,3)];

X = cell(1,length(P2));
for j = 1:length(P2)
    X{j} = DLT(x1n, x2n, P1, P2{j});
end
d1 = depth(P1, P2{1}, X{1});
d2 = depth(P1, P2{2}, X{2});
d3 = depth(P1, P2{3}, X{3});
d4 = depth(P1, P2{4}, X{4});
d = [d1 d2 d3 d4];
[~,best] = max(d);
P1 = K * P1;
P_best = K * P2{best};
X_best = X{best};
x ={x1(:,best_in==1),x2(:,best_in==1)};
X_best=X_best(:,best_in==1);

figure()
plot3(X_best(1,:), X_best(2,:), X_best(3,:), '.', 'Markersize', 2.5)
axis equal
hold on

[err, res] = ComputeReprojectionError({P1, P_best}, X_best, x);
RMS = sqrt(err/size(res,2));
figure()
hist(res,100)
%% Computer Exercise03 - Calibrated Structure from Motion and Local Optimization
P = {P1, P_best};
gammak = 10^-10;
U = X_best;
u = x;

[err, res] = ComputeReprojectionError(P,U,u);
figure();
plot(0, sum(res), 'r*')
hold on

Pnew = P;
Unew = U;
for i = 1:10
    % Computes the r and J matrices for the appoximate linear least squares problem .
    [r,J] = LinearizeReprojErr(Pnew, Unew, u);
    
    % Computes the LM update .
    deltav = -gammak * J' * r;
    
    % Updates the variabels
    [Pnew, Unew] = update_solution(deltav,Pnew,Unew);

    [err, res] = ComputeReprojectionError(Pnew,Unew,u);
    plot(i, sum(res), 'r*')
    hold on
end

RMS = sqrt(err/size(res,2))

%% Computer Exercise04 - Levenberg-Marquardt Method

P = {P1, P_best};
gammak = 10^-10;
lambda = 0.01;

U = X_best;
u = x;

[err, res] = ComputeReprojectionError(P,U,u);
figure();
plot(0, sum(res), 'r*')
hold on

Pnew = P;
Unew = U;
for i = 1:10
    % Computes the r and J matrices for the appoximate linear least squares problem .
    [r,J] = LinearizeReprojErr(Pnew, Unew, u);
    
    % Computes the LM update .
    C = J' * J + lambda * speye(size(J,2));
    c = J' * r;
    deltav = - C \ c ;
    
    % Updates the variabels
    [Pnew, Unew] = update_solution(deltav,Pnew,Unew);

    [err, res] = ComputeReprojectionError(Pnew,Unew,u);
    plot(i, sum(res), 'r*')
    hold on
end

RMS = sqrt(err/size(res,2))