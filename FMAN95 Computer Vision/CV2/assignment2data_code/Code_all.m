clear all
clc
%% Computer Exercise01-Calibrated vs. Uncalibrated Reconstruction
load('compEx1data.mat');

T_1 = [1 0 0 0;0 4 0 0;0 0 1 0;0.1 0.1 0 1];
T_2 = [1 0 0 0;0 1 0 0;0 0 1 0;1/16 1/16 0 1];

figure(1);
plot3(X(1,:),X(2,:),X(3,:),'.','Markersize',2,'Color','b');
axis equal
hold on
plotcams(P);

im = imread(imfiles{1});
xproj = pflat(P{1}*X);
% Determines which of the points are visible in image i
visible = isfinite(x{1}(1 ,:));
figure(2);
imagesc(im);
hold on
% Plots a *  at each point coordinate
plot(x{1}(1, visible), x{1}(2, visible),'*');
% Plots a red  o at each visible point in xproj
plot(xproj(1,visible), xproj(2,visible),'ro');
axis equal

% P' = PT^-1, X' = TX
P_1 = {};
P_2 = {};
for i=1:9
    P_1{i} = P{i}*inv(T_1);
    P_2{i} = P{i}*inv(T_2);
end
X_1 = pflat(T_1*X); %new projective solution
X_2 = pflat(T_2*X);

figure(3);
subplot(1,2,1);
plot3(X_1(1,:),X_1(2,:),X_1(3,:),'.','Markersize',2,'Color','b');
axis equal
hold on
plotcams(P_1);
subplot(1,2,2);
plot3(X_2(1,:),X_2(2,:),X_2(3,:),'.','Markersize',2,'Color','b');
axis equal
hold on
plotcams(P_2);


xproj_1 = pflat(P_1{1}*X_1);
xproj_2 = pflat(P_2{1}*X_2);

figure(4);
subplot(1,2,1);
imagesc(im);
hold on
% Plots a *  at each point coordinate
plot(x{1}(1, visible), x{1}(2, visible),'*');
% Plots a red  o’ at each visible point in xproj
plot(xproj_1(1,visible), xproj_1(2,visible),'ro');
axis equal;
subplot(1,2,2);
imagesc(im);
hold on
% Plots a 1*  at each point coordinate
plot(x{1}(1, visible), x{1}(2, visible),'*');
% Plots a red  o at each visible point in xproj
plot(xproj_2(1,visible), xproj_2(2,visible),'ro');
axis equal;

%% Computer Exercise02 - RQ Factorization and Computation of K

[K_1 R_1] = rq(P_1{2}(:,1:3));
K_1 = K_1./K_1(3, 3)
[K_2 R_2] = rq(P_2{2}(:,1:3));
K_2 = K_2./K_2(3, 3)

%% Computer Exercise03 - Direct Linear Transformation DLT
im_1 = imread("cube1.JPG");
im_2 = imread("cube2.JPG");

load("compEx3data.mat");
%creating the transformation N
N1 = norm_matrix(x{1});
N2 = norm_matrix(x{2});
x1_n = N1*x{1};
x2_n = N2*x{2};

figure(5);
subplot(2,2,1);
plot(x{1}(1,:), x{1}(2,:),'*');
subplot(2,2,2);
plot(x1_n(1,:), x1_n(2,:),'*');
subplot(2,2,3);
plot(x{2}(1,:), x{2}(2,:),'*');
subplot(2,2,4);
plot(x2_n(1,:), x2_n(2,:),'*');
axis equal;

% Compute matrix M
M_1 = zeros(3*length(Xmodel),3*4+length(Xmodel));
M_2 = zeros(3*length(Xmodel),3*4+length(Xmodel));
zero = zeros(1,4);
Xmodel = [Xmodel;ones(1,length(Xmodel))];
for i=1:length(Xmodel)
    M = [Xmodel(:,i)' zero zero;zero Xmodel(:,i)' zero; zero zero Xmodel(:,i)'];
    for j = 1:3
        M_1(3*i-(3-j),1:12) = M(j,:);
        M_2(3*i-(3-j),1:12) = M(j,:);
    end
    M_1(3*i-2:3*i,12+i) = -x1_n(:,i);
    M_2(3*i-2:3*i,12+i) = -x2_n(:,i);
end

% Computes the singular value decomposition of M
[U1,S1,V1] = svd(M_1);
[U2,S2,V2] = svd(M_2);
v_1 = V1(:,end);
v_2 = V2(:,end);

eig1 = S1'*S1; min_eig1 = eig1(end,end);
eig2 = S2'*S2; min_eig2 = eig2(end,end);
% min_eigvalue = min(eig(eigen_matrix>0));
min_mv1 = norm(M_1*v_1);
min_mv2 = norm(M_2*v_2);

% Extract camera matrix P
P{1} = inv(N1)*reshape(v_1(1:12),[4 3])';
P{2} = inv(N2)*reshape(v_2(1:12),[4 3])';

[R1, Q1] = rq(P{1});
[R2, Q2] = rq(P{2});

K1 = R1./R1(3, 3);
K2 = R2./R2(3, 3);

X1 = pflat(P{1}*Xmodel);
X2 = pflat(P{2}*Xmodel);

figure(6);
subplot(1,2,1);
imagesc(im_1);
hold on
plot(x{1}(1,:), x{1}(2,:),'*');
plot(X1(1,:), X1(2,:),'o');
axis equal;
subplot(1,2,2);
imagesc(im_2);
hold on
plot(x{2}(1,:), x{2}(2,:),'*');
plot(X2(1,:), X2(2,:),'o');
axis equal;

figure(7);
plot3(Xmodel(1,:),Xmodel(2,:),Xmodel(3,:),'.','Markersize',2,'Color','b');
hold on
plot3([Xmodel(1,startind );  Xmodel(1,endind )],[Xmodel(2,startind );  Xmodel(2,endind )],[Xmodel(3,startind );  Xmodel(3,endind)],'b-');
plotcams(P)
axis equal;
%% Computer Exercise04 - Feature Extraction and Matching using SIFT
im1 = imread("cube1.JPG");
im2 = imread("cube2.JPG");
[f1 d1] = vl_sift(single(rgb2gray(im1)), 'PeakThresh', 1);
[f2 d2] = vl_sift(single(rgb2gray(im2)), 'PeakThresh', 1);
vl_plotframe(f1);
[matches,scores] = vl_ubcmatch(d1,d2);
x1 = [f1(1, matches(1,:)); f1(2 ,matches(1 ,:))];
x2 = [f2(1, matches(2,:)); f2(2 ,matches(2 ,:))];
perm = randperm(size(matches ,2));
figure(8);
imagesc ([ im1 im2 ]);
hold on
plot([x1(1,perm (1:10)); x2(1,perm (1:10))+size(im1 ,2)] , ...
[x1(2, perm(1:10)); x2(2, perm(1:10))], '-');
%% 
im1 = imread("cube1.JPG");
im2 = imread("cube2.JPG");
% Detect keypoints and compute the descriptors
 f1 = detectSIFTFeatures(rgb2gray(im1));
 [d1, valid_points1] = extractFeatures(rgb2gray(im1), f1);
 f2 = detectSIFTFeatures(rgb2gray(im2));
 [d2, valid_points2] = extractFeatures(rgb2gray(im2), f2);
 % Perform the matching
 matches = matchFeatures(d1,d2,'MatchThreshold',30, 'MaxRatio',0.7);
 % Extract the matched keypoints
 matchedPoints1 = valid_points1(matches(:,1),:);
 matchedPoints2 = valid_points2(matches(:,2),:);
 x41 = matchedPoints1.Location';
 x42 = matchedPoints2.Location';

 perm = randperm(size(matches,1));
 figure;
 imagesc([im1 im2]);
 hold on;
 plot([x41(1,perm(1:10)); x42(1,perm(1:10))+size(im1,2)], ...
 [x41(2,perm(1:10)); x42(2,perm(1:10))],'-');
 hold off;
%% Computer Exercise05 - Triangulation using DLT

% DLT Equations for Triangulation
n = length(x41);
X = [];
P1 = inv(K1) * P{1};
P2 = inv(K2) * P{2};
x_41 = pflat(inv(K1) * [x41; ones(1,n)]);
x_42 = pflat(inv(K2) * [x42; ones(1,n)]);

for i = 1:n
    M = [P1 -x_41(:,i) zeros(3,1); P2 zeros(3,1) -x_42(:,i)];
    [U,S,V] = svd(M);
    v = V(:, end);
    X(:, i) = v(1:4);
end
X = pflat(X);

xp1 = pflat(P{1}*X);
xp2 = pflat(P{2}*X);

% Plot the projections
figure(9);
subplot(1,2,1);
imagesc(im_1);
hold on
plot(x41(1,:), x41(2,:), '*');
plot(xp1(1,:), xp1(2,:), 'ro');
subplot(1,2,2);
imagesc(im_2);
hold on
plot(x42(1,:), x42(2,:), '*');
plot(xp2(1,:), xp2(2,:), 'ro');
%% Without normalization
n = size(x41, 2);
X = [];
x_41 = pflat([x41; ones(1,n)]);
x_42 = pflat([x42; ones(1,n)]);
P1 = P{1};
P2 = P{2};

for i=1:n
    M = [P1 -x_41(:,i) zeros(3,1); P2 zeros(3,1) -x_42(:,i)];
    [U,S,V] = svd(M);
    v = V(:,end);
    X(:, i) = v(1:4);
end
X = pflat(X);

x1_proj = pflat(P1*X);
x2_proj = pflat(P2*X);

figure(10);
subplot(1,2,1)
imagesc(im_1)
hold on
plot(x41(1,:), x41(2,:), '*')
plot(x1_proj(1,:), x1_proj(2,:), 'ro')

subplot(1,2,2)
imagesc(im_2)
hold on
plot(x42(1,:), x42(2,:), '*')
plot(x2_proj(1,:), x2_proj(2,:), 'ro')

%% 
% Finds the points with reprojection error less than 3 pixels in both images
good_points = (sqrt(sum((x41-xp1 (1:2 ,:)).^2)) < 3 & sqrt(sum((x42-xp2 (1:2 ,:)).^2)) < 3);
X = pflat(X(:, good_points));

figure(11);
plot3(X(1, :), X(2, :), X(3, :), 'r.');
hold on;
plot3([Xmodel(1,startind); Xmodel(1,endind)], [Xmodel(2,startind); Xmodel(2,endind)], [Xmodel(3,startind); Xmodel(3,endind)],'b-');
plotcams(P);
