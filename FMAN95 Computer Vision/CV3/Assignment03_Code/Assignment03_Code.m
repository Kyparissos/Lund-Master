clear all
clc
%% Computer Exercise 01-The Fundamental Matrix
% Loading Data
load('compEx1data.mat');
im_1 = imread('kronan1.JPG');
im_2 = imread('kronan2.JPG'); 

% Compute Normalization Matrices
% N1 = norm_matrix(x{1});
% N2 = norm_matrix(x{2});
N2 = diag([1 1 1]); N1 = N2;

% Normalize image points
x1_n = N1*x{1};
x2_n = N2*x{2};

% Compute M matrix
M = zeros(length(x1_n), 9);

for i = 1:length(x1_n)
    x_ = x2_n(:,i)*x1_n(:,i)';
    M(i,:) = x_(:)';
end

% SVD
[U,S,V] = svd(M);
v = V(:, end);
norm_F = reshape(v,[3 3]);
% Normalized F from v
[U_,S_,V_] = svd(norm_F);
S_(3,3) = 0;
norm_F = U_*S_*V_';
figure(1)
plot(diag(x2_n'*norm_F*x1_n));
% Unnormalized F
E_unn = N2'*norm_F*N1;
E_unn = E_unn./E_unn(3,3);

% Computes the Epipolar Lines
l = E_unn*x{1};
l = l./sqrt(repmat(l(1,:).^2 + l(2,:).^2,[3 1]));

% Randomly Pick Up 20 points and their corresponding lines
rand = randi(length(x{2}),1,20);
x_2 = x{2}(:,rand);
l_ = l(:,rand);
figure(2)
imagesc(im_2)
hold on
plot(x_2(1,:),x_2(2,:),'r*')
rital(l_)

figure(3);
hist(abs(sum(l.*x{2})),100);
d = mean(abs(sum(l.*x{2})));
disp(d)
%% Computer Exercise 02-Compute Camera Matrices
P1 = [diag([1 1 1]) zeros(3,1)];
e_2 = null(E_unn');
e2x = [0 -e_2(3) e_2(2); e_2(3) 0 -e_2(1); -e_2(2) e_2(1) 0];
P2 = [e2x*E_unn e_2];

% Use traingulation (with DLT) to compute the 3D points
zero = zeros(3,1);
X = [];
for i = 1:length(x{1})
    M = [P1 -x1_n(:,i) zero; P2 zero -x2_n(:,i)];
    [U, S, V] = svd(M);
    v = V(:, end);
    X(:, i) = v(1:4, 1);
end
X = pflat(X);
% Projection with P1
x1 = pflat(P1*X);
x2 = pflat(P2*X);

% Plot the image points
figure(4);
subplot(1,2,1);
imagesc(im_1);
hold on;
plot(x{1}(1, :), x{1}(2, :), 'b*', 'MarkerSize',10);
hold on;
plot(x1(1, :), x1(2, :), 'r.', 'MarkerSize',10);
subplot(1,2,2);
imagesc(im_2);
hold on;
plot(x{2}(1, :), x{2}(2, :), 'b*', 'MarkerSize',10);
plot(x2(1, :), x2(2, :), 'r.', 'MarkerSize',10);
% axis equal

% Plot the 3D points
figure(5);
plot3(X(1, :), X(2, :), X(3, :), '.');
% axis equal
%% Computer Exercise 03-Normalize the image points using the inverse of K
load compEx1data.mat
load compEx3data.mat
im_1 = imread('kronan1.JPG');
im_2 = imread('kronan2.JPG');
% Normalize the image points using the inverse of K
x1_k = inv(K)*x{1};
x2_k = inv(K)*x{2};

% Set up the matrix M
M = [];
for i = 1:length(x{1})
    xx = x2_k(:,i)*x1_k(:,i)';
   M(i,:) = reshape(xx,[],1);
end

[U, S, V] = svd(M);
v = V(:, end);
Eapprox = reshape(v,[3 3]);
[U_,S_,V_] = svd(Eapprox);
if det(U_*V_') > 0 
    E = U_*diag ([1 1 0])*V_';
else
    V_ = -V_; 
    E = U_*diag ([1 1 0])*V_';
end
E = E./E(3,3);

E_unn = inv(K)'*E*inv(K);
E_unn = E_unn./E_unn(3,3);
l = E_unn*x{1};
l = l./sqrt(repmat(l(1,:).^2 +l(2,:).^2,[3 1]));

rand = randi(length(x{2}),1,20);
x_2 = x{2}(:,rand);
l_ = l(:,rand);
figure(6)
imagesc(im_2)
hold on
plot(x_2(1,:),x_2(2,:),'r*')
rital(l_)

figure(7);
hist(abs(sum(l.*x{2})),100);
%% Computer Exercise 04-Compute Four Camera normSolutions

P1 = [eye(3) zeros(3, 1)];

W = [0 -1 0; 1 0 0; 0 0 1];

P2 = {[U_*W*V_' U_(:, 3)], [U_*W*V_' -U_(:, 3)], [U_*W'*V_' U_(:, 3)], [U_*W'*V_' -U_(:, 3)]};

zero = zeros(3,1);
M_ = [];
X_ = {};
for i = 1:4
    for j = 1:length(x{1})
        M_ = [P1 -x1_k(:, j) zero; P2{i} zero -x2_k(:, j)];
        [U, S, V] = svd(M_);
        v = V(:, end);
        X_{i}(:, j) = v(1:4, 1);
    end
end
m_p = 0;
for i = 1:4
    X_{i} = pflat(X_{i});
    X1_p = P1*X_{i};
    X2_p = P2{i}*X_{i};
    sum_p = sum(X1_p(3,:) > 0) + sum(X2_p(3,:) >0);
    % Most points in front of the camera
    m_p = max(m_p,sum_p);
    if m_p == sum_p
        index = i;
    end
end

% Plot
projection1 = pflat(K*P1*X_{index});
projection2 = pflat(K*P2{index}*X_{index});
figure(8);
subplot(1,2,1);
imagesc(im_1)
hold on
plot(x{1}(1, :), x{1}(2, :), 'b*');
plot(projection1(1,:),projection1(2,:),'ro');
subplot(1,2,2);
imagesc(im_2)
hold on
plot(x{2}(1, :), x{2}(2, :), 'b*');
plot(projection2(1,:),projection2(2,:),'ro');

figure(9);
plot3(X_{index}(1,:),X_{index}(2,:),X_{index}(3,:),'.')
hold on
axis equal
plotcams({P1,P2{index}});
