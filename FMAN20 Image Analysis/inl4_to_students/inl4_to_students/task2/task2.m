%% Matlab stub for task 2 in assignment 4 in Image analysis
clear;clc;
load heart_data % load data
mean_background = mean(background_values);
mean_chamber = mean(chamber_values);
std_background = std(background_values);
std_chamber = std(chamber_values);
[M,N] = size(im); 
n = M*N; % Number of image pixels
Neighbours = edges4connected(M,N); % use 4-neighbours (or 8-neighbours with edges8connected)
i=Neighbours(:,1);
j=Neighbours(:,2);
A = sparse(i,j,1,n,n); % create sparse matrix of connections between pixels 
% Choose weights:
% Decide how important a short curve length is:
lambda = 0.1;
A =  lambda * A;
% set regularization term so  that A_ij = lambda
mu1 = mean_background;
mu2 = mean_chamber;
Ts = sparse((im(:)-mu1).^2)/(2*std_background); 
Tt = sparse((im(:)-mu2).^2)/(2*std_chamber); 
% create matrix of the full graph, adding source and sink as nodes n+1 and
% n+2 respectively
F = sparse(zeros(n+2,n+2));
F(1:n,1:n) = A; % set regularization weights
F(n+1,1:n) = Ts'; % set data terms 
F(1:n,n+1) = Ts; % set data terms 
F(n+2,1:n) = Tt'; % set data terms 
F(1:n,n+2) = Tt; % set data terms 
% make sure that you understand what the matrix F represents!
Fg = graph(F); % turn F into a graph Fg
help maxflow % see how Matlab's maxflow function works
[MF,GF,CS,CT] = maxflow(Fg,n+1,n+2); % run maxflow on graph with source node (n+1) and sink node (n+2)
disp(MF) % shows the optmization value (maybe not so interesting)
% CS contains the pixels connected to the source node (including the source
% node n+1 as final entry (CT contains the sink nodes).
% We can construct out segmentation mask using these indices
seg = zeros(M,N);
seg(CS(1:end-1)) = 1; % set source pixels to 1
imagesc(im);
imagesc(seg); % show segmentation 







