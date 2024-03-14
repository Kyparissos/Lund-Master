%% Matlab stub for task 2 in assignment 4 in Image analysis
load heart_data % load data

% Take mean and std of data sets
mean_background = mean2(background_values);
std_background = std2(background_values);
mean_chamber = mean2(chamber_values);
std_chamber = std2(chamber_values);

[M,N] = size(im); % height of image

n = M*N; % Number of image pixels

% create neighbour structure
Neighbours = edges4connected(M,N); % use 4-neighbours / 8-neighbours 

i=Neighbours(:,1); 
j=Neighbours(:,2);
A = sparse(i,j,1,n,n); % create sparse matrix of connections between pixels

% We can make A into a graph, and show it (test this for example 
% for M = 5, N = 6 to see. For the full image it's not easy to see structure)
% Ag = graph(A);
% plot(Ag);

% Choose weights:
% Decide how important a short curve length is using the pioir of the data
% set
lambda = length(chamber_values)/(length(chamber_values) ...
                               + length(background_values));
L = (- log(lambda));

A =  L * A; % set regularization term so  that A_ij = lambda

% Cauculate likelihood distribution function values
Pd_chamber = normpdf(im,mean_chamber,std_chamber);
Large_chamber = normpdf(mean_chamber,mean_chamber,std_chamber);
Pd_background = normpdf(im,mean_background,std_background);
Large_background = normpdf(mean_background,mean_background,std_background);

% % Tune points manually, uncomment for corresponding method
%% For 4-connected neighbours
% Pd_chamber(34,41) = Large_chamber; 
% Pd_chamber(24:28,60:67) = Large_chamber; 
% Pd_chamber(40:41,70:71) = Large_chamber; 
% Pd_background(48,90) = Large_background;
% Pd_background(63,77) = Large_background;
% Pd_background(1:16,4:8) = Large_background;
% Pd_background(77,54:57) = Large_background;
% Pd_background(87,45:50) = Large_background;
% Pd_background(84:91,45:50) = Large_background;
% Pd_background(85:96,72:96) = Large_background;
%% For 8-connected neighbours
% Pd_chamber(22:27,60:69) = Large_chamber; 
% Pd_background(1:12,4:9) = Large_background;
% Pd_background(84:91,45:50) = Large_background;
% Pd_background(85:96,72:96) = Large_background;

% Use likelihood PDFs as weights for edges in the graph
mu1 = sparse(- log(Pd_chamber));
mu2 = sparse(- log(Pd_background));
% Add offset to data set to ensure no negative weights
mu1 = mu1 + log(Large_chamber);
mu2 = mu2 + log(Large_background);

Ts = reshape(mu1,[],1); % set weights to source, according to assignment!
% Ts = sparse((im(:)-mu1).^2); for standard weighting without statistical weight
Tt = reshape(mu2,[],1); % set weights to sink, according to assignment!
% Tt = sparse((im(:)-mu2).^2); for standard weighting without statistical weight

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

% help maxflow % see how Matlab's maxflow function works

[MF,GF,CS,CT] = maxflow(Fg,n+1,n+2); % run maxflow on graph with source node (n+1) and sink node (n+2)

% disp(MF) % shows the optmization value (maybe not so interesting)

% CS contains the pixels connected to the source node (including the source
% node n+1 as final entry (CT contains the sink nodes).

% We can construct out segmentation mask using these indices
seg = zeros(M,N);
seg(CS(1:end-1)) = 1; % set source pixels to 1
% imagesc(seg); % show segmentation 

montage({im,seg},"size",[1,nan]);









