% This function is to be completed as a preparatory exercise.

% Your task is to finish this function which takes a LinearTankModel
% object, and square matrices Q1, Q2, M, N as inputs. 
% Notation is consitent with section 4 in the lab manual. 

function [CL_lqg, L, Kff] = setup_lqg(lt, Q1, Q2, R1, R2)
% Unpack state matrices
[Phi, Gam, C, D, Ts] = ssdata(lt.DtSys);

% Extract P 
P=lt.DtSys;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       Your code starts here        %%%%%

% Generate an LQG driven by white noise, 
% measuring only h_1 & h_2

%% Preparatory 2
% feedback gain
K = dlqr(Phi, Gam, Q1, Q2);
    
%% Preparatory 3
% Kalman estimator
G= eye(4);
L=dlqe(Phi, G, C, D, R1, R2);

% Generate LQG state-space object
% dreg assumes a negative feedback loop. 

[Ac, Bc, Cc, Dc]=dreg(Phi,Gam,C,D,K,L);
CL_lqg=ss(Ac,Bc,Cc,Dc,Ts);
    
%% Preparatory 4
% Calculate feed-forward gain. It should ensure DC gain = 1.

% Choose one of the Kff below. Motivate your choice

% Kff = eye(2)*evalfr(P,0);
% Kff = eye(2)*evalfr(P,1);
% Kff = [eye(2) zeros(2)]*evalfr(P,0);
% Kff = [eye(2) zeros(2)]*evalfr(P,1);
% Kff = [zeros(2) eye(2)]*evalfr(P,0);
% Kff = [zeros(2) eye(2)]*evalfr(P,1);
% Kff = [eye(2) zeros(2)]*evalfr(P,inf);
% Kff = [zeros(2) eye(2)]*evalfr(P,inf);

% Kff = inv(eye(2)*evalfr(P,0));
% Kff = inv(eye(2)*evalfr(P,1));
% Kff = inv([eye(2) zeros(2)]*evalfr(P,0));
Kff = inv([eye(2) zeros(2)]*evalfr(P,1));
% Kff = inv([zeros(2) eye(2)]*evalfr(P,0));
% Kff = inv([zeros(2) eye(2)]*evalfr(P,1));
% Kff = inv([eye(2) zeros(2)]*evalfr(P,inf));
% Kff = inv([zeros(2) eye(2)]*evalfr(P,inf));

% Kff = eye(2)*evalfr(P,0) + K*L;
% Kff = eye(2)*evalfr(P,1) + K*L;
% Kff = [eye(2) zeros(2)]*evalfr(P,0) + K*L;
% Kff = [eye(2) zeros(2)]*evalfr(P,1) + K*L;
% Kff = [zeros(2) eye(2)]*evalfr(P,0) + K*L;
% Kff = [zeros(2) eye(2)]*evalfr(P,1) + K*L;
% Kff = [eye(2) zeros(2)]*evalfr(P,inf) + K*L;
% Kff = [zeros(2) eye(2)]*evalfr(P,inf) + K*L;

% Kff = inv(eye(2)*evalfr(P,0)) + K*L;
% Kff = inv(eye(2)*evalfr(P,1)) + K*L;
% Kff = inv([eye(2) zeros(2)]*evalfr(P,0)) + K*L;
% Kff = inv([eye(2) zeros(2)]*evalfr(P,1)) + K*L;
% Kff = inv([zeros(2) eye(2)]*evalfr(P,0)) + K*L;
% Kff = inv([zeros(2) eye(2)]*evalfr(P,1)) + K*L;
% Kff = inv([eye(2) zeros(2)]*evalfr(P,inf)) + K*L;
% Kff = inv([zeros(2) eye(2)]*evalfr(P,inf)) + K*L;

% Hint: Instead of evalfr(P,?) you can also use dcgain(P)
%   Which evalfr(P,?) does this correspond to? Why?

%%%%%        Your code ends here         %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end