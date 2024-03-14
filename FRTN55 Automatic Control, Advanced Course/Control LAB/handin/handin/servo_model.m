% ***** servo_model.m *****
%
% Defines the linear process model for the flexible servo. Note that the
% km motor parameter is different from the one used in Lab 1.
% 
% Process parameters
m1 = 2.3;
m2 = 2.1;
d1 = 3.2;
d2 = 8.6;
k = 400;
km = 0.443;  % NOTE: this parameter is different compared to Lab 1
ky = 280;

% System matrices
A = [0 1 0 0; -k/m1 -d1/m1 k/m1 0; 0 0 0 1; k/m2 0 -k/m2 -d2/m2];
B = [0; km/m1; 0; 0];
C = [0 0 ky 0];
D = 0;

% Process model in transfer function form
P = zpk(ss(A,B,C,D));
