% ***** servo_model.m *****
%
% Output: After running the script you will have a 
% model of the process in the variable P.
%
% Calculates the process model according to the description in
% the lab instructions. The process input is the motor
% voltage u, process output is the position of mass 2.

% Process parameters, see lab instructions.
m1=2.3;
m2=2.1;
d1=3.2;
d2=8.6;
k=400;
km=2.95;
ky1=280;
ky2=280;

% System matrices
A=[0 1 0 0;-k/m1 -d1/m1 k/m1 0; 0 0 0 1; k/m2 0 -k/m2 -d2/m2];
B=[0 km/m1 0 0]';
C_proc=[ky1 0 0 0; 0 0 ky2 0];
C1=C_proc(1,:);
C2=C_proc(2,:);

% Process model on transfer function form
P=minreal(zpk(ss(A,B,C2,0)));
P_sim=P;
%% 
[p,z] = pzmap(P);
axis equal
grid on
margin(P);
%% 
%impulse
