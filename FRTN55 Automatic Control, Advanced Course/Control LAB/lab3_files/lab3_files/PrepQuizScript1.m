%% Script to use as a starting point in the Lab 3 Preparatory quiz

% Specify tank properties
lt = LinearTankModel;               % Create linearTankModel object
lt.Gammas = [0.3 0.3];              % Set gamma_1 and gamma_2
% Set the sample time
lt.SampleTime = 1;  %Ts

% Needed for the simulation
Ts = lt.SampleTime;

ContSys=lt.CtSys;
t_end=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      Your code starts here     %%%%%

t_0=5;

t= 0:lt.SampleTime:t_end;

u1  = 2*heaviside(t - t_0);
u2  = -2*heaviside(t - t_0);

u = [u1; u2];
lsim(ContSys,u,t)

