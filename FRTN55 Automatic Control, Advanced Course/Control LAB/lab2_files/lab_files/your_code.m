%% Initialization

% Add paths to RASPlib
addpath([fileparts(mfilename('fullpath')) '/util'])

% Time-step
Ts = 0.01;
h = Ts;
%% Assignment 4

% Alphadot
figure;
plot(alphadot_gyro)                     % plot
y1 = detrend(alphadot_gyro,'linear');   % remove linear trend
figure;
plot(y1)                                % plot again
y1_var = var(y1);                       % calculate stationary variance
figure;
pwelch(y1)                              % plot periodogram (estimate of spectrum)

% Alpha
figure;
plot(alpha_accel)                       % plot
y2 = detrend(alpha_accel,'linear');     % remove linear trend
figure;
plot(y2)                                % plot again
y2_var = var(y2);                       % calculate stationary variance
figure;
pwelch(y2)                              % plot periodogram (estimate of spectrum)

%% Assignment 5

Phi_imu = [1,0;h,1];      % Assignment 5
Gam_imu = 0;      % Assignment 5
Cd_imu  = [1,0;0,1];      % Assignment 5
Dd_imu  = 0;      % Assignment 5
Gd_imu  = [h;0];       % Assignment 5

R1_imu  = 0.001;      % Assignment 5
% R1_imu  = 100;      % Assignment 5

R2_imu  = eye(2);      % Assignment 5

[L_imu, ~, ~, Eig_imu] = dlqe(Phi_imu, Gd_imu, Cd_imu, R1_imu, R2_imu);

% Assignment 5: Analyse L_imu and Eig_imu for different values of R1_imu

%% Assignment 6

y1_var  = 2e-6;
y2_var  = 1e-5; 

Phi_imu = [1,0;h,1];      % Assignment 5
Gam_imu = 0;      % Assignment 5
Cd_imu  = [1,0;0,1];      % Assignment 5
Dd_imu  = 0;      % Assignment 5
Gd_imu  = [h;0];       % Assignment 5

R1_imu  = 0.005;      % Assignment 6
R2_imu  = [y1_var,0;0,y2_var];     % Assignment 6

[L_imu, ~, ~, Eig_imu] = dlqe(Phi_imu, Gd_imu, Cd_imu, R1_imu, R2_imu);
imu_kalman = ss( ...
    Phi_imu*(eye(2) - L_imu*Cd_imu), ...
    Phi_imu*L_imu, ...
    eye(2) - L_imu*Cd_imu, ...
    L_imu, ...
    Ts);

figure;
opt = bodeoptions;
opt.MagUnits = 'abs';
opt.MagScale = 'log';
bodemag(imu_kalman, opt); grid on;
bw_imu = bandwidth_mimo(imu_kalman)

%% Assignment 7-8

Phi_wheel = [1,0;h,1];    % Assignment 7
Gam_wheel = 0;    % Assignment 7
Cd_wheel  = [0,1];    % Assignment 7
Dd_wheel  = 0;    % Assignment 7
Gd_wheel  = [h;0];    % Assignment 7

R2_wheel  = 1;    % Assignment 8
R1_wheel  = 500000;    % Assignment 8

[L_wheel, ~, ~, Eig_wheel] = dlqe(Phi_wheel, Gd_wheel, Cd_wheel, R1_wheel, R2_wheel);
wheel_kalman = ss( ...
    Phi_wheel*(eye(2) - L_wheel*Cd_wheel), ...
    Phi_wheel*L_wheel, ...
    eye(2) - L_wheel*Cd_wheel, ...
    L_wheel, ...
    Ts);

figure;
opt = bodeoptions;
opt.MagUnits = 'abs';
opt.MagScale = 'log';
bodemag(wheel_kalman, opt); grid on;
bw_wheel = bandwidth_mimo(wheel_kalman)

%% Assignment 10

A = [-3.1,  58.4,   62.7,   0;
     1,     0,      0,      0;
     40.1,  -318,   -766,   0;
     0,     0,      1,      0];
B = [-148; 0; 1808; 0];
C = [1 0 0 0; 0 1 0 0; 0 0 0 1];
D = 0;

% Assignment 10: Discretize A and B into Phi and Gam
sysd = c2d(ss(A,B,C,D), Ts);
[Phi, Gam] = ssdata(sysd);

Eig_system = eig(Phi)   % Assignment 10: Analyze what the eigenvalues

%% Assignment 11
% m2 = 0.5*pi/180;   % Assignment 11
% m1 = m2/0.1;   % Assignment 11
% m4 = 10*pi/180;   % Assignment 11
% m3 = m4/10;   % Assignment 11
% 
% mu = 0.1;   % Assignment 11

m1 = 0.03;   % Assignment 11
m2 = 0.5*pi/180;   % Assignment 11
m3 = 10*pi;   % Assignment 11
m4 = 0.001;   % Assignment 11
mu = 0.00085;   % Assignment 11



Q1 = diag(1./[m1, m2, m3, m4].^2)
Q2 = 1/mu^2;

%% Assignment 12

feedback_gain = dlqr(Phi, Gam, Q1, Q2);
initial_sim_cl(Phi, Gam, feedback_gain, Ts);    % Assignment 12: Satisfy requirements
