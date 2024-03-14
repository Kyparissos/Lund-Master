addpath('functions')
% Specify tank properties
lt = LinearTankModel;               % Create linearTankModel object
lt.Gammas = [0.7 0.7];              % Set gamma_1 and gamma_2
lt.LinPointLowerTanks = [10 10];    % Specify h_10 and h_20
% Set the sample time
lt.SampleTime = 0.1;  %Ts

% To be used after Task 4
% lt.estimData = load('estimData.mat').estimData;
% lt.useEstimData=...;

% Needed for the simulation
Ts = lt.SampleTime;

% Set MPC parameters for cost function
predictionHorizon   = 1;   % H_p
firstPenaltySample  = 1;

% Specify costs, Q1 = diag(q1), Q2 = diag(q2)
q1 = [10 10 0 0];
q2 = [1 1];
R1 = eye(4);
R2 = eye(2);

% Specify Constraints
pumpConstraints = struct('Max', {inf inf},...
    'Min', {-inf -inf}, 'RateMin', {-inf -inf}, 'RateMax', {inf inf});

tankConstraints = struct('Max', {inf inf inf inf}, ...
    'Min', {-inf -inf -inf -inf});

% Create an LGQ controller with white noise models for process and
% measurement nosie.
[CL_lqg_wn, L, Kff] = setup_lqg(lt, diag(q1), diag(q2), R1, R2);
% Create an MPC object
mpcobj = setup_mpc(lt, pumpConstraints, tankConstraints,q1, q2,...
    predictionHorizon, firstPenaltySample);

% Manually set the disturbance model and estimator to match that of the
% LQG controller.
% Comment out these two lines to run with the standard, integrated white
% noise as output disturbance model.
setoutdist(mpcobj, 'Model', ss(zeros(4,1)));
setEstimator(mpcobj, ssdata(lt.DtSys)*L, L);

