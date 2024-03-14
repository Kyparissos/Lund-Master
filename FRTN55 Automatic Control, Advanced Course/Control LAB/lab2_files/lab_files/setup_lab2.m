% Setup Arduino
setup_arduino

% Add paths to RASPlib
addpath([fileparts(mfilename('fullpath')) '/RASPlib'],...
		[fileparts(mfilename('fullpath')) '/RASPlib/src'],...
		[fileparts(mfilename('fullpath')) '/RASPlib/include'],...
		[fileparts(mfilename('fullpath')) '/RASPlib/examples'],...
		[fileparts(mfilename('fullpath')) '/RASPlib/blocks'],...
        [fileparts(mfilename('fullpath')) '/util'])

% Simulink model time-step
Ts = 0.01;

% Set default (zero) values for controller variables
imu_kalman          = ss(zeros(2), zeros(2), zeros(2), 0, Ts);    % Kalman filter 1
wheel_kalman        = ss(zeros(2), zeros(2,1), zeros(2), 0, Ts);  % Kalman filter 2
feedback_gain       = [0 0 0 0];                              % State feedback gain
integral_gain       = 0;                                      % Integral gain
