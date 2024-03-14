% Example design of a controller for the linear servo.
% The user can specify the controller poles and zeros, and the
% desired crossover frequency, and controllers and closed-loop transfer
% functions are calculated.

s = tf('s');


%% Specify feedback controller here:
C_P= 0.12*(1+1/s);          % gain
omega_c = 6;
% -- Lag filter --
M= 10;            % M>1
a= 0.1*omega_c;
C_lag=(s+a)/(s+a/M);

% -- Notch filter --
gamma= 0.3;        % 0<gamma<1
b_n= 10;          % b_n>0
w_cr= 19;
C_notch=(s^2+gamma*b_n*s+w_cr^2)/(s^2+b_n*s+w_cr^2);

% -- Lead filter --
N= 2.7;           % N>1
bl= omega_c/sqrt(N);
C_lead=N*(s+bl)/(s+bl*N);


C= C_P*C_lead*C_lag*C_notch; %put together the controller e.g. C=C_P*C_lag;

%% Loop gain L(s), open-loop system
L=C*P;

T= L/(1+L);
Tr = 0.065;
d = 5;
%% Plot Bode diagram of L(s) with phase- and gain margins
figure(2)
hold off
margin(L);
grid

% Specify feedforward filter here:
F=(1/T)/((Tr*s + 1)^d); % Change to your desired feedforward filter!

% Define process model for simulation
P_sim=P;


% Uncomment these lines to run the Simulink model for the simulated plant
% and then evaluate the results
sim('servo_simulated')
specs;


% -----------------------------------------------------------------

% Calculate all relevant closed-loop transfer functions from the
% process and controller data.
S=1/(1+L); % Sensitivity function
G_yr=1-S; % r -> y
G_yd=feedback(P,C); % d -> y
G_yn=S; % n -> x
G_ur=feedback(C,P); % r -> u
G_un=-G_ur; % n -> u
G_ud=-G_yr; % d -> u