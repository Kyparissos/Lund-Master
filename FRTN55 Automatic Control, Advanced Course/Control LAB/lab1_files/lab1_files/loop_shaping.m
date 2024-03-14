% Example design of a controller for the linear servo.
% The user can specify the controller poles and zeros, and the
% desired crossover frequency, and controllers and closed-loop transfer
% functions are calculated.

s = tf('s');

%% Specify feedback controller here:
C_P= 0.045;          % gain
omega_c =3.4;
% -- Lag filter --
M= 10;            % M>1
a= 0.5*omega_c;
C_lag=(s+a)/(s+a/M);

% -- Notch filter --
gamma= 0.3;        % 0<gamma<1
b_n= 10;          % b_n>0
w_cr= 19;
C_notch=(s^2+gamma*b_n*s+w_cr^2)/(s^2+b_n*s+w_cr^2);

% -- Lead filter --
N= 3;           % N>1
bl= omega_c/sqrt(N);
C_lead=N*(s+bl)/(s+bl*N);


% C= C_P; %put together the controller e.g. C=C_P*C_lag;
C= C_P* C_lag* C_notch* C_lead;
%% Loop gain L(s), open-loop system
L=C*P;

T= L/(1+L);
Tr = 0.065;
d = 1;

%% Plot Bode diagram of L(s) with phase- and gain margins
figure(2)
hold off
margin(L);
grid

% Specify feedforward filter here:
F=(1/T)/((Tr*s + 1)^d); % Change to your desired feedforward filter!
% F= 1/((Tr*s + 1)^d);
% Define process model for simulation
P_sim=P;


% Uncomment these lines to run the Simulink model for the simulated plant
% and then evaluate the results
%sim('servo_simulated_matlab2012')
%specs;


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
%% 
sim('servo_simulated');
specs;
