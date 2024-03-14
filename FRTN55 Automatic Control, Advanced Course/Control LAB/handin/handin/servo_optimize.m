clear all

%%% Script for hand-in assignment on control design using convex optimization

%%% Part 1 %%%

%%% Task 0: Fill in your name(s) here: 

% Define the process model and sample it.
servo_model           % Defines the plant state-space matrices A,B,C,D

%%% Task 1: Sample the process model with the interval h = 0.05. The resulting
%%% sampled state and input matrices should be called Phi and Gam

h = 0.05;
sysd = c2d(ss(A,B,C,D),h);
[Phi,Gam] = ssdata(sysd);
%% 

%%% Task 2: Enter the generalized plant model below

% Formulate the generalized plant state-space model
Gamw = [zeros(4,2),Gam];
Gamu = Gam;
Cz = [C;zeros(1,4)];
Dzw = [0 0 0; 0 0 0];
Dzu = [0;1];
Cy = [zeros(1,4);C];
Dyw = [1 0 0;0 -1 0];
Dyu = [0;0];
Pgen = ss(Phi,[Gamw,Gamu],[Cz;Cy],[Dzw,Dzu;Dyw,Dyu],h);
% [Phip,Gamp] = ssdata(Pgen);
%% 

%%% Task 3: Design some stabilizing LQG controller for Pgen
G = eye(4);
R1 = eye(4);
R2 = eye(2);
L = dlqe(Phi,G,Cy,R1,R2);       % G matches x, R1 matches x, R2 matches y
%% 
Q1 = eye(4);
Q2 = 10e-6;
K = dlqr(Phi,Gam,Q1,Q2);       % Q1 matches x, Q2 matches u
%% 

% Formulate the stabilized system T (see Lecture 13)
Pyu = Pgen(3:4,4);      % Extract subsystem Pyu of generalized plant
ctrl_lqg = ss(Phi-Gam*K-Phi*L*Pyu.C,[Phi*L Gam],[-K; -Pyu.C], ... 
  [0 0 1; eye(2) zeros(2,1)],h);  % Formulate LQG controller
Tgen = lft(Pgen,ctrl_lqg,1,2);    % Connect generalized plant and LQG controller

% Extract subsystems of stabilized generalized plant
Tzw = Tgen(1:2,1:3);    % from inputs 1:3 (w) to outputs 1:2 (z)
Tzu = Tgen(1:2,4);      % from input 4 (u) to outputs 1:2 (z)
Tyw = Tgen(3:4,1:3);    % from inputs 1:3 (w) to outputs 3:4 (y)
Tyu = Tgen(3:4,4);      % from input 4 (u) to outputs 3:4 (y)
figure(1);clf
plotspecs(Tzw)          % Check result against time-domain specs

%%% Task 4: Try some different static values of Q. Comment on how it worked
%%% here:

Q = [1,1];        % Define static Q parameter (1x2 vector)
Gzw = Tzw + Tzu*Q*Tyw;  % Calculate closed loop
figure(1);clf
plotspecs(Gzw)          % Check result against Lab 1 time-domain specs

return % Remove this to activate Part 2 below

%% 

%%% Part 2 %%%

Nq = 30;                % Number of basis functions in the Q filter
Nt = 120;               % Number of time points in the step response
Nw = 120;               % Number of frequency points in the freq response

% Create Q filter basis functions (two parallel FIR filters)
z = tf('z',h);
Qb = {};
for i = 0:Nq-1
  Qb{i+1} = [1/z^i 0];
  Qb{Nq+i+1} = [0 1/z^i];
end

% Setup and solve the convex optimization problem using CVX
cvx_solver sdpt3        % Alternatively, use sedumi
cvx_begin

  variable q(2*Nq)      % Create vector of optimization variables

  t = 0:h:h*Nt;         % Vector of time points for evaluation
  w = 0:pi/h/Nw:pi/h;   % Vector of frequency points for evaluation

  % Compute closed-loop step and frequency responses as functions of q.
  % The results end up in the 3D matrices Gzw_stepresp and Gzw_freqresp
  % (These lines of code are quite difficult to dechiper, we know...)
  Tzw_step = shiftdim(step(Tzw,t),1);
  Tzw_freq = freqresp(Tzw,w);
  for i = 1:2*Nq
    TzuQTyw_step{i} = shiftdim(step(Tzu*Qb{i}*Tyw,t),1);
    TzuQTyw_freq{i} = freqresp(Tzu*Qb{i}*Tyw,w);
  end
  Gzw_stepresp = Tzw_step;
  Gzw_freqresp = Tzw_freq;
  for i = 1:2*Nq
    Gzw_stepresp = Gzw_stepresp + TzuQTyw_step{i}*q(i);
    Gzw_freqresp = Gzw_freqresp + TzuQTyw_freq{i}*q(i);
  end
  
  constr = get_constraints(t); % Get a struct with time-domain constraints

  %%%% Task 5: State the correct optimization objective below
  
  minimize sum_square(Gzw_stepresp(?,?,:))

  subject to:

  %%%% Task 6: Add the correct constraint on the response from n to u below
  
  % sum_square(Gzw_stepresp(?,?,:))*h <= LIMIT
  
  %%% Task 7: Add the correct time-domain constraints below
  
  % Upper and lower constraints on step response from r to z0:
  % constr.Gz0r_lower <= squeeze(Gzw_stepresp(?,?,:)) <= constr.Gz0r_upper;
  
  % Upper and lower constraints on step response from r to u:
  % constr.Gur_lower <= squeeze(Gzw_stepresp(?,?,:)) <= constr.Gur_upper;

  % Upper and lower constraints step response from d to z0:
  % constr.Gz0d_lower <= squeeze(Gzw_stepresp(?,?,:)) <= constr.Gz0d_upper;

  % Upper and lower constraints on step response from d to u:
  % constr.Gud_lower <= squeeze(Gzw_stepresp(?,?,:)) <= constr.Gud_upper;

  %%% Task 8: Experiment with a constraint on Ms and comment below
 
  % abs(1 - Gzw_freqresp(?,?,:)) <= MS;
  
  %%% Task 9: Experiment with Nq and comment on the results here:
 
  cvx_end                 % Run the optimization

% Check result of optimization
if isnan(cvx_optval) || isinf(cvx_optval)
  error('Optimization failed - problem infeasible or numerically unstable')
end

% Make a state-space realization of Q and compute the closed-loop system Gzw
if Nq == 1
  Q1 = ss(q(1)); Q2 = ss(q(2));
else
  Q1 = ss(diag(ones(Nq-2,1),-1),[1;zeros(Nq-2,1)],q(2:Nq)',q(1),h);
  Q2 = ss(diag(ones(Nq-2,1),-1),[1;zeros(Nq-2,1)],q(Nq+2:2*Nq)',q(Nq+1),h);
end
Q = [Q1 Q2];

Gzw = Tzw + Tzu*Q*Tyw;  % The optimized closed-loop system

% Compute the optimized controller (see Lecture 12) and plot its Bode diagram
[Phiq,Gamq,Cq,Dq] = ssdata(ss(Q));
Phic = [Phi-Gamu*K-Phi*L*Cy-Gamu*Dq*Cy Gamu*Cq; -Gamq*Cy Phiq];
Gamc = [Phi*L+Gamu*Dq; Gamq];
Cc = [-(K+Dq*Cy) Cq];
Dc = Dq;
ctrl_opt = ss(Phic,Gamc,Cc,Dc,h);

% Compute the maximum sensitivity and maximum complementary sensitivity
Ms = norm(feedback(1,-ctrl_opt*Pyu),inf)
Mt = norm(feedback(-ctrl_opt*Pyu,1),inf)

% Plot the controller Bode diagrams
figure(3);clf
plotctrl(ctrl_opt)

% Plot the closed-loop Bode magnitude diagrams
figure(2);clf
plotbodemag(Gzw)

% Plot the closed-loop step responses, with specifications
figure(1);clf
plotspecs(Gzw)
