function initial_sim_cl(Phi, Gamma, K, Ts) 
%INITIAL_SIM_CL
% Simulates the closed loop system response given x0 = [0, 0.04, 0, 0]
% 
%   initial_sim_cl(Phi, Gamma, K, Ts);
%
% First plot:   x1
% Second plot:  x2
% Third plot:   x3
% Fourth plot:  x4
% Fifth plot:   u

x0 = [0, 0.04, 0, 0];
sys_cl = ss(Phi - Gamma*K, [], [eye(4); -K], 0, Ts);
initial(sys_cl, x0, 1)

end
