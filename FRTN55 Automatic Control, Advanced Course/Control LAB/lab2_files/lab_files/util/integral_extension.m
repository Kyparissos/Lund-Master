function [feedback_gain, integral_gain] = integral_extension(Phi, Gamma, Q1, Q2, Qi, Ts) 
%INTEGRAL_EXTENSION
% extends the statespace with an extra integral state.
% Calculates the feedback gain and the integral states feedback gain and return
% them as feedback_gain and integral_gain

Q1e = blkdiag(Q1, Qi);                   % Extended Q1 matrix
sys = ss(Phi, Gamma, [0 0 0 1], 0, Ts);  % Define system with theta as the only output
Ke = lqi(sys, Q1e, Q2);                  % Calculate extended feedback gain vector
feedback_gain = Ke(1:4);                 % Extract K
integral_gain = Ke(5);                   % Extract ki

end
