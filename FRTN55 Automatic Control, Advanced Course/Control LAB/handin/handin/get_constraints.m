function constr = get_constraints(times)
% Return time-domain constraints for the flexible servo for the time points 
% specified in the vector times

% limits on step response from r to z0
t = [0 0.25 0.6 1.1 2 10];
z = [0 1.06 1.1 1.1 1.02 1.02];
constr.Gz0r_upper = interp1(t,z,times)';
t = [0 0.25 0.6 0.9 10];
z = [0 0 0.88 0.98 0.98];
constr.Gz0r_lower = interp1(t,z,times)';

% limits on step response from r to u
t = [0 0.2 0.5 1 10];
u = [0.2 1.0 1.0 0.2 0.2];
constr.Gur_upper = interp1(t,u,times)';
t = [0 0.2 1 1.5 10];
u = [-0.2 -0.4 -0.4 -0.2 -0.2];
constr.Gur_lower = interp1(t,u,times)';

% limits on step response from d to z0
t = [0 0.3 1.2 1.9 10];
z = [0.22 1.9 1.9 0.2 0.2];
constr.Gz0d_upper = interp1(t,z,times)';
t = [0 10];
z = [-0.2 -0.2];
constr.Gz0d_lower = interp1(t,z,times)';

% limits on step response from d to u
t = [0 10];
u = [0.1 0.1];
constr.Gud_upper = interp1(t,u,times)';
t = [0 0.2 1 1.6 10];
u = [0 -1.5 -1.5 -1.2 -1.2];
constr.Gud_lower = interp1(t,u,times)';
