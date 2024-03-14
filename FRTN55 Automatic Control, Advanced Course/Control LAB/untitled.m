% clc;clear;
s = tf('s');
% c = 1/(s^2-1*s+1);
p = -(s-1)/(s^2+1*s+1); % RHP zero
d = s/(s+2); % high pass filter
f = 1/(s*(s+50)); % integrator
% p = 1/(s*(s+6)*(s+3));
% c = 42*(4*s+6)/(s+6);
p = (s-1)/(s-2);
c = 30/(s-20);
s = (s-2)/(s+2);
% p = p*c;
% p = feedback(p*c,1);
p = feedback(1,c*s);
zpk(p)
% pole(p)
figure(1)
step(p)
figure(2)
margin(p)
% figure(3)
% margin(p*c)
% nyquist(p)
% hold all
% margin(d)
hold off
%% controllability & observability
A = [0];
B = [1];
C = [1];
co = ctrb(A,B);
ob = obsv(A,C);
r = rank(co);
unco = length(A) - rank(co)
unob = length(A)-rank(ob)
%% 

fplot(@(s) 2/(s+3));
hold off
%% 
p = (s-2)/(s+2);
bode(p)
%%
syms x
eqn = 1*x^2 + -4*x -1 == 0
s = vpasolve(eqn)