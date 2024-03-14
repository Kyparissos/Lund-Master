%% discritizing

syms s h k;
% h = 1; k = 1;
A = [2];
B = [1];
C = [1];
D = [0];
phi = expm(A*h)
gam = int(expm(A*s)*B,s,0,h)
% sys = ss(A,B,C,D,h);
% bode(sys)
g = C*inv(s*eye(2)-A)*B + D;
g = C*inv(s-phi)*gam
%% kalman riccati
syms a p1 p2 p3 p4
% phi = [0.5 0;-1 0.5];
% G = [3;0];
% p = dlyap(phi,G*G')
% d = (phi*phi'-eye(2))\(-G*G') % wrong!

A = [0.5 a; 0 0];
G = [1 0;1 1];
C = [1 0];
R1 = eye(2);
R2 = 1;
P = [p1 p2;p2 p3];
lyap = (A*P*A.')+(G*R1*G');
ricca = lyap - (A*P*C.')*inv(C*P*C.'+R2)*(A*P*C.').'
L = (P*C.')/(C*P*C.'+R2)


%% solve falman filter
A = [2];
G = [1];
C = [1];
Q = 1;
R = 1;
L = dlqe(A,G,C,Q,R)
[M,P,Z,E]= dlqe(A,G,C,Q,R)
%% solve LQR
A = 2;
B = 1;
Q1 = 1;
Q2 = 2;
[K,S,e] = dlqr(A,B,Q1,Q2)
A-K*B
%% 
g = (s-50)/((s-3)*(s+3));
b = 2/(s+2);
t = feedback(0.1*g,1);
bode(t,b)
legend('t','b')