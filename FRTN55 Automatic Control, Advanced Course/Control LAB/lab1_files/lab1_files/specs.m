%
% Plots results with specifications.
%

figure(21)
plot(y_out.time,y_out.signals.values);

xx =[ 1.0
    1.2
    1.6
    2.1
    3.0
    5.0
    10.0];

yy =[ 0.0
    5.3
    5.5
    5.5
    5.1
    5.1
    5.1];

x =[ 1.25
    1.6
    1.9
    5.0
    5.3
    6.2
    6.9
    10.0];


y = [0.0
    4.4
    4.9
    4.9
    4.4
    4.4
    4.9
    4.9];

%xspec1 = [0.1500  0.5500    0.8000    3.0000];
%yspec1 = [0       0.8500    0.9500    0.9500];
%xspec2 = [0       0.3500    0.5500    0.8000    3.0000];
%yspec2 = [0       1.2000    1.2000    1.0500    1.0500];

hold on
plot(xx,yy,'b--',x,y,'b--')
hold off
title('Specifications for y(t).')
xlabel('t [s]')
ylabel('y [V]')

figure(22)
plot(u_out.time,u_out.signals.values);

x =[0.0
    0.85
    0.90
    1.1
    1.8
    5.0
    5.2
    6.0
    6.6
    10.0];

y =[0.15
    0.15
    0.8000
    0.8000
    0.15
    0.15
    0.30
    0.30
    0.25
    0.25];

xx =[0.0
    1.00
    1.1
    2.0
    2.5
   10.0];

yy =[-0.15
   -0.15
   -0.25
   -0.25
   -0.15
   -0.15];



hold on
plot(x,y,'b--',xx,yy,'b--')
hold off
title('Specifications for u(t).')
xlabel('t [s]')
ylabel('u [V]')
