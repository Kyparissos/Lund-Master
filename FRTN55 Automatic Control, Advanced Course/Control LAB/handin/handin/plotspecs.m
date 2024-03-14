function plotspecs(G)

h = 0.05;
t = 0:h:3;

if G.Ts ~= h
  error('The system does not have the correct sampling time')
end 

constr = get_constraints(t);

y = step(G,t);

subplot(2,3,1)
stairs(t,y(:,1,1),'b','Linew',1);
hold on
plot(t,constr.Gz0r_upper,'r--','Linew',1)
plot(t,constr.Gz0r_lower,'r--','Linew',1)
xlabel('Time (s)')
grid on
title('G_{11}')

subplot(2,3,2)
stairs(t,y(:,1,2),'b','Linew',1);
grid on
title('G_{12}')

subplot(2,3,3)
stairs(t,y(:,1,3),'b','Linew',1);
hold on
plot(t,constr.Gz0d_upper,'r--','Linew',1)
plot(t,constr.Gz0d_lower,'r--','Linew',1)
xlabel('Time (s)')
grid on
title('G_{13}')

subplot(2,3,4)
stairs(t,y(:,2,1),'b','Linew',1);
hold on
plot(t,constr.Gur_upper,'r--','Linew',1)
plot(t,constr.Gur_lower,'r--','Linew',1)
xlabel('Time (s)')
grid on
title('G_{21}')

subplot(2,3,5)
stairs(t,y(:,2,2),'b','Linew',1);
xlabel('Time (s)')
grid on
title('G_{22}')

subplot(2,3,6)
stairs(t,y(:,2,3),'b','Linew',1);
hold on
plot(t,constr.Gud_upper,'r--','Linew',1)
plot(t,constr.Gud_lower,'r--','Linew',1)
xlabel('Time (s)')
xlabel('Time (s)')
grid on
title('G_{23}')
