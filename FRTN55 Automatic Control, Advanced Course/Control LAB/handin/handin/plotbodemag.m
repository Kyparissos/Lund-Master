function plotbodemag(G)

w = logspace(-2,log10(pi/G.Ts));
mag = bode(G,w);

subplot(2,3,1)
loglog(w,abs(squeeze(mag(1,1,:))),'b','Linew',1)
grid on
axis([0.01 100 0.001 10])
xlabel('Frequency (rad/s)')
title('G_{11}')

subplot(2,3,2)
loglog(w,abs(squeeze(mag(1,2,:))),'b','Linew',1)
grid on
axis([0.01 100 0.001 10])
xlabel('Frequency (rad/s)')
title('G_{12}')

subplot(2,3,3)
loglog(w,abs(squeeze(mag(1,3,:))),'b','Linew',1)
grid on
axis([0.01 100 0.001 10])
xlabel('Frequency (rad/s)')
title('G_{13}')

subplot(2,3,4)
loglog(w,abs(squeeze(mag(2,1,:))),'b','Linew',1)
grid on
axis([0.01 100 0.001 10])
xlabel('Frequency (rad/s)')
title('G_{21}')

subplot(2,3,5)
loglog(w,abs(squeeze(mag(2,2,:))),'b','Linew',1)
grid on
axis([0.01 100 0.001 10])
xlabel('Frequency (rad/s)')
title('G_{22}')

subplot(2,3,6)
loglog(w,abs(squeeze(mag(2,3,:))),'b','Linew',1)
grid on
axis([0.01 100 0.001 10])
xlabel('Frequency (rad/s)')
title('G_{23}')

end
