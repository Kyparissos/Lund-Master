function plotctrl(C)

subplot(1,2,1)
bodeplot(C(1,1),'b')
grid on
title('C_{11}')

subplot(1,2,2)
bodeplot(C(1,2),'b')
grid on
title('C_{12}')

end
