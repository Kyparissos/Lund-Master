clear;clc;
x=-3:0.01:3;
y=0*(x>2&x<-2)+(1-x.^2).*((-1<=x)&(x<=1))+(-2*abs(x).^3+10*x.^2-16*abs(x)+8).*((1<abs(x))&(abs(x)<=2));
dy=diff(y)/0.01;
diff(x)
plot(x(2:end),dy)
daspect([1 1 1])

% 
% daspect([1 1 1])
% f=[1 4 6 8 7 5 3];
% 
% for i=1:7
%     y(x==x-i)
% end
% plot(x,y)
% %F=conv(y,f,"same");
% %plot(x,F)
% 
% %%
% 
% x=0:0.01:1;
% mu1 = 0.4;
% mu2 = 0.32;
% mu3 = 0.55;
% sigma1 = 0.01;
% sigma2 = 0.05;
% sigma3 = 0.2;
% y1 = normpdf(x,mu1,sigma1);
% y2 = normpdf(x,mu2,sigma2);
% y3 = normpdf(x,mu3,sigma3);
% %j1 = integral(y1,0,1);
% %j2 = integral(y2,0,1);
% %j3 = integral(y3,0,1);
% y3(x==0.44)
% 
% plot(x,y1,x,y2,x,y3)
% %{
% f1 = cdf('Normal',x,mu1,sigma1);
% f2 = cdf('Normal',x,mu2,sigma2);
% f3 = cdf('Normal',x,mu3,sigma3);
% %plot(x,f1,x,f2,x,f3)
% 
% 0.3*0.8^10*0.2^6
% 0.2*0.8^12*0.2^4
% 0.2*0.8^10*0.2^6
% 0.3*0.8^8*0.2^8
% %}
% %% 
% 0.2^5*0.7^5*0.8^5;
% 0.2^5*0.7^5*0.8^3*0.3^2;
% 0.2^3*0.3*0.7^7*0.8^4;
% 1.7623e-05+2.4783e-06+8.0958e-05
% 1.7623e-05/1.0106e-04
% 2.4783e-06/1.0106e-04
% 8.0958e-05/1.0106e-04
% %% 
% 
% %%
% 0.2^5*0.7^5*0.8^5*0.35;
% 8.0958e-05*0.25;
% 2.4783e-06*0.4;
% 6.1682e-06+2.0240e-05+9.9132e-07
% 6.1682e-06/2.7400e-05
% %% 
