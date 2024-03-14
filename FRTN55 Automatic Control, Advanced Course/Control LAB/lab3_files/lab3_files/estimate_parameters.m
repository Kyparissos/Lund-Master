addpath('functions')
umax = 10;
ustat = 7;
stepstop = 10;

statstart1 = 30;
statstop1 = 100;
statstart2 = 100;
statstop2 = 170;

out = sim('estimate_parameters_sim');
%%
% figure(1)
% subplot(411)
% plot(out.y_out.time, out.y_out.data(:,1))
% title('y1')
% subplot(412)
% plot(out.y_out.time, out.y_out.data(:,2))
% title('y2')
% subplot(413)
% plot(out.y_out.time, out.y_out.data(:,3))
% title('y3')
% subplot(414)
% plot(out.y_out.time, out.y_out.data(:,4))
% title('y4')
% 
% figure(2)
% subplot(211)
% plot(out.u_out.time, out.u_out.data(:,1))
% title('u1')
% subplot(212)
% plot(out.u_out.time, out.u_out.data(:,2))
% title('u2')
%% Plot heights
t = out.y_out.time;
y = out.y_out.data;
lt = LinearTankModel;
kc = lt.MeasurementConstant;
h = y/kc;
% figure(3)
% plot(t,h)
% title('h')
% legend('h1','h2','h3','h4')
%% Denoise the heights
n = 50;
MACoeff = ones(1,n)/n;
MA = filter(MACoeff, 1, h);
fDelay = (length(MACoeff)-1)/2;
clear n MACoeff
% figure(4)
% plot(t-fDelay*.05,MA)
% legend('h1','h2','h3','h4')
% title('Denoised h')

%% Get derivatives
dh = (MA(2:end, :) - MA(1:end-1, :))./(t(2:end) - t(1:end-1));
% figure(5)
% plot(t(2:end)*[1 1], dh(:, 3:4));
% legend('dh3', 'dh4')
% title('dh')
%%
% First fill the upper tanks, and calculate a3/A3 and a4/A4
t1 = find(t>stepstop+5 & t<stepstop+9);
abyA3s = -dh(t1, 3)./sqrt(2*lt.Gravity*MA(t1,3));
abyA3 = mean(abyA3s);

abyA4s = -dh(t1, 4)./sqrt(2*lt.Gravity*MA(t1,4));
abyA4 = mean(abyA4s);

% using u1 = u1stat, wait for stationarity and calculate B41 and a2/A2
t2 = find(t > statstop1-5 & t< statstop1);
B41s =abyA4*sqrt(2*lt.Gravity*MA(t2,4))/ustat;
B41 = mean(B41s);
abyA2s = B41s*ustat./sqrt(2*lt.Gravity*MA(t2, 2));
abyA2 = mean(abyA2s);

t3 = find(t > statstop2-5 & t< statstop2);
B32s = abyA3*sqrt(2*lt.Gravity*MA(t3,3))/ustat;
B32 = mean(B32s);
abyA1s = B32s*ustat./sqrt(2*lt.Gravity*MA(t3, 1));
abyA1 = mean(abyA1s);

B11s = abyA1s.*sqrt(2*lt.Gravity*MA(t2, 1))/ustat;
B11 = mean(B11s);

B22s = abyA2s.*sqrt(2*lt.Gravity*MA(t3,2))/ustat;
B22 = mean(B22s);

clear abyA1s abyA2s abyA3s abyA4s...
    B11s B22s B32s B41s t1 t2 t3
clear umax ustat stepstop statstart1 statstop1...
    statstart2 statstop2
clear dh fDelay h kc lt MA out t y
%% Save data
estimData.abyA = [abyA1 abyA2 abyA3 abyA4];
estimData.B = [B11 0; 0 B22; 0 B32; B41 0];

save('estimData.mat','estimData');

clear abyA1 abyA2 abyA3 abyA4 B B11...
    B22 B32 B41 estimData
