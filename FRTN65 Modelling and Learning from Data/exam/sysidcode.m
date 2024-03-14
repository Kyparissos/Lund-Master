% Exam 2023 sysid-problem
close all
clear

load sysiddata.mat

z = iddata(y,u,h);

figure(1)
subplot(211)
plot(t,y)
ylabel('y')
subplot(212)
plot(t,u)
ylabel('u')
xlabel('t')

%%
% This step response was taken on the real system
% Can be used as test data for model validation
% 

figure(2)
plot(t,ystep)
title('Measured unit step response (to be used as test data)')

%%
% This evaluates different model orders

%NN = [(1:15)'  (1:15)' ones(15,1)];
%V = arxstruc(z,z,NN);
%Nbest = selstruc(V);

%% Study some ARX models of different orders
% You will use the step response as test data
% So you do not need to split y and u into train and test here


model_arx3 = arx(z,[3,3,1]);
model_arx5 = arx(z,[5 5 1]);
model_arx10 = arx(z,[10 10 1]);
oemodel = oe(z,[4,4,2]);
figure(3)
pred_horizon = inf;
compare(z,oemodel, model_arx3,model_arx5,model_arx10,pred_horizon)

%% The ARX models
%present(model_arx3)
%present(model_arx5)
%present(model_arx10)
present(oemodel)
%%
figure(4)
bode(oemodel,model_arx3,model_arx5,model_arx10,{0.1,pi/h})
legend('oe','arx3','arx5','arx10')
grid on
%fixfig

%% Spectrum of u and y

figure(5)
pwelch(u)
hold on
pwelch(y)

%%  Checking the residuals

figure(6)
%resid(z,model_arx3)
%resid(z,model_arx5)
resid(z,oemodel)

%% Testing the estimated models on the step response data
figure(7)
y3 = step(model_arx3,t);
y5 = step(model_arx5,t);
y10 = step(model_arx10,t);
yoe = step(oemodel,t);
plot(t,ystep,t,y3,t,y5,t,yoe)
[~, hobj] = legend('ystep','arx3','arx5','oe','FontSize',16,'Location','SouthEast')
set(findobj(hobj,'type','line'),'LineWidth',2.5);
title('Step Response Evaluation')
axis([0 1 -1 3])

% Hmm, this does not work well !




