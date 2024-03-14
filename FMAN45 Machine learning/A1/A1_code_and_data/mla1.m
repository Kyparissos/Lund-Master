%% Task 4
clear
clc
close all
load A1_data.mat

lambda1 = 0.1;
omega_hat1 = skeleton_lasso_ccd(t, X, lambda1);
y_hat1 = Xinterp * omega_hat1;
y_data1 = X * omega_hat1;

lambda2 = 10;
omega_hat2 = skeleton_lasso_ccd(t, X, lambda2);
y_hat2 = Xinterp * omega_hat2;
y_data2 = X * omega_hat2;

lambda3 = 2;
omega_hat3 = skeleton_lasso_ccd(t, X, lambda3);
y_hat3 = Xinterp * omega_hat3;
y_data3 = X * omega_hat3;

figure
hold on
scatter(n, t, 20, 'b')
scatter(n, y_data1, 20, 'filled')
plot(ninterp, y_hat1, 'r')
legend('Original Data Points', 'Reconstructed Data Points', 'Interpolatation')
xlabel('Time')
title('\lambda = 0.1')

figure
hold on
scatter(n, t, 20, 'b')
scatter(n, y_data2, 20, 'filled')
plot(ninterp, y_hat2, 'r')
legend('Original Data Points', 'Reconstructed Data Points', 'Interpolatation')
xlabel('Time')
title('\lambda = 10')

figure
hold on
scatter(n, t, 20, 'b')
scatter(n, y_data3, 20, 'filled')
plot(ninterp, y_hat3, 'r')
legend('Original Data Points', 'Reconstructed Data Points', 'Interpolatation')
xlabel('Time')
title('\lambda = 2')

%% 
non_zero_coordinates1 = sum(omega_hat1~=0);
non_zero_coordinates2 = sum(omega_hat2~=0);
non_zero_coordinates3 = sum(omega_hat3~=0);

%% Task 5
clear;clc
close all
load A1_data.mat

lambda_min = 0.01;
lambda_max = max(abs(X'*t));
N_lambda = 100;
N_folds = 10;

lambda_grid = exp( linspace( log(lambda_min), log(lambda_max), N_lambda));
[wopt, lambdaopt, RMSEval, RMSEest] = skeleton_lasso_cv(t, X, lambda_grid, N_folds);

t_optimal = Xinterp*wopt;
t_data = X*wopt;
%%
legends = [];
figure
hold on
legends = [legends scatter(log(lambda_grid), RMSEval, 'bx')];
plot(log(lambda_grid), RMSEval, 'b')
legends = [legends scatter(log(lambda_grid), RMSEest, 'gx')];
plot(log(lambda_grid), RMSEest, 'g')
legends = [legends xline(log(lambdaopt), '--r')];
legend(legends, 'RMSE Validation', 'RMSE Estimate', 'Optimal \lambda')
xlabel('log(\lambda)')

omega_hat_opt = skeleton_lasso_ccd(t, X, lambdaopt);
y_hat_opt = Xinterp * omega_hat_opt;
y_data_opt = X * omega_hat_opt;

figure
hold on
scatter(n, t, 20, 'b')
scatter(n, y_data_opt, 20, 'filled')
plot(ninterp, y_hat_opt, 'r')
legend('Original Data Points', 'Reconstructed Data Points', 'Interpolatation')
xlabel('Time')
title('\lambda = 1.1')

%% Task 6
clear;clc 
load('A1_data.mat')
%soundsc(Ttrain,fs)

lambda_min = 0.0001;
N_lambda = 100;
N_folds = 3;
frame_length = size(Xaudio, 1);
N_frames = floor(length(Ttrain)./frame_length);
lambda_maxes = zeros(N_frames, 1);

for frame_index=1:N_frames
    lambda_maxes(frame_index) = max(abs(Xaudio'*Ttrain(1 + frame_length*(frame_index-1) : frame_index*frame_length)));
end

lambda_max = max(lambda_maxes);
lambda_grid = exp( linspace( log(lambda_min), log(lambda_max), N_lambda));

[wopt, lambdaopt, RMSEval, RMSEest] = skeleton_multiframe_lasso_cv(Ttrain, Xaudio, lambda_grid, N_folds);
save("multiframe_results.mat", 'wopt', 'lambdaopt', "RMSEest", 'RMSEval');

%% PLOTS
load A1_data.mat
load multiframe_results.mat
legends = [];

figure
hold on
legends = [legends scatter(log(lambda_grid), RMSEval, 'bx')];
plot(log(lambda_grid), RMSEval, 'b')
legends = [legends scatter(log(lambda_grid), RMSEest, 'gx')];
plot(log(lambda_grid), RMSEest, 'g')
legends = [legends xline(log(lambdaopt), '--r')];
legend(legends, 'RMSE Validation', 'RMSE Estimate', 'Optimal \lambda')
xlabel('log(\lambda)')


%% TASK 7
clear;clc
load A1_data.mat
load multiframe_results.mat

Ytest = lasso_denoise(Ttest, Xaudio, lambdaopt);
soundsc(Ytest, fs);

save('denoised_audio','Ytest','fs');
%% 
Ytrain = lasso_denoise(Ttrain, Xaudio, 0.01);
%%
soundsc(Ttrain,fs)
%%
soundsc(Ytrain,fs)
%%
soundsc(Ttest,fs)
%%
soundsc(Ytest,fs)
