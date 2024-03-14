clear all
load coal_mine_disasters.mat
%% 1.b & 1.c
tStart = 1658;
tEnd = 1980;
tsIni = 1690;
teIni = 1980;


n = 751;
N = 10000;
d = 6;
tau = T;
burnIn = 1000;
psi = 20;
%%
for i = 2:d
    q = i;
    steps = (teIni - tsIni) / i;
    tMid = tsIni : steps : teIni;
    t = [tsIni, tMid(2:end - 1), teIni];
    tl = length(t);
    rhos = 0.01 * ones(i, 1);
    theta = gamrnd(2, 1/psi);
    lambda = gamrnd(2, 1/theta, 1, d);
    bp = zeros(N, length(t));
    
    for j = 1:burnIn
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        temp = zeros(1, tl - 1);
        for k = 1:tl-1
            temp(k) = sum((t(k) <= tau) & (tau < t(k + 1)));
        end
        lambda = gamrnd(temp' + 2, 1./(theta + (t(2:end) - t(1:end - 1))'));
        accept = zeros(1, tl-2);

        for a = 2:tl-1
            R = rhos(1) * (t(a + 1) - t(a - 1));
            Xstar =  t(a) - R + 2 * R * rand;

            while (Xstar < t(a - 1) || Xstar >  t(a + 1))
                Xstar =  t(a) - R + 2 * R * rand;
            end

            num = posterior(lambda, [t(1:a - 1), Xstar, t(a + 1:end)], tau);
            den = posterior(lambda, t, tau);
            A = min(1, num / den);
            U = rand(1);

            if (U <= A)
                t(a) = Xstar;
                accept(a - 1) = accept(a - 1) + 1;
            end    
        end
    end



    for j = 1:N
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        temp = zeros(1, tl - 1);
        for k = 1:tl-1
            temp(k) = sum((t(k) <= tau) & (tau < t(k + 1)));
        end
        lambda = gamrnd(temp' + 2, 1./(theta + (t(2:end) - t(1:end - 1))'));
        accept = zeros(1, tl-2);

        for a = 2:tl-1
            R = rhos(1) * (t(a + 1) - t(a - 1));
            Xstar =  t(a) - R + 2 * R * rand;

            while (Xstar < t(a - 1) || Xstar >  t(a + 1))
                Xstar =  t(a) - R + 2 * R * rand;
            end

            num = posterior(lambda, [t(1:a - 1), Xstar, t(a + 1:end)], tau);
            den = posterior(lambda, t, tau);
            A = min(1, num / den);
            U = rand(1);

            if (U <= A)
                t(a) = Xstar;
                accept(a - 1) = accept(a - 1) + 1;
            end    
        end

        bp(j, :) = t;
    end

%     figure
%     plot(bp)
%     title(['d = ' num2str(q-1) ' breakpoints'])
end

% figure
% plot(T, 1:length(T))
% hold on
% for j = 2:d
%     line([mean(bp(:, j)) mean(bp(:, j))], [0 length(T)], 'Color', [rand rand rand])
% end
% axis([tStart T(end) 0 length(T)])
%% 1.d
psi = 50;
rhos = 0.01 * ones(d, 1);

thetaMean = zeros(psi, 1);
thetaVar = zeros(psi, 1);
lambdaMean = zeros(psi, d);
lambdaVar = zeros(psi, d);
steps = (teIni - tsIni) / d;
tMid = tsIni : steps : teIni;
t = [tStart, tMid(2:end - 1), tEnd];
tl = length(t);
tArr = zeros(N, tl);
tMean = zeros(psi, tl);
tVar = zeros(psi, tl);
for i = 1:psi
    theta = gamrnd(2, 1 / i);
    lambda = gamrnd(2, 1 / theta, 1, d);
    
    thetaTemp = zeros(N, 1);
    lambdaTemp = zeros(N, d);


    for j = 1:burnIn
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (i + sum(lambda)));
        temp = zeros(1, tl - 1);
        for k = 1:tl-1
            temp(k) = sum((t(k) <= tau) & (tau < t(k + 1)));
        end
        lambda = gamrnd(temp' + 2, 1./(theta + (t(2:end) - t(1:end - 1))'));
        accept = zeros(1, tl-2);

        for a = 2:tl-1
            R = rhos(1) * (t(a + 1) - t(a - 1));
            Xstar =  t(a) - R + 2 * R * rand;

            while (Xstar < t(a - 1) || Xstar >  t(a + 1))
                Xstar =  t(a) - R + 2 * R * rand;
            end

            num = posterior(lambda, [t(1:a - 1), Xstar, t(a + 1:end)], tau);
            den = posterior(lambda, t, tau);
            A = min(1, num / den);
            U = rand(1);

            if (U <= A)
                t(a) = Xstar;
                accept(a - 1) = accept(a - 1) + 1;
            end    
        end
    end



    for j = 1:N
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (i + sum(lambda)));
        temp = zeros(1, tl - 1);
        for k = 1:tl-1
            temp(k) = sum((t(k) <= tau) & (tau < t(k + 1)));
        end
        lambda = gamrnd(temp' + 2, 1./(theta + (t(2:end) - t(1:end - 1))'));
        accept = zeros(1, tl-2);

        for a = 2:tl-1
            R = rhos(1) * (t(a + 1) - t(a - 1));
            Xstar =  t(a) - R + 2 * R * rand;

            while (Xstar < t(a - 1) || Xstar >  t(a + 1))
                Xstar =  t(a) - R + 2 * R * rand;
            end

            num = posterior(lambda, [t(1:a - 1), Xstar, t(a + 1:end)], tau);
            den = posterior(lambda, t, tau);
            A = min(1, num / den);
            U = rand(1);

            if (U <= A)
                t(a) = Xstar;
                accept(a - 1) = accept(a - 1) + 1;
            end    
        end

        thetaTemp(j) = theta;
        lambdaTemp(j, :) = lambda';
        tArr(j, :) = t;
    end
    
    thetaMean(i) = mean(thetaTemp);
    thetaVar(i) = var(thetaTemp);
    lambdaMean(i, :) = mean(lambdaTemp);
    lambdaVar(i, :) = var(lambdaTemp);
    tMean(i, :) = mean(tArr);
    tVar(i, :) = var(tArr);
end

figure()
subplot(1,2,1)
plot(thetaMean, "-x");
title("The Mean of \theta with different \psi");
grid on;
hold on

subplot(1,2,2)
plot(thetaVar, "-x");
title("The Variance of \theta with different \psi");
grid on;

figure()
subplot(1,2,1)
plot(lambdaMean, "-x");
title("The Mean of \lambda with different \psi");
legend("\lambda_1", "\lambda_2", "\lambda_3", "\lambda_4", "\lambda_5", "\lambda_6");
hold on

subplot(1,2,2)
plot(lambdaVar, "-x");
title("The Variance of \lambda with different \psi");
legend("\lambda_1", "\lambda_2", "\lambda_3", "\lambda_4", "\lambda_5", "\lambda_6");

figure()
subplot(1,2,1)
plot(tMean(:, 2:end - 1), "-x");
title("The Mean of the breakpoints with different \psi");
legend("t_1", "t_2", "t_3", "t_4", "t_5");
hold on

subplot(1,2,2)
plot(tVar(:, 2:end - 1), "-x");
title("The Variance of the breakpoints with different \psi");
legend("t_1", "t_2", "t_3", "t_4", "t_5");
%% Acceptance Ratio
rhoPlot = linspace(.001, .1);
acceptedPlot = zeros(length(rhoPlot), d - 1);
accSamp = 100;

for r = 1:length(rhoPlot)
   for i = 1:accSamp
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end-1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rhoPlot(r) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = posterior(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = posterior(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
       
       acceptedPlot(r, :) = acceptedPlot(r, :) + accepted;
   end
end

ratio = sum(acceptedPlot, 2) / (accSamp * (d - 1));

figure();
plot(rhoPlot, ratio);
hold on;
plot([.001, .1],[.3, .3]);
title("Acceptance Rate with respect to \rho");
%% 2.e - theta & lambda
nor = 30;
psi = 20;

rhoTemp = (1:nor) * 0.01;
rho = zeros(d, nor);

for i = 1:d
   rho(i, :) = rhoTemp; 
end

thetaMean = zeros(nor, 1);
lambdaMean = zeros(nor, d);

for n = 1:nor
    theta = gamrnd(2, 1 / psi);
    lambda = gamrnd(2, 1 / theta, 1, d);
    thetaTemp = zeros(N, 1);
    lambdaTemp = zeros(N, d);
    
    for i = 1:burnIn
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end-1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rho(1, n) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = posterior(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = posterior(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
    end
    
    for i = 1:N
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end-1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rho(1, n) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = posterior(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = posterior(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
        
        thetaTemp(i) = theta;
        lambdaTemp(i, :) = lambda';
    end
    
    thetaMean(n) = mean(thetaTemp);
    lambdaMean(n, :) = mean(lambdaTemp);
end

figure()
subplot(1,2,1)
plot(rho(1, :), thetaMean, "-x");
title("The Mean of \theta with Different \rho");

subplot(1,2,2)
plot(rho(1, :), lambdaMean, "-x");
title("The Mean of \lambda with Different \rho");
legend("\lambda_1", "\lambda_2", "\lambda_3", "\lambda_4", "\lambda_5", "\lambda_6");

%% 2.e - t
rArr = [.02, .03, .04, .05];
nor = length(rArr);
rho = zeros(d, nor);

for i = 1:d
   rho(i, :) = rArr; 
end

tArr = zeros(N, length(t));

for n = 1:nor
    tArr = zeros(N, length(t));
    theta = gamrnd(2, 1 / psi);
    lambda = gamrnd(2, 1 / theta, 1, d);
    
    for i = 1:burnIn
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end-1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rho(1, n) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = posterior(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = posterior(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
    end
    
    for i = 1:N
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end-1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rho(1, n) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = posterior(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = posterior(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
        
        tArr(i, :) = t;
    end
    
    figure();
    subplot(2, 2, 1);
    autocorr(tArr(:, 2), 5e2)
    title("The Correlation Function for t_1");
    xlabel("Time lag");
    ylabel("Dependency on \rho = " + rho(1, n));
    
    subplot(2, 2, 2)
    autocorr(tArr(:, 3), 5e2)
    title("The Correlation Function for t_2");
    xlabel("Time lag");
    ylabel("Dependency on \rho = " + rho(1, n));
    
    subplot(2, 2, 3)
    autocorr(tArr(:, 4), 5e2)
    title("The Correlation Function for t_3");
    xlabel("Time lag");
    ylabel("Dependency on \rho = " + rho(1, n));
    
    subplot(2, 2, 4)
    autocorr(tArr(:, 5), 5e2)
    title("The Correlation Function for t_4");
    xlabel("Time lag");
    ylabel("Dependency on \rho = " + rho(1, n));
    
    if (n == 3)
        figure();
        plot(tArr);
        title("The Behavior of the Chain For " + (d - 1) + " Breakpoints with \rho = " + rho(1, n));
        xlabel("N");
        ylabel("Year");
    end
end
%%
clear
clc
load atlantic.txt
% F inverse
f_inv = @(u, mu, beta) mu - beta * log(-log(u));
%Parameter estimates
[betaHat,muHat] = est_gumbel(atlantic);
n = length(atlantic);
B = 1000;
boot_beta = zeros(B,1);
boot_mu = zeros(B,1);
% bootstrap
for b = 1:B
    U = rand(n,1);
    boot = f_inv(U, muHat, betaHat);
    [beta_boot, mu_boot]  = est_gumbel(boot);
    boot_beta(b) = beta_boot;
    boot_mu(b) = mu_boot;
end
delta_beta = sort(betaHat - boot_beta);
delta_mu = sort(muHat - boot_mu);
% compute 95% confidence interval
alpha = 0.05;
Beta_lower = betaHat - delta_beta(ceil((1 - alpha / 2) * B));
Beta_upper = betaHat - delta_beta(ceil(alpha * B / 2));
Mu_lower = muHat - delta_mu(ceil((1 - alpha / 2) * B));
Mu_upper = muHat - delta_mu(ceil(alpha * B / 2));
%% 
T = 3 * 14 * 100;
boot_wave = zeros(B,1);
wave_hat = f_inv(1 - 1 / T, muHat, betaHat);
for b = 1:B
      boot_wave(b,1) = f_inv(1 - 1 / T, boot_mu(b), boot_beta(b));
end

delta_wave = sort(boot_wave - wave_hat);
disp(wave_hat - delta_wave(ceil(alpha * B)))

