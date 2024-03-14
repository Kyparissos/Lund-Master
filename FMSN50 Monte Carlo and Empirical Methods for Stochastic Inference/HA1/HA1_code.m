clear all;
load powercurve_V164.mat
month = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"];
lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6];
k = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];
N = 10000;
steps = 100;
sqN = sqrt(N);
se = 1.96;
%% 2.a -- Standard Monte Carlo
disp("------------Standard Monte Carlo------------")
tau = zeros(12, N/steps);
sd = zeros(12, N/steps);
LB = zeros(12, N/steps);
UB = zeros(12, N/steps);
width = zeros(1,12);
for i = 1:12
    q = 1; 
    for j = 100:steps:N
        number = wblrnd(lambda(i), k(i), j, 1);
        power = P(number);
        tau(i,q) = mean(power);
        sd(i,q) = std(power);
        LB(i,q) = tau(i,q) - (se * sd(i,q) / sqrt(j));
        UB(i,q) = tau(i,q) + (se * sd(i,q) / sqrt(j));
        width(1,i) = UB(i,end) - LB(i,end);
        q = q + 1;
    end

    disp(month(i) + ": LB = " + LB(i,end) + ", UB = " + UB(i,end) + "; Width = " + width(1,i));
end
%% 2.a -- Truncated version
disp("------------Truncated version------------");
a = 3.5;
b = 25;
tau_T = zeros(12, N/steps);
sd_T = zeros(12, N/steps);
LB_T = zeros(12, N/steps);
UB_T = zeros(12, N/steps);
width_T = zeros(1,12);
for i = 1:12
    q = 1;
    for j = 100:steps:N
        u = rand(j,1);
        Fxa = wblcdf(a, lambda(i), k(i));
        Fxb = wblcdf(b, lambda(i), k(i));
        temp = Fxa + u * (Fxb - Fxa);
        speed = wblinv(temp, lambda(i), k(i));
        power = P(speed) * (Fxb - Fxa);
        tau_T(i,q) = mean(power);
        sd_T(i,q) = std(power);
        LB_T(i,q) = tau_T(i,q) - (se * sd_T(i,q) / sqrt(j));
        UB_T(i,q) = tau_T(i,q) + (se * sd_T(i,q) / sqrt(j));
        width_T(1,i) = UB_T(i,end) - LB_T(i,end);
        q = q + 1;
    end
    disp(month(i) + ": LB = " + LB_T(i,end) + ", UB = " + UB_T(i,end) + "; Width = " + width_T(1,i));
end
%% 2.a -- comparison
figure(1)
x = 100:steps:N;
y = mean(tau(1,:));
const = @(x)(y).*x.^(0);
plot(x,const(x),'green');
hold on
plot(x, LB(1,:),'--','color','blue');
plot(x, UB(1,:),'--','color','blue');
plot(x, LB_T(1,:),'--','color','red');
plot(x, UB_T(1,:),'--','color','red');
%% 2.b
disp("------------Control Variate------------");
tau_CV = zeros(12, N/steps);
sd_CV = zeros(12, N/steps);
LB_CV = zeros(12, N/steps);
UB_CV = zeros(12, N/steps);
width_CV = zeros(1,12);
for i = 1:12
    q = 1;
    for j = 100:steps:N
        X = wblrnd(lambda(i), k(i), j, 1);
        Y = X;
        c = cov(P(X), Y);
        beta = -c(1,2)./var(Y);
        m = gamma(1+1/k(i))*lambda(i);
        Z = P(X) + beta * (Y - m);
        tau_CV(i,q) = mean(Z);
        sd_CV(i,q) = std(Z);
        LB_CV(i,q) = tau_CV(i,q) - (se * sd_CV(i,q) / sqrt(j));
        UB_CV(i,q) = tau_CV(i,q) + (se * sd_CV(i,q) / sqrt(j));
        width_CV(1,i) = UB_CV(i,end) - LB_CV(i,end);
        q = q + 1;
    end
        disp(month(i) + ": LB = " + LB_CV(i,end) + ", UB = " + UB_CV(i,end) + "; Width = " + width_CV(1,i));
end
%% 2.c
disp("------------Important Sampling------------");
tau_IS = zeros(12, N/steps);
sd_IS = zeros(12, N/steps);
LB_IS = zeros(12, N/steps);
UB_IS = zeros(12, N/steps);
width_IS = zeros(1,12);

for i = 1:12
    q = 1;
    for j = 100:steps:N
        x = 0:0.5:30;
        fx = @(x) wblpdf(x, lambda(i), k(i));
        phi = @(x) P(x);
        temp = phi(x) .* fx(x)';
        [M,I] = max(temp);
        mu = x(I);
        sigma = M/1e5;
        gx = @(x) normpdf(x,mu,sigma);
        IS = @(x) (fx(x) .* P(x)') ./ gx(x);
        X = sigma.*randn(1,j) + mu;
        ZZ = IS(X);
        tau_IS(i,q) = mean(ZZ);
        sd_IS(i,q) = std(ZZ);
        LB_IS(i,q) = tau_IS(i,q) - (se * sd_IS(i,q) / sqrt(j));
        UB_IS(i,q) = tau_IS(i,q) + (se * sd_IS(i,q) / sqrt(j));
        width_IS(1,i) = UB_IS(i,end) - LB_IS(i,end);
        q = q + 1;
    end
    disp(month(i) + ": LB = " + LB_IS(i,end) + ", UB = " + UB_IS(i,end) + "; Width = " + width_IS(1,i));
end
%% 2.d
disp("------------Antithetic Sampling------------");
tau_AS = zeros(12, N/steps);
sd_AS = zeros(12, N/steps);
LB_AS = zeros(12, N/steps);
UB_AS = zeros(12, N/steps);
width_AS = zeros(1,12);

for i = 1:12
    q = 1;
    for j = 100:steps:N
        u = rand(j/2, 1); % Half the sample size
        V = P(wblinv(u, lambda(i), k(i)));
        VTilde = P(wblinv(1 - u, lambda(i), k(i)));
        W = (V + VTilde) ./ 2;
        tau_AS(i,q) = mean(W);
        sd_AS(i,q) = std(W);
        LB_AS(i,q) = tau_AS(i,q) - (se * sd_AS(i,q) / sqrt(j));
        UB_AS(i,q) = tau_AS(i,q) + (se * sd_AS(i,q) / sqrt(j));
        width_AS(1,i) = UB_AS(i,end) - LB_AS(i,end);
        q = q + 1;
    end
        disp(month(i) + ": LB = " + LB_AS(i,end) + ", UB = " + UB_AS(i,end) + "; Width = " + width_AS(1,i));
end
%% 2.e
disp("------------Probability the Turbine Delivers Power------------");
for i = 1:12
    v = wblrnd(lambda(i), k(i), N, 1);
    power = P(v);
    probP = nnz(power) / numel(power); 
    disp(month(i) + ": Probability of power = " + probP)
end
%% 2.f
disp("------------Average Ration------------");
d = 164;
rho = 1.225;

for i = 1:12
    v = wblrnd(lambda(i), k(i), N, 1);
    power = P(v);
    powerAct = mean(power);
    powerTot = 0.5 * rho * pi * 0.25 * d^2 * gamma(1 + 3 / k(i))*lambda(i)^3;
    powerRat = powerAct / powerTot;
    sd = std(power / powerTot);
    LB = powerRat - (se * sd / sqN);
    UB = powerRat + (se * sd / sqN);
    width = UB - LB;
    disp(month(i) + ": LB = " + LB + ", UB = " + UB + "; Width = " + width)
end
%% 2.e
disp("------------Capacity Factor & Availability Factor------------");
capacityF = zeros(1,12);
avaliaF = zeros(1,12);
for i = 1:12
    v = wblrnd(lambda(i), k(i), N, 1);
    power = P(v);
    capacityF(1,i) = mean(power) / 9.5e6;
    avaliaF(1,i) = nnz(power) / numel(power);
end
CF = sum(capacityF)/12;
AF = sum(avaliaF)/12;
disp("Average Capacity Factor: " + CF)
disp("Average Availability Factor: " + AF)
%% 3.a
disp("------------Expected Power Generation------------");

alpha = 0.638;
p = 3;
q = 1.5;
k = 1.96;
lambda = 9.13;

fx = @(x) wblpdf(x, lambda, k);
Fx = @(x) wblcdf(x, lambda, k);
Fxy = @(x, y) Fx(x).*Fx(y) .*(1+alpha*(1-Fx(x).^p).^q.*(1-Fx(y).^p).^q);
fxy = @(x, y) fx(x).*fx(y) .*(1+alpha*(1-Fx(x).^p).^(q-1).*(1-Fx(y).^p).^(q-1).*((Fx(x).^p).*(1+p*q)-1).*((Fx(y).^p).*(1+p*q)-1));

x = 0:0.5:30;
phi = @(x) P(x);
temp = phi(x) .* fx(x)';
[M,I] = max(temp);
mu = x(I);
sigma = M/1e5;
gx = @(x) normpdf(x,mu,sigma);
IS = @(x) (fx(x) .* P(x)') ./ gx(x);
v1 = sigma.*randn(1,N) + mu;
v2 = sigma.*randn(1,N) + mu;
z1 = IS(v1);
z2 = IS(v2);
tau1 = mean(z1);
tau2 = mean(z2);
tauSum = tau1 + tau2;
disp("Expected Value: " + tauSum)
%% 3.b
disp("------------Covariance of Produced Power------------");

x1 = 0:0.5:30;
x2 = 0:0.5:30;
phi = @(x) P(x);
temp = phi(x1) .* phi(x2) .* fxy(x1,x2)';
[M,I] = max(temp);
mu = [x(I), x(I)];
sigma = eye(2) * (M/1e10);

num = mvnrnd(mu,sigma,N);
X = num(:,1);
Y = num(:,2);
gx = mvnpdf([X Y], mu, sigma);
P1 = P(X);
P2 = P(Y);
P12Exp = mean((P1.*P2).*(fxy(X, Y)./gx));
cov = P12Exp - tau1 * tau2;
disp("Covariance: " + cov)
%% 3.C
disp("------------Variability & Standard Deviation------------");
speed = wblrnd(lambda, k, N, 1);
power = P(speed);
va = var(power);
vaxy = va * 2 + 2 * cov;
sdvaxy = sqrt(vaxy);
disp("Variance: " + vaxy);
disp("Standard Deviation: " + sdvaxy)
%% 3.d
disp("------------Confidence Interval------------");
% sigma0 = eye(2) * 55;
% gx = mvnpdf([x1',x1'],mu,sigma0);
% plot(x1, temp/1e10);
% hold on

sigma1 = eye(2)*(M/1e10); 
sigma2 = eye(2)*55;
N=10000;

X1 = mvnrnd(mu,sigma1,N);
g1 = mvnpdf(X1,mu,sigma1);
X2 = mvnrnd(mu,sigma2,N);
g2 = mvnpdf(X2,mu,sigma2);

P1 = P(X1(:,1));
P2 = P(X1(:,2));

P1_2 = P(X2(:,1));
P2_2 = P(X2(:,2));

tau_b = mean(((P1+P2)<9.5e6).*(fxy(X1(:,1),X1(:,2))./g1));
tau_a = mean(((P1_2+P2_2)>9.5e6).*(fxy(X2(:,1),X2(:,2))./g2));

above_U = tau_a + abs(se * std((P1_2+P2_2)>9.5e6)/sqN)
above_L = tau_a - abs(se * std((P1_2+P2_2)>9.5e6)/sqN)
above_W = above_U - above_L
below_U = tau_b + abs(se * std((P1+P2)<9.5e6)/sqN)
below_L = tau_b - abs(se * std((P1+P2)<9.5e6)/sqN)
below_W = below_U - below_L