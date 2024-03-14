%10
load('population_2023.mat')
N = 1000;
n = 50; %number of observation
tau = zeros(1,n+1); %filter expectation
w = zeros(N,1); %weights 
p = @(x,y) unifpdf(y, 0.7*x, 1.2*x); %observation density, for weights
part = unifrnd(0.6, 0.99, [N 1]); %initialization
w = p(part,Y(1)); %weighting
tau(1) = sum(part.*w)/sum(w); %estimation

%confidence interval
[x_sort,index]=sort(part); % sort data
cum_weight=cumsum(w(index))/sum(w); %cumulative normalized sum of the weights
%for sorted data
indexLower=find(cum_weight>=0.025,1); % index for lower 2.5% quantile
indexUpper=find(cum_weight>=0.975,1); % index upper 2.5% quantile
tauLower(1)=x_sort(indexLower); % lower 2.5% quantile
tauUpper(1)=x_sort(indexUpper); % upper 2.5% quantile
%% 

ind = randsample(N,N,true,w); %selection
part = part(ind);
B = unifrnd(0.9, 3.9, [N,n+1]);  %stochastic repoduction rate
for k = 1:n % main loop
    part = B(:,k+1).* part.*(1-part); %mutation
    w = p(part,Y(k + 1)); %weighting
    tau(k + 1) = sum(part.*w)/sum(w); %estimation

    ind = randsample(N,N,true,w); %selection
    part = part(ind);
    
    %confidence interval
    [x_sort,index]=sort(part); % sort data
    cum_weight=cumsum(w(index))/sum(w); %cumulative normalized sum of the weights
    %for sorted data
    Ilower=find(cum_weight>=0.025,1); % index for lower 2.5% quantile
    Iupper=find(cum_weight>=0.975,1); % index upper 2.5% quantile
    tauLower(k+1)=x_sort(indexLower); % lower 2.5% quantile
    tauUpper(k+1)=x_sort(indexUpper); % upper 2.5% quantile
end
%% 

figure()
hold on
plot(1:n+1, tau, '--o','LineWidth', 1);
plot(1:n+1, X, '--x', 'LineWidth', 1);
hold off
xlabel('generation')
ylabel('relative population size')
legend('Estimation','X')
title('Filter estimate')

figure()
hold on
plot(1:n+1, X, '--x','LineWidth', 1);
plot(1:n+1, tauLower, '--','LineWidth', 1);
plot(1:n+1, tauUpper, '--','LineWidth', 1);
hold off
xlabel('generation')
ylabel('relative population size')
legend('X','lower bound','upper bound')
title('Confidence interval')