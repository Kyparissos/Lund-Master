load sysid02.mat
N = length(y);
Ntrain = round(N/2);
ytrain = y(1:Ntrain);
ytest = y(Ntrain+1:N);
utrain = u(1:Ntrain);
utest = u(Ntrain+1:N);

ztrain = iddata(ytrain,utrain);
ztest = iddata(ytest,utest);

% plot and check that there are no outliers, and that signals have zero
% mean, otherwise fix outliers and remove means

NN = struc(1:10,1:10,1:10);
V = arxstruc(ztrain, ztest,NN),

Nbest = selstruc(V,0)
arxbest = arx(ztrain,Nbest)

oemodel = oe(ztrain,[4,4,1])

figure(1)
compare(ztest,arxbest,oemodel)
print -depsc prob2_01.eps

figure(2)

p = bodeoptions;
set(p,'ConfidenceRegionNumberSD',5)
h = bodeplot(arxbest,oemodel,p)
showConfidence(h);
print -depsc prob2_02.eps

%% 
figure(3)
pzmap(oemodel)
print -depsc prob2_03.eps
%% 

figure(4)
resid(ztest, arxbest,'b',oemodel,'r')
legend('arxbest','oemodel')
print -depsc prob2_04.eps
