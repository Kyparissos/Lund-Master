clear
%% 
% Gaussian mean and standard devition
m1 = 0;
m2 = 1;
% s1= 0.4;
s1 = 4;

% Generate traning data and testing data
train1 = normrnd(m1,s1,10,1);
train2 = normrnd(m2,s1,10,1);
train = cat(1,train1,train2);
test1 = normrnd(m1,s1,1000,1);
test2 = normrnd(m2,s1,1000,1);
testdata = cat(1,test1,test2);
% Parzen window estimation
width1 = 0.02;
width2 = 1;

x = linspace(-20,20,1000);
x = x';
%% 
% pdf1 = kde(train,Bandwidth= widths);
[pdf11] = ksdensity(train1,x,'Bandwidth',width1);
[pdf21] = ksdensity(train2,x,'Bandwidth',width1);
[pdfw1] = ksdensity(train,x,'Bandwidth',width1);

[pdf12] = ksdensity(train1,x,'Bandwidth',width2);
[pdf22] = ksdensity(train2,x,'Bandwidth',width2);
[pdfw2] = ksdensity(train,x,'Bandwidth',width2);

%% width1
figure(1)
plot(x,pdf11)
hold on
plot(x,pdf21)
plot(x,pdfw1)
plot(train1,0,'bo')
plot(train2,0,'ro')
hold off
legend('p(x|y=1)','p(x|y=2)','p(x)')
%% width2
figure(2)
plot(x,pdf12)
hold on
plot(x,pdf22)
plot(x,pdfw2)
plot(train1,0,'bo')
plot(train2,0,'ro')
hold off
legend('p(x|y=1)','p(x|y=2)','p(x)')
%% true model
figure(3)
y1true = normpdf(x,m1,s1);
hold on
y2true = normpdf(x,m2,s1);
ytrue = y1true + y2true;
plot(x,normpdf(x,m1,s1))
plot(x,normpdf(x,m2,s1))
plot(x,ytrue)
legend('y=1 true','y=2 true','y true')
%% 
y1c = normcdf(x,m1,s1);
y2c = normcdf(x,m2,s1);
plot(x,y1c)
hold on
plot(x,y2c)
%% width1
% py1|x = px|y1 * py1/ px ;
py1x = pdf11.*0.5./pdfw1;
py2x = pdf21.*0.5./pdfw1;
figure(4)
plot(x,py1x)
hold on
plot(x,py2x)
hold off
legend('p(y=1|x)','p(y=2|x)')
%% width2
py1x = pdf12.*0.5./pdfw2;
py2x = pdf22.*0.5./pdfw2;
figure(4)
plot(x,py1x)
hold on
plot(x,py2x)
hold off
legend('p(y=1|x)','p(y=2|x)')
%% true model
py1x = y1true.*0.5./ytrue;
py2x = y2true.*0.5./ytrue;
figure(4)
plot(x,py1x)
hold on
plot(x,py2x)
hold off
legend('p(y=1|x)','p(y=2|x)')
%% find the threshold
% for i = 1:100
%     if py1x(i)==py2x(i)
%         threshold = x(i);
%     elseif py1x(i)<py2x(i)
%             threshold = x(i-1);
%             break
%     end
% end
%% plug-in
ctright = 0;
ctwrong = 0;
for i = 1:1000
    for j = 1:1000
        if test1(i)>x(j)&&test1(i)<x(j+1)
            pt = j;
            break
        end
    end
    if py1x(pt)>=py2x(pt)
        ctright = ctright+1;
    else
        ctwrong = ctwrong+1;
    end
end
%% 
ctright = 0;
ctwrong = 0;
for i = 1:1000
    for j = 1:1000
        if test2(i)>x(j)&&test2(i)<x(j+1)
            pt = j;
            break
        end
    end
    if py1x(pt)<=py2x(pt)
        ctright = ctright+1;
    else
        ctwrong = ctwrong+1;
    end
end
%% cross validation
sumloglikelihood = 0;
partition = cvpartition(20,"KFold",5);
%% 

for k = 1:5
    loglikelihood = 0;
    idxTrain = training(partition,k);
    cvTrain = train(idxTrain,:);
    idxTest = test(partition,k);
    cvTest = train(idxTest,:);
    [pdf] = parzen(cvTrain,x);
    for i = 1:length(cvTest)
        for j = 1:1000
            if cvTest(i)>x(j)&&cvTest(i)<x(j+1)
                pt = j;
                break
            end
        end
        likelihood = -log10(pdf(pt));
        loglikelihood = loglikelihood+likelihood;
    end
    sumloglikelihood = sumloglikelihood + loglikelihood;
end
ht = matlabFunction(sumloglikelihood);
h = fminbnd(ht,0.01,1)
%% 
sumloglikelihood = 0;
partition = cvpartition(20,"KFold",5);

maxloglikelihood = zeros(100,1);
hs = linspace(0.01,5,100);
for s = 1:100    
    sumloglikelihood = 0;
    for k = 1:5
        logp = 0;
        idxTrain = training(partition,k);
        cvTrain = train(idxTrain,:);
        idxTest = test(partition,k);
        cvTest = train(idxTest,:);
        pdf = my_parzen(x,cvTrain,hs(s));
        for i = 1:length(cvTest)
            for j = 1:1000
                if cvTest(i)>x(j)&&cvTest(i)<x(j+1)
                    pt = j;
                    break
                end
            end
            p = log10(pdf(pt));
            logp = logp + p;
        end
        sumloglikelihood = sumloglikelihood + logp;
    end
    maxloglikelihood(s,1) = sumloglikelihood;
end
[M,h] = max(maxloglikelihood);
optimalwidth = hs(h)
