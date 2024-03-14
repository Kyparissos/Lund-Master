% clear;clc;
% I = beach1;
% In = beach2;
% I(1:size(In,1),size(I,2)+1:size(I,2)+size(In,2),:) = In;
% In = beach3;
% I(1:size(In,1),size(I,2)+1:size(I,2)+size(In,2),:) = In;
% In = beach4;
% I(1:size(In,1),size(I,2)+1:size(I,2)+size(In,2),:) = In;
% In = beach5;
% I(1:size(In,1),size(I,2)+1:size(I,2)+size(In,2),:) = In;
% In = beach6;
% I(1:size(In,1),size(I,2)+1:size(I,2)+size(In,2),:) = In;
% imshow(I);
%% Q1.1 RGB
load('beachall.mat')
I = double(I);
A = reshape(I,194*1667,3);
% redChannel=I(:, :, 1);
% greenChannel=I(:, :, 2);
% blueChannel=I(:, :, 3);
% A=double([redChannel(:), greenChannel(:), blueChannel(:)]);
k = 10;
[mc,n] = kmeans(A,k);
m=reshape(mc,size(I,1),size(I,2));
n=n/255;
clusteredImage=label2rgb(m,n);
imshow(clusteredImage)
% classim = reshape(classification_data,194,[],1);
% x = label2rgb(classim);
% imshow(x)
new = imread("newbeach.jpg");
B = double(new);
C = reshape(B,720*1080,3);
nnmodel = fitcknn(A,mc);
predict = nnmodel.predict(C);
y = reshape(predict,720,[],1);
y = label2rgb(y,n);
%imshow(y)
%% Q1.2 CIE LAB
load('beachall.mat')
I = double(I);
Ilab = rgb2lab(I);
A = reshape(Ilab,194*1667,3);
k = 10;
[mc,n] = kmeans(A,k);
m=reshape(mc,size(Ilab,1),size(Ilab,2));
n=n/100;
for i=1:size(m,1)
    for j=1:size(m,2)
        num = m(i,j,1);
        m(i,j,1) = n(num,1);
        m(i,j,2) = n(num,2);
        m(i,j,3) = n(num,3);
    end
end
clusteredImage=lab2rgb(m);
imshow(clusteredImage)
%%
new = imread("newbeach.jpg");
newd = double(new);
P = rgb2lab(new);
B = rgb2lab(newd);
C = reshape(B,720*1080,3);
nnmodel = fitcknn(A,mc);
predict = nnmodel.predict(C);
y = reshape(predict,720,[],1);
for i=1:size(y,1)
    for j=1:size(y,2)
        num = y(i,j,1);
        y(i,j,1) = n(num,1); 
        %P(i,j,1);
        y(i,j,2) = n(num,2);
        y(i,j,3) = n(num,3);
    end
end
y1 = lab2rgb(y);
imshow(y1)

%% Q4
clear;clc;
Irgb = imread("sunset.jpg");
I = im2gray(Irgb);
[mserRegions, mserConnComp] = detectMSERFeatures(I, ... 
    'RegionAreaRange',[5000 12000],'ThresholdDelta',4);

figure
imshow(I)
hold on
plot(mserRegions, 'showPixelList', true,'showEllipses',false)
title('MSER regions')
hold off
%% red channel
% I = imread("sunset.jpg");
% redChannel=I(:, :, 1);
% for i=1:size(I,1)
%     for j=1:size(I,2)
%         if  redChannel(i,j)>235
%             In(i,j)=1;
%         else
%             In(i,j)=0;
%         end
%     end
% end
% imshow(In)
%% lab channel
I = imread("sunset.jpg");
%I = double(I);
Ilab = rgb2lab(I);
aChannel=Ilab(:, :, 2);
for i=1:size(I,1)
    for j=1:size(I,2)
        if  aChannel(i,j)>40
            In(i,j)=1;
        else
            In(i,j)=0;
        end
    end
end
imshow(In)
%% 
mask = logical(In);
% se=strel('disk', 4); 
% mask=imdilate(mask,se);
mask= imgaussfilt(double(mask),1);
J = inpaintCoherent(I,logical(mask));
%montage({I,J});
imshow(J)
%% 
% Inew = I;
% for i=1:size(I,1)
%     for j=1:size(I,2)
%         if  mask(i,j)>0
%             Inew(i,j,:)=I(i-260,j,:);
%         end
%     end
% end
%imshow(Inew)
%% Q3.1
r = imread("butterfly_color.jpg");
T = imread("butterfly_blobs.jpg");
r = double(r)/255;
T = double(T)/255;
% redchannelco = double(r(:,:,1))/255;
% greenchannelco = double(r(:,:,2))/255;
% bluechannelco = double(r(:,:,3))/255;
% redchannelbl = double(T(:,:,1))/255;
% greenchannelbl = double(T(:,:,2))/255;
% bluechannelbl = double(T(:,:,3))/255;
%b=[b1;b2;b3];
% r=(redchannelco;greenchannelco;bluechannelco);
% T=(redchannelbl;greenchannelbl;bluechannelbl);
% xm = (r(:,:,1)-r(:,:,1).^2);
% xm = xm(:);
% ym = (T(:,:,1)-r(:,:,1).^2);
% ym = ym(:);
% b1 = xm\ym;
% xm = (r(:,:,2)-r(:,:,2).^2);
% xm = xm(:);
% ym = (T(:,:,2)-r(:,:,2).^2);
% ym = ym(:);
% b2 = xm\ym;
% xm = (r(:,:,3)-r(:,:,3).^2);
% xm = xm(:);
% ym = (T(:,:,3)-r(:,:,3).^2);
% ym = ym(:);
% b3 = xm\ym;
for i= 1:3
    xm = (r(:,:,i)-r(:,:,i).^2);
    xm = xm(:);
    ym = (T(:,:,i)-r(:,:,i).^2);
    ym = ym(:);
    b(i) = xm\ym;
    new(:,:,i)=b(i)*r(:,:,i)+(1-b(i))*r(:,:,i).^2;
end
% b1 = (r(:,:,1)-r(:,:,1).^2).\(T(:,:,1)-r(:,:,1).^2);
% b2 = (r(:,:,2)-r(:,:,2).^2).\(T(:,:,2)-r(:,:,2).^2);
% b3 = (r(:,:,3)-r(:,:,3).^2).\(T(:,:,3)-r(:,:,3).^2);
% new(:,:,1)=b1*r(:,:,1)+(1-b1)*r(:,:,1).^2;
% new(:,:,2)=b2*r(:,:,2)+(1-b2)*r(:,:,2).^2;
% new(:,:,3)=b3*r(:,:,3)+(1-b3)*r(:,:,3).^2;
imshow(new)
%% Q3.2 ransac wrong
clear;clc;
r = imread("butterfly_color.jpg");
T = imread("butterfly_blobs.jpg");
r = double(r)/255;
T = double(T)/255;
for i=1:3
    xm = (r(:,:,i)-r(:,:,i).^2);
    xm = xm(:);
    ym = (T(:,:,i)-r(:,:,i).^2);
    ym = ym(:);
    sampleSize = 2; % number of points to sample per trial
    maxDistance = 2; % max allowable distance for inliers
%     points = [(1:762048)',ym];
    points = [xm,ym];
    fitLineFcn = @(xm,ym) (xm\ym);
    evalLineFcn = ...   % distance evaluation function
        @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);

    [modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn,sampleSize,maxDistance);
    % modelInliers = polyfit(points(inlierIdx,1),points(inlierIdx,2),1);
%     inlierPts = inlierIdx.*ym;
    b(i)=inlierIdx\ym;
    new(:,:,i)=b(i)*r(:,:,i)+(1-b(i))*r(:,:,i).^2;
end
imshow(new)
%% outlier wrong
clear;clc;
r = imread("butterfly_color.jpg");
T = imread("butterfly_blobs.jpg");
r = double(r)/255;
T = double(T)/255;
for i=1:3
    xm = (r(:,:,i)-r(:,:,i).^2);
    xm = xm(:);
    ym = (T(:,:,i)-r(:,:,i).^2);
    ym = ym(:);
    TF = isoutlier(ym);
    inlierIdx = (~TF);
    k = find(TF);
    inlierPts = [xm.*inlierIdx,ym.*inlierIdx];
    b(i)=inlierPts(:,1)\inlierPts(:,2);
    new(:,:,i)=b(i)*r(:,:,i)+(1-b(i))*r(:,:,i).^2;
end
imshow(new)
%% wrong
clear;clc;
r = imread("butterfly_color.jpg");
T = imread("butterfly_blobs.jpg");
r = double(r)/255;
Td = double(T)/255;
Td = Td(:);
TF = isoutlier(Td);
inlierIdx = (~TF);
k = find(TF);
T = inlierIdx.*Td;
for i=1:3
    xm = (r(:,:,i)-r(:,:,i).^2);
    xm = xm(:);
    ym = (T(:,:,i)-r(:,:,i).^2);
    ym = ym(:);       
    b(i)=xm\ym;
    new(:,:,i)=b(i)*r(:,:,i)+(1-b(i))*r(:,:,i).^2;
end
imshow(new)
%% threshold wrong
clear;clc;
T = imread("butterfly_blobs.jpg");
green = T(:,:,2);
for i=1:size(T,1)
    for j=1:size(T,2)
        if  green(i,j)<100
            T(i,j,:)=0;
        end
    end
end
imshow(T)
%% ransac right
clear;clc;
r = imread("butterfly_color.jpg");
T = imread("butterfly_blobs.jpg");
r = double(r)/255;
T = double(T)/255;
for i=1:3
    xm = (r(:,:,i)-r(:,:,i).^2);
    xm = xm(:);
    ym = (T(:,:,i)-r(:,:,i).^2);
    ym = ym(:);
    bestins = 0;
    bnd = 0.1;
    iters = 100;
    nn = length(xm);
    for ii = 1:iters
        id = randi(nn);
        a_test = xm(id)\ym(id);        
        ins = sum(abs(ym-(a_test*xm+(1-a_test)*xm.^2))<bnd);
        if ins>bestins
            bestins = ins;
            a_rs = a_test;            
        end
    end
    b(i)=a_rs;
    % best inlier percentage
    disp(bestins/nn)
    new(:,:,i)=b(i)*r(:,:,i)+(1-b(i))*r(:,:,i).^2;
end
imshow(new);
%% Q6
clear;clc;
load('mite_gt.mat')
I1=imread("mite_1.jpg");
I2=imread("mite_2.jpg");
I3=imread("mite_3.jpg");
I4=imread("mite_4.jpg");
I5=imread("mite_5.jpg");
I6=imread("mite_6.jpg");
I7=imread("mite_7.jpg");
I8=imread("mite_8.jpg");
I9=imread("mite_9.jpg");
I10=imread("mite_10.jpg");
I11=imread("mite_11.jpg");
I12=imread("mite_12.jpg");
%I = [I1,I2,I3,I4,I5,I6];
x1=X{1};
m1=insertMarker(I1,x1);
x2=X{2};
m2=insertMarker(I2,x2);
x3=X{3};
m3=insertMarker(I3,x3);
x4=X{4};
m4=insertMarker(I4,x4);
x5=X{5};
m5=insertMarker(I5,x5);
x6=X{6};
m6=insertMarker(I6,x6);
% I = [m1,m2,m3;
%     m4,m5,m6];
% figure(1)
% montage({m1,m2,m3,m4,m5,m6})
I = [I7,I8,I9;
    I10,I11,I12];
figure(1)
imshow(I)
%%
g1 = rgb2gray(I);
BW = imbinarize(g1,0.4);
bw = (~BW);
% for i = 1:size(I,1)
%     for j = 1:size(I,2)
%         if I(i,j)<160
%             bw(i,j) = 0;
%         else
%             bw(i,j) = 1;
%         end
%     end
% end
%bw = imfill(bw,'holes');
figure(2)
imshow(bw)
 a1 = imfindcircles(bw,[15 50],"Sensitivity",0.85,"EdgeThreshold",0.4);
 b1 = insertMarker(I,a1);
figure(3)
imshow(b1)
%% Q5
clear;clc;
I = imread("color_tomb.png");
I1 = I(1:1800,:,:);
I2 = I(1801:2400,:,:);
I1g = rgb2gray(I1);
I1bi = imbinarize(I1g);
I2g = rgb2gray(I2);
% color1 = I2(:475,120:150,:);
% mean1 = mean(color1,[1,2]);
% color2 = I2(501:530,120:150,:);
for k=1:6
    color(:,:,:,k)=I2((501-25*(k-1)):(510-25*(k-1)),121:130,:);
    meanx(:,:,:,k)=mean(color(:,:,:,k));
    meany(:,:,k)=mean(meanx(:,:,:,k));
end
% figure(1)
% imshow(I1g)
% figure(2)
% imshow(color1)
%I = double(I);
L = bwlabel(I1bi,4);
% max = max(L,[],'all');
for i=1:size(I1,1)
    for j=1:size(I1,2)
        for k=1:6
        if L(i,j)==k+1
            I(i,j,1)=meany(1,:,k);
            I(i,j,2)=meany(2,:,k);
            I(i,j,3)=meany(3,:,k);
        end
        end
    end
end
figure(3)
imshow(I)
%% 
I = I2;
Ibw = ~im2bw(I,graythresh(I));
Ifill = imfill(Ibw,'holes');
Iarea = bwareaopen(Ifill,100);
Ifinal = bwlabel(Iarea);
stat = regionprops(Ifinal,'boundingbox');
imshow(I); hold on;
for cnt = 1 : numel(stat)
    bb = stat(cnt).BoundingBox;
    rectangle('position',bb,'edgecolor','r','linewidth',2);
end
%% 
color(:,:,:,2)=I2((500-55*2):(530-55*2),120:150,:);
mean(color,[1,2]);





