function S = im2segment(im)
%applying gaussian filter
%im = imread('im8.jpg')



%ost = 0
imgauss = imgaussfilt(im);
m = mean(imgauss);
mm = mean(m);
mcol = mean(imgauss);
mrow = mean(imgauss,2);
mrow = mrow';
ks = double(mcol > 22);
ksj = double(mrow > 11);
symd = find(ks(1,:)==1);
syme = find(ksj(1,:)==1);
olool = symd(end);
olooh = syme(end);
imgauss(:,olool:end) = 0;
imgauss(olooh:end,:) = 0;
ks2 = double(mcol > 17);
symd2 = find(ks2(1,:)==1);
for out = 1:(length(symd2)-1)
    kls = symd2(out+1) - symd2(out);
    rko = round(length(im) * 0.035);
    if kls > rko
        jj = symd2(out+1) - 2;
        ksaa = symd2(out) + 2;
        imgauss(:, ksaa:jj) = 0;
    end
end




ss = imgauss;
chk = (zeros(1,size(ss,2)));%blank matrix
%establishing a critical value for thresholding
%two methods were used based on the average pixel intensity
if mm <= 20
    critic = mm*4;
elseif mm <= 33 && mm > 20
    critic = mm*1.8;
else
    critic = 41;
end
%thresholding process
for i=1:size(imgauss,1)
    for j=1:size(imgauss,2)
        if imgauss(i,j) >= critic
            ss(i,j) = 1;
        else
            ss(i,j) = 0;
        end
    end
end
%deleting any bright pixels on the border
%relevant for next step
ss(1,:) = 0;
ss(end,:) = 0;
ss(:,1) = 0;
ss(:,end) = 0;
%reducing 'noise' pixels without any 4-neighbors will be removed
for nn=2:(size(ss,1)-1)
    for mm=2:(size(ss,2)-1)
        if ss(nn+1,mm)== 0 && ss(nn-1,mm) == 0 && ss(nn,mm+1) == 0 && ss(nn,mm-1) == 0
            ss(nn,mm) = 0;
        end
    end
end
%check to see which columns have '1' pixels present
for c = 1:size(ss,2)
    chk(c) = any(ss(:,c));
end
%estabishing the exact column number wth 1 pixels present
symb = find(chk(1,:)==1);
lisp = 0; %counter
syma = find(chk(1,:)==1);%size check
symb(size(symb,2)+1) = 9999999;%increasing size for looping syntax purposes

%determining the number of segments
for k = 1:(size(syma,2))
    kmj = symb(k+1)-symb(k);
    if kmj >= 8
        lisp = lisp+1;
    end
end

S = cell(1,lisp); %creating the cell array space
KJJ = zeros(size(imgauss)); %blank matrix size of image to fill with segments
LL = 1;%counter
kmj= [];%variable to signify segment break
%segmentation process
for k = 1:(size(syma,2))
    kmj = symb(k+1)-symb(k);
    ff = round(length(im) * 0.025);
    if kmj <= ff
        KJJ(:,symb(k)) = ss(:,symb(k));
    else
        KJJ(:,symb(k)) = ss(:,symb(k));
        S{LL} = KJJ;
        LL = LL + 1;
        KJJ = zeros(size(imgauss));
    end
end
for kjl = 1:length(S)
    if ((size(S{kjl},1) == 0) && (size(S{kjl},2) == 0))
        S(kjl) = [];
    end
end
end
