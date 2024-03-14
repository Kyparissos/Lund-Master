function S = im2segment(im)
llkl = imgaussfilt(im)
vf = llkl;
avg = mean(vf)
avgg = mean(avg)
ss = vf;
chk = (zeros(1,size(ss,2)));
if avgg <= 23
    critic = avgg*1.3
else
    critic = 41
end

for i=1:size(vf,1)
    for j=1:size(vf,2)
        if vf(i,j) >= critic
            ss(i,j) = 1;
        else
            ss(i,j) = 0;
        end
    end
end

ss(1,:) = 0;
ss(end,:) = 0;
ss(:,1) = 0;
ss(:,end) = 0;

for nn=2:(size(ss,1)-1)
    for mm=2:(size(ss,2)-1)
        if ss(nn+1,mm)== 0 && ss(nn-1,mm) == 0 && ss(nn,mm+1) == 0 && ss(nn,mm-1) == 0
            ss(nn,mm) = 0;
        end
    end
end

for c = 1:size(ss,2)
    chk(c) = any(ss(:,c));
end

symb = find(chk(1,:)==1);
lisp = 0
syma = find(chk(1,:)==1);
symb(size(symb,2)+1) = 99999999;
for k = 1:(size(syma,2))
    kmj = symb(k+1)-symb(k);
    if kmj >= 8
        lisp = lisp+1
    end
end

S = cell(1,lisp);
KJJ = zeros(size(vf));
LL = 1;
kmj= [];

for k = 1:(size(syma,2))
    kmj = symb(k+1)-symb(k);
    if kmj <= 8
        KJJ(:,symb(k)) = ss(:,symb(k));
    else
        KJJ(:,symb(k)) = ss(:,symb(k))
        S{LL} = KJJ
        LL = LL + 1;
        KJJ = zeros(size(vf));
    end
end
end
