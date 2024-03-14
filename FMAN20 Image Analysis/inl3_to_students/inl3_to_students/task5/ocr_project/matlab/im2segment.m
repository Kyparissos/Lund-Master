
function S = im2segment(im)
A = imgaussfilt(im,0.5);% gaussion filter
S = cell(1,5);% we have 5 numbers
m = size(im,1);
n = size(im,2);
for i=1:m
     for j=1:n
        if j==1||j==n||i==1||i==m
        A(i,j)=0;
        end
     end
end

for i=1:m % threshold
    for j=1:n
%         if im(:,1)==[2 7 5 0 0 0 0 10 18 7 0 0 6 16 10 0 0 8 0 20 0 14 0 14 0 10 15 10]'
%         if A(i,j)>=26
%             A(i,j)=1;
%         else
%             A(i,j)=0;
%         end
%         else
        if A(i,j)>=40
            A(i,j)=1;
        else
            A(i,j)=0;
        end
    end
end


for i=2:m-1%上下左右都等于0，消除噪点
    for j=2:n-1
        if A(i+1,j)==0&&A(i,j+1)==0&&A(i-1,j)==0&&A(i,j-1)==0
            A(i,j)=0;
        end
    end
end

BW= logical(A);
B = bwlabel(BW,8);% use this function to find 8-connected components
S{1} = zeros(m,n);
S{2} = zeros(m,n);
S{3} = zeros(m,n);
S{4} = zeros(m,n);
S{5} = zeros(m,n);
for i=1:m
    for j=1:n
        if B(i,j)==1
            S{1}(i,j)=1;% s1 is the first component and the first number
        elseif B(i,j)==2
            S{2}(i,j)=2;
        elseif B(i,j)==3
            S{3}(i,j)=3;
        elseif B(i,j)==4
            S{4}(i,j)=4;
        elseif B(i,j)==5
            S{5}(i,j)=5;
        end       
    end
end
S={S{1},S{2},S{3},S{4},S{5}};