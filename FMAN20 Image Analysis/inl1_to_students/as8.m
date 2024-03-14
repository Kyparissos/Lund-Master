clear;clc;
load("assignment1bases.mat")
for k = 1:2
    for i = 1:3
        for j = 1:400
            r(i,j) = project(stacks,bases,k,i,j);
        end
        rm(i,k) = mean(r(i,:));
       
    end
end
rm

for k = 1:2
    for j = 1:15
    figure;
    imagesc(stacks{k}(:,:,j))
    colormap(gray)
    end
end
%{
for i = 1:3
    for p = 1:4
        figure;
        colormap(gray)
        imagesc(bases{i}(:,:,p))
    end
end
%}
