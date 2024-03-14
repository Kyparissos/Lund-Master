clear;clc;
im = imread('im8.jpg');
S = im2segment(im);
BW = S{5};
%s = [];
[row,col] = find(BW);
rowmin = min(row);
rowmax = max(row);
colmin = min(col);
colmax = max(col);
si = BW(rowmin-1:rowmax+1,colmin-3:colmax+3);
%si = imresize(si, [23, 23]);
m = size(si,1);
n = size(si,2);
for i=1:m 
    for j=1:n
        if si(i,j)>1
           si(i,j)=1;
        end
    end
end
a1 = regionprops(si,'Perimeter');
s(1) = (cat(1,a1.Perimeter))/100;
a2 = regionprops(si,'EulerNumber');
s(2) = cat(1,a2.EulerNumber);
a3 = regionprops(si,'Eccentricity');
s(3) = cat(1,a3.Eccentricity);
a4 = regionprops(si,'Extent');
s(4) = cat(1,a4.Extent);
a5 = regionprops(si,'FilledArea');
s(5) = (cat(1,a5.FilledArea))/100;
a6 = regionprops(si,'Centroid');
s(6:7) = (cat(1,a6.Centroid))/10;
a7 = regionprops(si,'Extrema');
s(8:23) = (cat(1,a7.Extrema))/10;

%features = s'
%[hog_8x8, vis8x8] = extractHOGFeatures(si,'CellSize',[8 8]);
cellSize = [4 4];
%hogFeatureSize = length(hog_8x8);
%trainingFeatures = zeros(si,hogFeatureSize,'single');
trainingFeatures = extractHOGFeatures(si,'CellSize',cellSize);  
fre = double(trainingFeatures');
features = fre;
imshow(si)