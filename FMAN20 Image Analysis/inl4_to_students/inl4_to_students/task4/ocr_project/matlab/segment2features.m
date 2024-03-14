function features = segment2features(I)
% features = segment2features(I)
%features = [];
BW = I;
[row,col] = find(BW);
rowmin = min(row);
rowmax = max(row);
colmin = min(col);
colmax = max(col);
si = BW(rowmin:rowmax,colmin-1:colmax+1);
si = imresize(si, [23, 18]);
m = size(si,1);
n = size(si,2);
for i=1:m 
    for j=1:n
        if si(i,j)>1
           si(i,j)=1;
        end
    end
end

fill = imfill(si,8,'holes');
hole = fill - si;
if hole == 0
 s(1) = 0;
 s(2) = 0;
else
o = ones(m,1)*[1:n];
p = [1:m]'*ones(1,n);
area = sum(sum(hole));
meanx = sum(sum(hole.*o))/area;
meany = sum(sum(hole.*p))/area;
s(1) = meanx/10;
s(2) = meany/10;
end
a1 = regionprops(si,'Perimeter');
%s(8) = (cat(1,a1.Perimeter))/100;
a2 = regionprops(si,'EulerNumber');
%s(2) = cat(1,a2.EulerNumber);
a3 = regionprops(si,'Eccentricity');
%s(3) = cat(1,a3.Eccentricity);
a4 = regionprops(si,'Extent');
%s(4) = cat(1,a4.Extent);
a5 = regionprops(si,'FilledArea');
%s(5) = (cat(1,a5.FilledArea))/100;
a6 = regionprops(si,'Centroid');
%s(6:7) = (cat(1,a6.Centroid))/10;
a7 = regionprops(si,'Extrema');
%s(3:18) = (cat(1,a7.Extrema))/10;
%features = s';
%[hog_8x8, vis8x8] = extractHOGFeatures(si,'CellSize',[8 8]);
%hogFeatureSize = length(hog_8x8);
%trainingFeatures = zeros(numImages,hogFeatureSize,'single');
trainingFeatures = extractHOGFeatures(si,'Cellsize',[4 4]);  
fre = double(trainingFeatures');
features = fre;
%features(433) = meanx/10;
%features(434) = meany/10;
end