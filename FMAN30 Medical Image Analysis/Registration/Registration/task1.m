clear;
HE = dir('Collection 1/HE/*.bmp');
p63AMACR = dir('Collection 1/p63AMACR/*.bmp');
 
for i=1:length(HE)
    im1{i} = imread(fullfile(HE(i).folder, HE(i).name));
    im1{i} = imresize(im1{i},0.5);
    im1{i} = rgb2gray(im1{i});
    im2{i} = rgb2gray(imresize(imread(fullfile(p63AMACR(i).folder, ...
        p63AMACR(i).name)),0.5));
end
%% 
close all
T1 = [];
for i= 1:length(im1)
% im1 = im1{i};
% im2 = im2{i};
points1 = detectSIFTFeatures(im1{i});
points2 = detectSIFTFeatures(im2{i});
[features1, valid_points1] = extractFeatures(im1{i}, points1);
[features2, valid_points2] = extractFeatures(im2{i}, points2);
indexPairs = matchFeatures(features1,features2, ...
    "MatchThreshold",8,"MaxRatio",0.6);
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);
matchedLocations1 = matchedPoints1.Location';
matchedLocations2 = matchedPoints2.Location';
x= matchedLocations1;
y= matchedLocations2;
% figure; 
% showMatchedFeatures(im1,im2,matchedPoints1,matchedPoints2);
% [R,t] = rigid_registration(x,y);
 
[R, t,percent] = myransac(x,y,3); 
% pp(i,3)= percent;
td = sqrt(t(1)^2+t(2)^2);
Rd = acosd(R(1,1));
s = 1;
T1 = cat(1,T1,[Rd,td]);
Rs{i} = R;
ts{i} = t;
T = [s*R, t ;0, 0, 1];
tform = affine2d(T');
im1_warp = imwarp(im1{i}, tform, 'OutputView', imref2d(size(im2{i})));
% figure;
% imshow(imfuse(im2{i}, im1_warp, 'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
end
%% Task 2
clear;
part2 = dir("Collection 2\HE\*.jpg");
part1 = dir("Collection 2\TRF\*.tif");

for i = 1:length(part1)
    im1{i} = im2gray(imresize(imread(fullfile(part1(i).folder, ...
        part1(i).name)),0.5));
    im2{i} = im2gray(imresize(imread(fullfile(part2(i).folder, ...
        part2(i).name)),0.5));
    im1{i} = imcomplement(histeq(im1{i}));
end

%% 
close all
T1 = [];
for i= 1:length(im1)
% im1 = im1{i};
% im2 = im2{i};
points1 = detectSIFTFeatures(im1{i});
points2 = detectSIFTFeatures(im2{i});
[features1, valid_points1] = extractFeatures(im1{i}, points1);
[features2, valid_points2] = extractFeatures(im2{i}, points2);
indexPairs = matchFeatures(features1,features2, ...
    "MatchThreshold",8,"MaxRatio",0.7);

matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);
matchedLocations1 = matchedPoints1.Location';
matchedLocations2 = matchedPoints2.Location';
x= matchedLocations1;
y= matchedLocations2;
% figure; 
% showMatchedFeatures(im1,im2,matchedPoints1,matchedPoints2);

[R, t, s,percent] = myransac2(x,y,4);
pp(i,1)= percent;
td = sqrt(t(1)^2+t(2)^2);
Rd = acosd(R(1,1));
T1 = cat(1,T1,[Rd,td,s]);
T = [s*R, t ;0, 0, 1];
tform = affine2d(T');
im1_warp = imwarp(im1{i}, tform, 'OutputView', imref2d(size(im2{i})));
figure; 
% imshow(imfuse(im2{i}, im1_warp, 'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
imshow(imfuse(im2{i}, im1_warp, "blend"));

end

