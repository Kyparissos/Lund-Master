clear;clc;
close all
clear all
load(fullfile('databases','hep_proper_mask'));
X1_masks = Y1;
X2_masks = Y2;
load(fullfile('databases','hep_proper'));
%% 
i = 400;
I  = X1(:,:,1,i);
%% 
close all
Im = imread("test.png");
I = rgb2gray(Im);
% Bi = imbinarize(I,48);
I  = imgaussfilt(I,8);
% I = uint8(I);
ed = edge(I,"canny");
SE = strel('disk',3);
ed = imdilate(ed,SE);
CC.NumObjects = bwconncomp(ed,4);
J = grayconnected(I,100,200);
% J = grayconnected(I,150,240);
figure(1)
imagesc(labeloverlay(Im,J))

figure
imshow(ed)

% subplot(1,3,1);
% imshow(I,[0 255])
% 
% subplot(1,3,2);
% imshow(ed)

% figure(2)
% % subplot(1,3,3);
% imagesc(gau)
% colormap gray