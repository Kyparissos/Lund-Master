clear;clc;
I = imread("arcimboldo_low.jpg");
HSV = rgb2hsv(I);
m = size(HSV,1);
n = size(HSV,2);
HSV(:, :, 2) = HSV(:, :, 2) * 2;
HSV(HSV > 1) = 1;
new = hsv2rgb(HSV);
imshow(new)
%% 
I = imread("michelangelo_colorshift.jpg");
grayImage = rgb2gray(I);
redchannel = I(:,:,1);
greenchannel = I(:,:,2);
bluechannel = I(:,:,3);
meanr = mean2(redchannel);
meang = mean2(greenchannel);
meanb = mean2(bluechannel);
meanGray = mean2(grayImage);
redchannel = uint8(double(redchannel)* meanGray / meanr);
greenchannel = uint8(double(greenchannel)* meanGray / meang);
bluechannel = uint8(double(bluechannel)* meanGray / meanb);
rgbImage = cat(3, redchannel, greenchannel, bluechannel);
imshow(rgbImage)