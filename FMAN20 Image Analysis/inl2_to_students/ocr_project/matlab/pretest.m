clear; clc; close;
im = imread('im1.jpg');
S = im2segment(im);
BW = S{2};
processedImage = imbinarize(im2gray(BW));
figure;
subplot(1,2,1)
imshow(BW)
subplot(1,2,2)
imshow(processedImage)