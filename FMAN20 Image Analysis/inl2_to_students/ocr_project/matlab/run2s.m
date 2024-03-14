clear; clc; close;
im = imread('im1.jpg');
S = im2segment(im);
BW = S{1};
chk = [];
for c = 1:size(BW,2)
    chk(c) = any(BW(:,c));
end
%estabishing the exact column number wth 1 pixels present
corn = find(chk(1,:)==1);

MP = zeros(size(BW,1),size(corn,2));
MP = BW(:,(corn(1)-3):(corn(size(corn,2))+3));

p = regionprops(MP,'EulerNumber')
%o = (cat(1,p.EulerNumber))/100;

imshow(MP)