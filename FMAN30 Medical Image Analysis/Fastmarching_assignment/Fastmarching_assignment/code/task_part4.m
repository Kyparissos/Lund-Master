clear;
%% Loading MR-carotid-coronal
path = 'C:\Users\jingm\Desktop\Fastmarching_assignment\Fastmarching_assignment\images\MR-carotid-coronal';
[im, info] = mydicomreadfolder(path);
%% 
imagesc(im(:,:,30))
colormap gray
%% fastmarch
rows = size(im,1);
columns = size(im,2);
depth = size(im,3);
left_seed = [290,130,30];
right_seed = [275,368,30];
seedpoint = right_seed;
startind = sub2ind([rows columns depth],seedpoint(1),seedpoint(2),seedpoint(3));
sigma = 2500;
SPEED = exp(-(im - im(seedpoint(1),seedpoint(2),seedpoint(3))).^2/(2*sigma^2));
SPEED = exp(-(im - im(seedpoint(1),seedpoint(2),seedpoint(3))).^2/sigma);

im_cost = 1./SPEED;
maxendcost = 1e9;
arrivalmap = fastmarch(single(im_cost),single(startind),maxendcost);
% imshow3D(arrivalmap)
% figure()
% imagesc(arrivalmap(:,:,31))
% colorbar
%% 3d
[x,y,z] = meshgrid(1:rows,1:columns,1:depth);
thresholds = 1e6;
isosurface(x,y,z,arrivalmap<=thresholds);
cameratoolbar
lighting flat