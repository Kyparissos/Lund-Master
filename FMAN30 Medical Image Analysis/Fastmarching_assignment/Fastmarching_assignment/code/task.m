clear;
%% 
heart_path = 'C:\Users\jingm\Desktop\Fastmarching_assignment\Fastmarching_assignment\images\MR-heart-single.dcm';
thorax_path = 'C:\Users\jingm\Desktop\Fastmarching_assignment\Fastmarching_assignment\images\CT-thorax-single.dcm';
%% 
[info_heart_test] = mydicominfo(heart_path);
[info_thorax_test] = mydicominfo(thorax_path);

%% 
[info_heart,im_heart] = mydicomread(heart_path);
[info_thorax,im_thorax] = mydicomread(thorax_path);
% images become zeros since str2double converts scalar to NaN...
%% 
imagesc(im_heart)
colormap gray
colorbar
title('Heart image')
axis off
%% 
imagesc(im_thorax)
colormap gray
colorbar
title('Thorax image')
axis off
%% 2.1
clear;clc

transversal_path = 'C:\Users\jingm\Desktop\Fastmarching_assignment\Fastmarching_assignment\images\MR-thorax-transversal';
[im, info] = mydicomreadfolder(transversal_path);

spacing = strsplit(info.PixelSpacing, '\');
[sx, sy, sz] = size(im);
depth = [info.Rows*str2double(spacing(1)), 0];
width = [0, info.Columns*str2double(spacing{2})];
height = [0, int8(str2double(info.SpacingBetweenSlices)*(sz-1) + str2double(info.SliceThickness)*sz)];

% Transverse plane (x=width, y=depth)
figure(1)

im_transverse = im(:,:,sz/2);
im_transverse = reshape(im_transverse, [sx, sy]);
% im_transverse = transpose(im_transverse);
im_transverse = flip(im_transverse);
imagesc(width, depth, im_transverse)
colormap gray
colorbar
set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0]);
set(gca,'ticklength',[0.05 0.05]);
axis image
title('Transverse Plane')
xlabel('width [mm]')
ylabel('depth [mm]')

% Coronal plane (x=width, y=height)
figure(2)

im_coronal = im(sx/2,:,:);
im_coronal = reshape(im_coronal, [sy, sz]);
im_coronal = transpose(im_coronal);
% im_coronal = flip(im_coronal);
imagesc(width, height, im_coronal)
colormap gray
colorbar
set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0]);
set(gca,'ticklength',[0.05 0.05]);
axis image
title('Coronal Plane')
xlabel('width [mm]')
ylabel('height [mm]')

% Sagital plane (x=depth, y=height)
figure(3)

im_sagital = im(:,sy/2,:);
im_sagital = reshape(im_sagital, [sx, sz]);
im_sagital = transpose(im_sagital);
im_sagital = flip(im_sagital,2);
imagesc(depth, height, im_sagital)
colormap gray
colorbar
set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0]);
set(gca,'ticklength',[0.05 0.05]);
axis image
title('Sagital Plane')
xlabel('depth [mm]')
ylabel('height [mm]')

%% 2.2
clear;clc
path = 'C:\Users\jingm\Desktop\Fastmarching_assignment\Fastmarching_assignment\images\MR-carotid-coronal';
[im, info] = mydicomreadfolder(path);

spacing = strsplit(info.PixelSpacing, '\');
[sx, sy, sz] = size(im);
height = [info.Rows*str2double(spacing(1)), 0];
width = [0, info.Columns*str2double(spacing{2})];
depth = [0, int8(str2double(info.SpacingBetweenSlices)*(sz-1) + str2double(info.SliceThickness)*sz)];

% Transverse plane (x=width, y=depth)
figure(4)
im_transverse = max(im,[],1);
im_transverse = reshape(im_transverse,[sx,sz]);
im_transverse = transpose(im_transverse);
% im_transverse = flip(im_transverse);
imagesc(width, depth, im_transverse)
colormap gray
colorbar
set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0]);
set(gca,'ticklength',[0.05 0.05]);
axis image
title('Transverse Plane')
xlabel('width [mm]')
ylabel('depth [mm]')

% Coronal plane (x=width, y=height)
figure(5)
im_coronal = max(im,[],3);
% im_coronal = reshape(im_coronal, [sy, sz]);
% im_coronal = transpose(im_coronal);
im_coronal = flip(im_coronal);
imagesc(width, height, im_coronal)
colormap gray
colorbar
set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0]);
set(gca,'ticklength',[0.05 0.05]);
axis image
title('Coronal Plane')
xlabel('width [mm]')
ylabel('height [mm]')

% Sagital plane (x=depth, y=height)
figure(6)
im_sagital = max(im,[],2);
im_sagital = reshape(im_sagital, [sx, sz]);
% im_sagital = transpose(im_sagital);
im_sagital = flip(im_sagital,1);
imagesc(depth, height, im_sagital)
colormap gray
colorbar
set(gca,'xcolor',[0 0 0],'ycolor',[0 0 0]);
set(gca,'ticklength',[0.05 0.05]);
axis image
title('Sagital Plane')
xlabel('depth [mm]')
ylabel('heitght [mm]')


