% This function plot the difference between adjacent images
% Input: im is generated from dir('mhc.tif')
% Output: tif file series saved in file directory

function difim(im,timeinterval,pixelsize)
close all
colormap default
bg=getbg(im.name);
figure
set(0,'defaultfigureposition',[100 200 1600 800]);
firstframe=double(imread(im.name,1));
firstf=double(imread(im.name,2))-double(imread(im.name,1)); %need to convert to double than do subtraction
for i=1:numel(imfinfo(im.name))-1
%% Dif image
subplot(1,2,1);
imagesc( imgaussfilt(double(imread(im.name,i+1))-double(imread(im.name,i)),'FilterSize',3));   % note a Gaussian blur, radius
colorbar;
caxis([-quantile(firstf(:),0.95) quantile(firstf(:),0.95)]);
xlabel(sprintf('pixel, size=%g um/pixel',pixelsize));
ylabel(sprintf('pixel, size=%g um/pixel',pixelsize));
title({'First-order derivative of image series','Difim(t)=im(t+1)-im(t)'});
%% Original image
subplot(1,2,2)
imagesc(double(imread(im.name,i))-bg(i));
title('Original image (constant background subtraction)')
xlabel(sprintf('pixel, size=%g um/pixel',pixelsize));
ylabel(sprintf('pixel, size=%g um/pixel',pixelsize));
colorbar 
% lim = caxis;
caxis([0 quantile(firstframe(:),0.95)]);
%%
suptitle(sprintf('Frame=%g; Time=%g s',i,i*timeinterval));
%%
    if i==1
    imwrite(frame2im(getframe(gcf)),'Dif.tif','WriteMode','overwrite');
    else
    imwrite(frame2im(getframe(gcf)),'Dif.tif','WriteMode','append'); 
    end
end
end