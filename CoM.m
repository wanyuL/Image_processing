% Input im is image info from dir('*.tif'); scale is the multipler of sigma
% for binarization
% Out is .tif with center of mass of thresholded patches labeld
% 7/23/18 Wanyu
function CoM(im,conne,scale,timeinterval,pixelsize)
close all
bg=getbg(im.name);
figure;
firstframe=double(imread(im.name,1));
set(0,'defaultfigureposition',[100 200 1600 800]);  % change image size, left bottom width height
imwrite(frame2im(getframe(gcf)),sprintf('WeightedCentroid_threshold_scale_%g.tif',scale),'WriteMode','overwrite');
colormap parula; % colormap jet
for i=1:numel(imfinfo(im.name))
    disp(i)
  %  [J,~] = wiener2(imread(im.name,i),[4 4]);
J=imbgsubtr(imread(im.name,i),20);   % can change the parameter here
BW(:,:,i)=bwareaopen(imbimask(double(J),scale),conne);
% BW(:,:,i)=imbimask(double(J),2);
%     if i==1
%     imwrite(BW(:,:,i),'BW.tif','WriteMode','overwrite');
%     else
%     imwrite(BW(:,:,i),'BW.tif','WriteMode','append');   
%     end    
s = regionprops(BW(:,:,i),mat2gray(BW(:,:,i).*double(imread(im.name,i))),'WeightedCentroid');  % Calculate centroids for connected components in the image
centroids = cat(1, s.WeightedCentroid);
suptitle(sprintf('Frame=%g; Time=%g s',i,i*timeinterval));
subplot(1,2,1)
imagesc(BW(:,:,i).*(double(imread(im.name,i))-bg(i)));
colorbar
% lim = caxis

caxis([0 quantile(firstframe(:),0.95)])  % can change the parameter here
xlabel(sprintf('pixel, size=%g um/pixel',pixelsize));
ylabel(sprintf('pixel, size=%g um/pixel',pixelsize));
title({sprintf('Centers of mass of hotspots (*) (> %g connected pixels)',conne),sprintf( 'in thresholded (intensity > mean + %g * std) image',scale)});
hold on
plot(centroids(:,1),centroids(:,2), 'r*');
hold off
subplot(1,2,2)
imagesc(double(imread(im.name,i))-bg(i));
xlabel(sprintf('pixel, size=%g um/pixel',pixelsize));
ylabel(sprintf('pixel, size=%g um/pixel',pixelsize));
title('Original image (constant background subtraction)')
colorbar 
% lim = caxis;
caxis([0 quantile(firstframe(:),0.95)]);
% pause(1);
imwrite(frame2im(getframe(gcf)),sprintf('WeightedCentroid_threshold_scale_%g.tif',scale),'WriteMode','append'); 
end
% fullfile(folder,'binaryeye.jpg')
% quiver centroid; subtract background; 