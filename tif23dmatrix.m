% This function convert a series of .tif to 3d matrix in matlab
% usually the .tif is time-lapse movie after maximum projection in z
% direction. 
function out3dmatrix=tif23dmatrix(tifname)
info = imfinfo(tifname); 
n=numel(info);   % # of frames
test=imread(tifname,numel(info));  
out3dmatrix=zeros(size(test,1),size(test,2),n); % initialize images
for i=1:n
out3dmatrix(:,:,i)=double(imread(tifname,i));  % read in original image stacks
end