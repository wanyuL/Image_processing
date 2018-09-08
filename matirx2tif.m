% This function converts 3D matrix to tif file ; Wanyu 7/5/18
% Input: originalfilename is the '.... . tif'; matlabfig is matlab figure after image
% processing 
function matirx2tif(originalfilename,matlabfig) 
num = strfind(originalfilename,'.');
prefix = originalfilename(1:num-1);
for i = 1 : size(matlabfig, 3)
    imwrite(matlabfig(:,:,i) , strcat(prefix,'_',inputname(2),'.tif') , 'WriteMode' , 'append') ; % write a tif file, appending each image as a new page
end
end