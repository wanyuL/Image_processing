% tifstack: 3D matrix (first need to convert tif stack to 3D matrix)
function visual3dmatrix(tifstack,originalfilename)
num = strfind(originalfilename,'.tif');
prefix = originalfilename(1:num-1);
colormap jet;
for i=1:size(tifstack,3)   % 3rd dimension
imagesc(tifstack(:,:,i)); 
colorbar;
% caxis([min(tifstack(:)) max(tifstack(:))]);
% caxis([quantile(tifstack(:),0.001) quantile(tifstack(:),0.98)]);   % get rid of outlier for visualization
caxis([min(tifstack(:)) quantile(tifstack(:),0.999)]);
title([inputname(1) sprintf(' frame: %d ',i)]);   % inputname get inputname as a string
% pause(0.2);
M=frame2im(getframe(gcf));


    if i==1
        imwrite(M,strcat(prefix,'_',inputname(1),'.tif') , 'WriteMode' , 'overwrite');
    else
        imwrite(M,strcat(prefix,'_',inputname(1),'.tif') , 'WriteMode' , 'append');
    end
    
end