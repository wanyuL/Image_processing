% This function is still under construction


function visualtif(tifname)
info = imfinfo(tifname); 
n=numel(info);   % # of frames
test=imread(tifname,numel(info));  
for i=1:n
imagesc(imread(tifname,i)); 
colorbar;
% caxis([min(tifstack(:)) max(tifstack(:))]);
% caxis([quantile(tifstack(:),0.001) quantile(tifstack(:),0.98)]);   % get rid of outlier for visualization
caxis([min(tifstack(:)) quantile(tifstack(:),0.999)]);
title([inputname(1) sprintf(' frame: %d ',i)]);   % inputname get inputname as a string
pause(0.2);
end