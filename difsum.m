% This function calculates the intensity derivative of each pixel than sum
% them up as a metric to characterize dynacmis. 
% Input: '.tif figure time-lapse video' Wanyu 7/5/18
% Output
function difsum(tifstack)
num = strfind(tifstack,'.tif');
prefix = tifstack(1:num-1);
difsum=zeros(numel(imfinfo(tifstack))-1);   % initialize matrix
for i=1:numel(imfinfo(tifstack))-1
    temdif=imread(tifstack,i+1)-imread(tifstack,i);
    difsum(i)=squeeze(sum(abs(temdif(:))));   % squeeze remove singleton dimension 
end 
plot(1:numel(imfinfo(tifstack))-1,difsum);
xlabel('Number of Frame'); ylabel('Sum of Intensity Derivative (A.U.)');
title({tifstack(1:end),''},'Interpreter', 'none');
saveas(gcf,strcat(prefix,'_difsum.png'),'png');
end