% this function calculate the mean intensity of self-drawn background area
% modify 7/23/18 for stack images; 
% Input: the name of the image like 'mhc.tif'
% Output: bg vector; assume for each stack, the background is constant

function bg=getbg(im)
close all
            imshow(imread(im,numel(imfinfo(im))),[]);
            title('Please use your mouse to draw a circle as background')
            e = imellipse;  % impoly can also be considered
            BW_nu = createMask(e);
            for i=1:numel(imfinfo(im))
            bg_m=double(imread(im,i)).*BW_nu;
            bg(i)=mean2(nonzeros(bg_m)); 
            end
end

