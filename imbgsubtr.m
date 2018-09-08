function imout=imbgsubtr(im,filter_box_size)

im=im2double(im);
ws= filter_box_size;
% background = imopen(im,strel('disk',10)); % remove background
% im=imsubtract(im,background);
% Ib = medfilt2(im,[ws,ws])-medfilt2(im,[3,3]);
Ib=imfilter(im,fspecial('average',ws),'replicate')-imfilter(im,fspecial('average',3),'replicate'); % difference filter
% replicate: Input array values outside the bounds of the array are assumed to equal the nearest array border value.
% for background subtraction
% Ib=imgaussfilt(im,ws)-imgaussfilt(im,ws/50);
im = imsubtract(im,Ib);
imout=imgaussfilt(im,1); % smooth rough edges, filter the image with a Gaussian filter with standard deviation of 1.
% imshow(imout,[]);

%% other code from Maijia & Sabya
% th=mean(im(:))+params.scale*std(im(:));
% BW=imbinarize(im,th);
% BW=bwareaopen(BW,300);%% remove small unwanted areas, increase threshold gives less branch points
% %BW=bwmorph(BW,'clean');%gives more branch points
% BW = sgolayfilt(double(BW),3,21, [], 1); %% Savitzky-Golay filtering to smooth the edges, for column vectors
% BW=sgolayfilt(double(BW), 3,21, [], 2); %% Savitzky-Golay filtering to smooth the edges, for row vectors
% % change framlen to smooth boundaries more and less branch points obtained,
% % set polynormial order to 3 first
% imgaussfilt(A,'FilterSize',3);  Gaussian blur 
end 
