function BW=imbimask(im,scale)
th=mean(im(:))+scale*std(im(:));
BW=imbinarize(im,th);
% imshow(BW,[]);
end