function checkoffset(m,ch)
% basic statistics
stdm=std2(m)
meanm=mean2(m)
maxm=max(max(m))
minm=min(min(m))


figure;
s=surf(m);
s.EdgeColor = 'none';
colorbar
title(['Intensity plot of camera offset in the  ' ch ' channel (dual camera)'])
zlabel('Intensity (A.U.)')
figure;
imagesc(m,[min(min(m)),max(max(m))]);colorbar; 
title(['Intensity plot of camera offset in the  ' ch ' channel (dual camera)'])


%imhist(uint8(m)) % imhist only accept uint8 data
figure;
histogram(m);
title(['Histogram of camera offset in the ' ch ' channel (Dual camera)'])
xlabel('Intensity (A.U.)')
ylabel('Counts')


end
