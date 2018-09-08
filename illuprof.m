function illuprof(im,offset)
figure;
m=double(imread(im))-offset; % camera offset subtract
imshow(m,[]); title('Grayscale original image');


figure;
imagesc(m);cb=colorbar; lim=caxis;
title({'Illumination profile on 40x objective, entire FOV shown','Camera offset (Intensity 100 A.U.) subtracted'});
title(cb,'Intensity (A.U.)');xlabel('pixel number (1 pixel=0.16\mum)');ylabel('pixel number (1pixel=0.16\mum)');
figure;
mesh(m-100);zlabel('Intensity (A.U.)');title('3D plot of the illumination profile (Camera offset subtracted)');
xlabel('pixel number (1 pixel=0.16\mum)');ylabel('pixel number (1 pixel=0.16\mum)')




BW=zeros(size(m));BW1=BW;BW2=BW;
BW1(:,[floor(median(1:size(m,2))) ceil(median(1:size(m,2)))])=1;

figure;
BW2([floor(median(1:size(m,2))),ceil(median(1:size(m,2)))],:)=1;
subplot(1,2,1);
imshow(BW2,[]); title('Horizontal line ROI (white, width=2 pixels)')
subplot(1,2,2);

plot(1:size(m,1),sum(BW2.*m,1)/2);  
xlabel('Distance along the ROI line (pixel)');
ylabel('Intensity (A.U.)');xlim([1 size(m,1)]);
title('Intensity profile along the ROI line');

figure;

subplot(1,2,1);
imshow(BW1,[]); title('Vertical line ROI (white, width=2 pixels)')
subplot(1,2,2);
plot(1:size(m,2),sum(BW1.*m,2)/2);  
xlabel('Distance along the ROI line (pixel)');
ylabel('Intensity (A.U.)');xlim([1 size(m,2)]);
title('Intensity profile along the ROI line');



figure;
BW3=BW;
BW3(logical(eye(size(m))))=1;
subplot(1,2,1);
imshow(BW3,[]); title('Line ROI (white, width=1 pixels)')
subplot(1,2,2);
plot(1:size(m,1),sum(BW3.*m,1));  
xlabel('Distance along the ROI line (pixel)');
ylabel('Intensity (A.U.)');xlim([1 size(m,1)]);
title('Intensity profile along the ROI line');
end