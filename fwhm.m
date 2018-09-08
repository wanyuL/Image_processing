function fwhm(r2,pixelsize)
% Input: r2=csvread('red_line2.csv',1,0); pixelsize: 0.16um/pixel
% data from ImageJ; ROI data save
% from intenisty profile along Z axis 
figure;
r2(:,1)=r2(:,1)*pixelsize;
[f,gof] = fit(r2(:,1),r2(:,2),'gauss1');  
plot(f,r2(:,1),r2(:,2));
sigmar=f.c1/sqrt(2);
fwhm_r=2*sqrt(2*log(2)*sigmar); % Full Width at Half Maximum
suptitle({' ',[sprintf('FWHM=%0.5g', fwhm_r) sprintf(', Sigma=%0.5g', sigmar)],sprintf('R^2=%0.5g', gof.rsquare)});
xlabel('Length (\mum)'); ylabel('Intensity (A.U.)');
end