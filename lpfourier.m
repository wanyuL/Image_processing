% This function first transfer the image to its Fourier domain, then pass a
% circular low pass filter to the fourier domain to remove the high
% frequency noise. At last, transfer the image back to its spatial domain.
% Input:m: image; i: defines the diameter of the low-pass filter;
% plot_op=1, plot all figure.
% Output: H: low-pass filter; flp: low-pass filtered image; fum: m-flp;

function [H,flp,fum]=lpfourier(m,i,plot_op)

f=double(m);
F=fftshift(fft2(double(f)));
[X1,Y1]=size(f);
X=X1;
Y=Y1;
x0=0:X-1;
y0=0:Y-1;
[y,x]=meshgrid(y0,x0); % meshgrid(col,row)
 D=sqrt((x-X/2).^2+(y-Y/2).^2);

H=ones(X,Y); %transfer function
 D0=i; % transfer function
 H(D>D0)=0;

flp=(ifft2(ifftshift(H.*F)));


%% upsharp mask
fum=double(f)-flp; % unsharp mask
%% figure
switch plot_op
    case 1
subplot(1,3,1)
imshow(H); title(sprintf('Lowpass filter, D0=%d' ,D0));
subplot(1,3,2)
imshow(real(flp),[]);title('Lowpass filtered image' )
subplot(1,3,3)
imshow(real(fum),[]);title('Unsharp mask');         
end
end
