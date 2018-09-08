function trig_profile(laserlength,laserpower,ROI,plotwhich)



all=dir('*.csv');

% before 8_18_18
% r=csvread(all(2).name,1,1);  % only remain intensity info, delete first row ( variable names), first column (frame numbers)
% g=csvread(all(1).name,1,1);

%modify 8_18_18; include background subtraction
r=csvread(all(4).name,1,1)-csvread(all(2).name,1,1);  
g=csvread(all(3).name,1,1)-csvread(all(1).name,1,1);
% r=csvread(all(4).name,1,1);  
% g=csvread(all(3).name,1,1);
f=1:numel(r); % total frame number 
nl=1:numel(r)/5;   % number of loops

%% Extract frames for the green and red channel 


% method 1: only one frame 
% rfn=sort([5*(nl'-1)+2;5*(nl'-1)+4]); % frame numbers of the green signal, sort it in ascending order
% gfn=sort([5*(nl'-1)+1;5*(nl'-1)+3]); % frame numbers of the red signal
% rr=r(rfn); % red signal when 561 nm is on
% gg=g(gfn);  % green signal when 561 nm is on
% plot(rr,'o')
% plot(gg,'o')
% 
rfn1=5*(nl'-1)+2; % frame numbers of the green signal, sort it in ascending order
rfn2=5*(nl'-1)+4; % frame numbers of the red signal, sort it in descending order

gfn1=5*(nl'-1)+1; % frame numbers of the red signal
gfn2=5*(nl'-1)+3;
rr=(r(rfn1)+r(rfn2))/2; % red signal when 561 nm is on
gg=(g(gfn1)+g(gfn2))/2;  % green signal when 561 nm is on
% rr=r(rfn1); % red signal when 561 nm is on
% gg=g(gfn2);  % green signal when 561 nm is on
switch plotwhich
    case 1
figure;
roim=dir('*.png');
imshow(roim.name);
title([ROI]);
    case 2
figure;


% method 2: average within loop

plot(nl-1,gg,'-o','Color','g')

xlabel('Number of 875ms 405nm pulses shone on neuron')
ylabel('Average intensity of ROI in the green channel (A.U.) ')

title({'The response of tdEos to 405 nm flashes in the green channel', ['Laser pulse length=' laserlength ', laser power= ' laserpower ', ROI= ' ROI] })
    case 3
figure;
plot(nl-1,rr,'-o','Color','r');

xlabel('Number of 875ms 405nm pulses shone on neuron')
ylabel('Average intensity of ROI in the red channel (A.U.) ')
title({'The response of tdEos to 405 nm flashes in the red channel', ['Laser pulse length=' laserlength ', laser power= ' laserpower ', ROI= ' ROI] })

end
% fit an exponential