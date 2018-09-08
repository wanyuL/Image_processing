function trig_compare(laserlength,laserpower,ROI,plotwhich)



all=dir('*.csv');
r=csvread(all(4).name,1,1)-csvread(all(2).name,1,1);  
g=csvread(all(3).name,1,1)-csvread(all(1).name,1,1);


%% Extract frames for the green and red channel 



switch plotwhich
    case 1
figure;
roim=dir('*.png');
imshow(roim.name);
title([ROI]);
    case 2
figure;


% method 2: average within loop

plot(1:2,g,'-o','Color','g')
xticks([1 2])
xticklabels({'before','after'})

xlabel('before and after 875ms 405nm pulses shone on neuron')
ylabel('Average intensity of ROI in the green channel (A.U.) ')

title({'The response of tdEos to 405 nm flashes in the green channel', ['Laser pulse length=' laserlength ', laser power= ' laserpower ', ROI= ' ROI] })
    case 3
figure;
plot(1:2,r,'-o','Color','r')
xticks([1 2])
xticklabels({'before','after'})

xlabel('before and after 875ms 405nm pulses shone on neuron')
ylabel('Average intensity of ROI in the red channel (A.U.) ')

title({'The response of tdEos to 405 nm flashes in the red channel', ['Laser pulse length=' laserlength ', laser power= ' laserpower ', ROI= ' ROI] })

end