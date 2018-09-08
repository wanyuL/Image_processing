% 8/4/18; used for compare line scan in green and red channels used for \\wcsfs00.its.yale.internal\axonemes-702103-mbb\Wanyu\51314_Tm\4_13_18
function compareline(roi,k)       % roi=line1
roiname=dir(['*' roi '*']);

% {ans =
% 
%     'after_405_488_100%_1st_frames_ROI_line1_dim_green.csv'
% 
% 
% ans =
% 
%     'after_405_488_100%_1st_frames_ROI_line1_dim_red.csv'
% 
% 
% ans =
% 
%     'before_405_488_100%_ROI_line1_dim_green.csv'
% 
% 
% ans =
% 
%     'before_405_488_100%_ROI_line1_dim_red.csv'
  
br=csvread(roiname(4).name,1,0);
bg=csvread(roiname(3).name,1,0);
ar=csvread(roiname(2).name,1,0);
ag=csvread(roiname(1).name,1,0);
figure;
plot(br(:,1),br(:,2),'--','LineWidth',2,'Color','red');
hold on 
plot(bg(:,1),bg(:,2),'--','LineWidth',2,'Color','green');
plot(ar(:,1),ar(:,2),'-','LineWidth',2,'Color','red');
plot(ag(:,1),ag(:,2),'-','LineWidth',2,'Color','green')
hold off
xlabel('Distance (\mum)');
ylabel('Intensity (A.U.)');
xlim([0 max(ar(:,1))]);
l=legend('before; red','before; green','after 405; red','after 405; green','Location','best');
legend('boxoff');
title(['Intensity profile of ' roi],'Color','w');
%% Save for blackground power_point
if k==1
set(gca,'Color','k')
set(gca,'YColor','w','FontSize',12)
set(gca,'XColor','w','FontSize',12)
set(l,'TextColor',[1 1 1],'EdgeColor',[0 0 0],'Location','best')  
end
mkdir linescan_fig_matlab
 % cd linescan_fig_matlab
saveas(gcf,[roi '.fig'])
% saveas(gcf,[roi '.png'])
%%
if k==1
    clf
else
    figure
end
plot(br(:,1),(bg(:,2)-100)./(br(:,2)-100),'--','LineWidth',2,'Color','white');
hold on 
plot(ar(:,1),(ag(:,2)-100)./(ar(:,2)-100),'-','LineWidth',2,'Color','white');
hold off
xlabel('Distance (\mum)');
ylabel('(Green-100)/(Red-100)');
l=legend('before','after 405','Location','best');
legend('boxoff');
title(['Intensity profile of ' roi],'Color','w');
set(gca,'Color','k')
set(l,'TextColor',[1 1 1]);
xlim([0 max(ar(:,1))]);
if k==1
set(gca,'YColor','w','FontSize',12)
set(gca,'XColor','w','FontSize',12)
set(l,'TextColor',[1 1 1],'EdgeColor',[0 0 0],'Location','best');
saveas(gcf,[roi '_norm.fig'])
end

end