function linesum(roiname)
roi=imread([roiname '.tif']);
t=0:size(roi,1)-1;
plot(t*10,sum(roi,2),'Color','white');
set(gca,'Color','k')
set(gca,'YColor','w','FontSize',12)
set(gca,'XColor','w','FontSize',12)
xlabel('Time (s)');
ylabel('Sum of Intensity along ROI (A.U.)');
xlim([10*min(t) 10*max(t)])
saveas(gcf,[roiname 'sum_along_line.fig'])
end