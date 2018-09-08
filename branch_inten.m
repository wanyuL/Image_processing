function out=branch_inten
%% 
bggname=dir('*bg_green.csv');
bgg=mean(csvread(bggname.name,1,1));
bgrname=dir('*bg_red.csv');
bgr=mean(csvread(bgrname.name,1,1));

%% 
mgname=dir('*m_green.csv');
mgt=csvread(mgname.name,1,1)-bgg;   % import data 
mg=mean(mgt(mgt>=max(mgt)/2));

mrname=dir('*m_red.csv');
mrt=csvread(mrname.name,1,1)-bgr;   % import data 
mr=mean(mrt(mrt>=max(mrt)/2));

%%
d1gname=dir('*d1_green.csv');
d1gt=csvread(d1gname.name,1,1)-bgg;   % import data 
d1g=mean(d1gt(d1gt>=max(d1gt)/2));   

d1rname=dir('*d1_red.csv');
d1rt=csvread(d1rname.name,1,1)-bgr;   % import data 
d1r=mean(d1rt(d1rt>=max(d1rt)/2));

%%
d2gname=dir('*d2_green.csv');
d2gt=csvread(d2gname.name,1,1)-bgg;   % import data 
d2g=mean(d2gt(d2gt>=max(d2gt)/2));

d2rname=dir('*d2_red.csv');
d2rt=csvread(d2rname.name,1,1)-bgr;   % import data 
d2r=mean(d2rt(d2rt>=max(d2rt)/2));

out=[mg  d1g  d2g mr  d1r  d2r mg/mr d1g/d1r d2g/d2r];


%% plot intensity profile 

subplot(1,2,1)
plot(mgt,'-o','Color','g','LineWidth',2);
hold on
plot(d1gt,'-*','Color','g','LineWidth',1);
plot(d2gt,'-*','Color','g','LineWidth',1);
legend(sprintf('mg; Int=%d',round(mg)),sprintf('d1g; Int=%d',round(d1g)),sprintf('d2g; Int=%d',round(d2g)),'Location','NorthEast')
xlabel('distance (pixel)')
ylabel('Intensity along the ROI line with background subtracted (A.U.)')
subplot(1,2,2)
plot(mrt,'-o','Color','r','LineWidth',2);
hold on
plot(d1rt,'-*','Color','r','LineWidth',1);
plot(d2rt,'-*','Color','r','LineWidth',1);
legend(sprintf('mr; Int=%d',round(mr)),sprintf('d1r; Int=%d',round(d1r)),sprintf('d2r; Int=%d',round(d2r)),'Location','NorthEast')
xlabel('distance (pixel)')
ylabel('Intensity along the ROI line with background subtracted (A.U.)')
end
