%Downsample 1:4:end before var projection. 161122mjs

% This is Step 3 for processing calcium imaging data.
%
% After Step 2: downsample/stitching .tiff -->
% Graphical user interface for selecting ROI
% --> in preparation for Step 4: extract dF/F from ROIs
%
% Kwan Lab, 6/15/2016

function varargout = cellROI(varargin)
%To edit/open this program, use "File --> New --> GUI"

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cellselect_OpeningFcn, ...
                   'gui_OutputFcn',  @cellselect_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
 
% --- Executes just before cellselect is made visible.
function cellselect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cellselect (see VARARGIN)

% Choose default command line output for cellselect
clear global;
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = cellselect_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in loadimagestack.
function pushbutton_loadimagestack_Callback(hObject, eventdata, handles)
% hObject    handle to loadimagestack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global;
global pic
global rightpanel_pic leftpanel_pic %what is being displayed on left/right
global mean_pic max_pic var_pic red_pic %different types of projections
global bw 
global savedroi_pic
global nX nY nZ
global cellf
global pathname

cellf=[];           %time-lapse fluorescence for the current roi

[filename, pathname] = uigetfile('*.tif', 'Pick a TIF file');
cd(pathname);
disp(['loading ' pathname filename '...'])

%load the tiff
pic=loadtiffseq(pathname,filename);
[nY nX nZ]=size(pic);
set(handles.edit8,'String',filename(1:end-4));
set(handles.text1,'String',[num2str(nZ) ' frames of ' num2str(nX) 'x' num2str(nY) ' loaded.']);
savedroi_pic=false(nY,nX);  %mask of all the saved roi
bw=false(nY,nX);            %mask of the currently selected roi

%set up the sub-folder for saving files
mkdir(['ROI_' filename]);
cd(['ROI_' filename]);
if (isunix)   %it is a Mac/Unix
    pathname=[pwd '/'];
else    %it is a PC
    pathname=[pwd '\'];
end

disp('Getting mean projection...')
temp = double(pic(:,:,1:4:end));
mean_pic=mean(double(pic),3);
disp('Getting variance projection...')
%var_pic=var(double(pic),[],3); 
var_pic=var(temp,[],3); 
disp('Getting maximum projection...')
%max_pic=max(double(pic),[],3);  
max_pic=max(temp,[],3); 
red_pic=ones(size(mean_pic));

%plot the right panel
if get(handles.radiobutton4,'Value')==1
    rightpanel_pic=mean_pic;
elseif get(handles.radiobutton5,'Value')==1
    rightpanel_pic=max_pic;
else
    rightpanel_pic=var_pic;
end
zprojWhiteLev=str2num(get(handles.edit11,'String'));
temp_pic=round(zprojWhiteLev*253)*rightpanel_pic/nanmax(rightpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
temp_pic(temp_pic>252)=252;
cla(handles.axes1);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes1,graymap);
image(temp_pic,'Parent',handles.axes1);
axis(handles.axes1, 'off', 'equal');

%plot the left panel
if get(handles.radiobutton10,'Value')==1
    leftpanel_pic=mean_pic;
elseif get(handles.radiobutton11,'Value')==1
    leftpanel_pic=max_pic;
elseif get(handles.radiobutton12,'Value')==1
    leftpanel_pic=var_pic;
end
leftWhiteLev=str2num(get(handles.edit10,'String'));
temp_pic=uint16(round(leftWhiteLev*253)*leftpanel_pic/nanmax(leftpanel_pic(:)));
temp_pic(temp_pic>252)=252;
cla(handles.axes4);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes4,graymap);
image(temp_pic,'Parent',handles.axes4);
axis(handles.axes4, 'off', 'equal');

function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global leftpanel_pic rightpanel_pic %what is being displayed on left/right
global mean_pic max_pic var_pic red_pic %different types of projections
global savedroi_pic

rightpanel_pic=mean_pic;
zprojWhiteLev=str2num(get(handles.edit11,'String'));
roi_pic=round(zprojWhiteLev*253)*rightpanel_pic/nanmax(rightpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roi_pic(roi_pic>252)=252;
roi_pic(savedroi_pic)=254;
cla(handles.axes1);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes1,graymap);
image(roi_pic,'Parent',handles.axes1);
axis(handles.axes1, 'off', 'equal');

function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global leftpanel_pic rightpanel_pic %what is being displayed on left/right
global mean_pic max_pic var_pic red_pic %different types of projections
global savedroi_pic

rightpanel_pic=max_pic;
zprojWhiteLev=str2num(get(handles.edit11,'String'));
roi_pic=round(zprojWhiteLev*253)*rightpanel_pic/nanmax(rightpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roi_pic(roi_pic>252)=252;
roi_pic(savedroi_pic)=254;
cla(handles.axes1);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes1,graymap);
image(roi_pic,'Parent',handles.axes1);
axis(handles.axes1, 'off', 'equal');

function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global leftpanel_pic rightpanel_pic %what is being displayed on left/right
global mean_pic max_pic var_pic red_pic %different types of projections
global savedroi_pic

rightpanel_pic=var_pic;
zprojWhiteLev=str2num(get(handles.edit11,'String'));
roi_pic=round(zprojWhiteLev*253)*rightpanel_pic/nanmax(rightpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roi_pic(roi_pic>252)=252;
roi_pic(savedroi_pic)=254;
cla(handles.axes1);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes1,graymap);
image(roi_pic,'Parent',handles.axes1);
axis(handles.axes1, 'off', 'equal');

function radiobutton10_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global leftpanel_pic rightpanel_pic %what is being displayed on left/right
global mean_pic max_pic var_pic red_pic %different types of projections
global savedroi_pic

leftpanel_pic=mean_pic;
leftWhiteLev=str2num(get(handles.edit10,'String'));
roi_pic=round(leftWhiteLev*253)*leftpanel_pic/nanmax(leftpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roi_pic(roi_pic>252)=252;
roi_pic(savedroi_pic)=254;
cla(handles.axes4);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes4,graymap);
image(roi_pic,'Parent',handles.axes4);
axis(handles.axes4, 'off', 'equal');

function radiobutton11_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global leftpanel_pic rightpanel_pic %what is being displayed on left/right
global mean_pic max_pic var_pic red_pic %different types of projections
global savedroi_pic

leftpanel_pic=max_pic;
leftWhiteLev=str2num(get(handles.edit10,'String'));
roi_pic=round(leftWhiteLev*253)*leftpanel_pic/nanmax(leftpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roi_pic(roi_pic>252)=252;
roi_pic(savedroi_pic)=254;
cla(handles.axes4);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes4,graymap);
image(roi_pic,'Parent',handles.axes4);
axis(handles.axes4, 'off', 'equal');

function radiobutton12_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global leftpanel_pic rightpanel_pic %what is being displayed on left/right
global mean_pic max_pic var_pic red_pic %different types of projections
global savedroi_pic

leftpanel_pic=var_pic;
leftWhiteLev=str2num(get(handles.edit10,'String'));
roi_pic=round(leftWhiteLev*253)*leftpanel_pic/nanmax(leftpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roi_pic(roi_pic>252)=252;
roi_pic(savedroi_pic)=254;
cla(handles.axes4);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';

image(roi_pic,'Parent',handles.axes4);
colormap(handles.axes4,graymap);
axis(handles.axes4, 'off', 'equal');

function radiobutton13_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global leftpanel_pic rightpanel_pic %what is being displayed on left/right
global mean_pic max_pic var_pic red_pic shift_pic %different types of projections
global savedroi_pic

shift = [-str2num(get(handles.edit16,'String')) str2num(get(handles.edit15,'String'))];
shift_pic = circshift(red_pic, shift);

leftpanel_pic=shift_pic;
leftWhiteLev=str2num(get(handles.edit10,'String'));
roi_pic=round(leftWhiteLev*253)*leftpanel_pic/nanmax(leftpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roi_pic(roi_pic>252)=252;
roi_pic(savedroi_pic)=254;
cla(handles.axes4);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';

image(roi_pic,'Parent',handles.axes4);
colormap(handles.axes4,graymap);
axis(handles.axes4, 'off', 'equal');

% --- Executes on button press in pushbutton_selectcircle.
function pushbutton_selectcircle_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectcircle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global pic rightpanel_pic leftpanel_pic shift_pic
global  bw cellmask savedroi_pic
global nX nY nZ
global cellf
global pathname

%draw the circle
r=str2num(get(handles.edit7,'String'));
[cellx celly]=ginput(1);

if (cellx>0) && (cellx<nX) && (celly>0) && (celly<nY)
    axes(handles.axes1);
    hold(handles.axes1,'on');
    viscircles([cellx celly],r,'EdgeColor','r');
    axes(handles.axes4);
    hold(handles.axes4,'on');
    viscircles([cellx celly],r,'EdgeColor','r');
    
    %determine the circle mask
    [rr cc]=meshgrid(1:nX);
    bw=sqrt((rr-cellx).^2+(cc-celly).^2)<=r;
    
    for i=1:1:nZ
        cellf(i)=sum(sum(pic(:,:,i).*uint16(bw)));
    end
    cellf=cellf/sum(sum(bw));   %per-pixel fluorescence
    
    %plot the dF/F of the selected ROI on the bottom panel
    set(handles.text3,'String',[num2str(sum(sum(bw))) ' pixels. F(t=0)=' num2str(cellf(1)) '.']);
    cla(handles.axes3);
    plot(handles.axes3,[1:1:nZ],(cellf-median(cellf))/median(cellf),'k');
    axis(handles.axes3, 'tight');
    
    %if selecting circle, then turn off the auto-select poly function
    set(handles.checkbox7, 'Value',0);
end

% --- Executes on button press in pushbutton_selectPoly.
function pushbutton_selectPoly_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectPoly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pic rightpanel_pic leftpanel_pic shift_pic
global  bw cellmask savedroi_pic
global nX nY nZ
global cellf
global pathname

%are we selecting from left or right panel
if get(handles.radiobutton13,'Value')==1
    axes(handles.axes4);
else
    axes(handles.axes1);
end

%let user draw the polygon to get ROI and associated fluorescence
[bw,xi,yi] = roipoly;

axes(handles.axes1);
hold(handles.axes1,'on');
plot(handles.axes1,xi,yi,'r','LineWidth',1); 
axes(handles.axes4);
hold(handles.axes4,'on');
plot(handles.axes4,xi,yi,'r','LineWidth',1); 

for i=1:1:nZ    
    cellf(i)=sum(sum(pic(:,:,i).*uint16(bw)));
end
cellf=cellf/sum(sum(bw));   %per-pixel fluorescence

%plot the dF/F of the selected ROI on the bottom panel
set(handles.text3,'String',[num2str(sum(sum(bw))) ' pixels. F(t=0)=' num2str(cellf(1)) '.']);
cla(handles.axes3);
plot(handles.axes3,[1:1:nZ],(cellf-median(cellf))/median(cellf),'k');
axis(handles.axes3, 'tight');

%if selecting polygon, then turn off the auto-select circle function
set(handles.checkbox6, 'Value',0);

% --- Executes on button press in savetraces.
function pushbutton_savetraces_Callback(hObject, eventdata, handles)
% hObject    handle to savetraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rightpanel_pic leftpanel_pic red_pic shift_pic
global bw cellmask savedroi_pic
global xcorrpic
global nX nY nZ
global cellf
global pathname

filenameheader=get(handles.edit8,'String');
cellnumber=get(handles.edit2,'String');
cellmask(:,:,str2num(cellnumber))=bw;

%if red image is opened, then save a duplicate file in another directory
if get(handles.radiobutton13, 'Value')==1
    isRedCell=1;
else
    isRedCell=0;
end

temp = sprintf('%03d',str2num(cellnumber));
save(strcat(pathname,filenameheader,'_cell',temp,'.mat'),'cellf','bw','isRedCell');

%if red image is opened, then save a duplicate file in another directory
if get(handles.radiobutton13, 'Value')==1
    if (isunix)   %it is a Mac/Unix
        pathname_specCell=strcat(pathname,'SpecialCells/');
    else    %it is a PC
        pathname_specCell=strcat(pathname,'SpecialCells\');
    end
    if exist(pathname_specCell,'dir')~=7
        mkdir(pathname_specCell);
    end
    save(strcat(pathname_specCell,filenameheader,'_cell',temp,'.mat'),'cellf','bw','isRedCell');
end

set(handles.edit2,'String',num2str(str2num(cellnumber)+1));

zprojWhiteLev=str2num(get(handles.edit11,'String'));
savedroi_pic=savedroi_pic | bw;
roi_pic=round(zprojWhiteLev*253)*rightpanel_pic/nanmax(rightpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roi_pic(roi_pic>252)=252;
roi_pic(savedroi_pic)=254;
cla(handles.axes1);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';

image(roi_pic,'Parent',handles.axes1);
colormap(handles.axes1,graymap);
axis(handles.axes1, 'off', 'equal');

leftWhiteLev=str2num(get(handles.edit10,'String'));
roileftpanel_pic=round(leftWhiteLev*253)*leftpanel_pic/nanmax(leftpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roileftpanel_pic(roileftpanel_pic>252)=252;
roileftpanel_pic(savedroi_pic)=254;
cla(handles.axes4);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';

image(roileftpanel_pic,'Parent',handles.axes4);
colormap(handles.axes4,graymap);
axis(handles.axes4, 'off', 'equal');

if(get(handles.checkbox6, 'Value')==1)
    pushbutton_selectcircle_Callback(handles.pushbutton_selectcircle,[],handles);
elseif(get(handles.checkbox7, 'Value')==1)
    pushbutton_selectPoly_Callback(handles.pushbutton_selectPoly,[],handles);
end

% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
global rightpanel_pic leftpanel_pic red_pic shift_pic
global bw cellmask savedroi_pic
global xcorrpic
global nX nY nZ
global cellf
global pathname

% determine the key that was pressed
keyPressed = eventdata.Key;
if strcmpi(keyPressed,'s') %save
    % set focus to the button
    uicontrol(handles.savetraces);
    % call the callback
    pushbutton_savetraces_Callback(handles.savetraces,[],handles);
elseif strcmpi(keyPressed,'a') %undo
    %reset right panel
    zprojWhiteLev=str2num(get(handles.edit11,'String'));
    roi_pic=round(zprojWhiteLev*253)*rightpanel_pic/nanmax(rightpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
    roi_pic(roi_pic>252)=252;
    roi_pic(savedroi_pic)=254;
    cla(handles.axes1);
    graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
    colormap(handles.axes1,graymap);
    image(roi_pic,'Parent',handles.axes1);
    axis(handles.axes1, 'off', 'equal');
    % set focus to the button
    uicontrol(handles.pushbutton_selectcircle);
    % call the callback
    pushbutton_selectcircle_Callback(handles.pushbutton_selectcircle,[],handles);
end

% --- Executes on button press in selectXCorr.
function selectXCorr_Callback(hObject, eventdata, handles)
% hObject    handle to savetraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pic rightpanel_pic savedroi_pic shift_pic
global leftpanel_pic bw cellmask
global xcorrpic
global nX nY nZ
global cellf
global pathname

[cellx celly]=ginput(1);    

boxLen=str2num(get(handles.edit9,'String'));
if isempty(boxLen)
    x1=1;
    x2=nX;
    y1=1;
    y2=nY;
else
    x1=round(cellx)-round(boxLen/2);
    x2=round(cellx)+round(boxLen/2);
    y1=round(celly)-round(boxLen/2);
    y2=round(celly)+round(boxLen/2);
    if x1<1
        x1=1;
    end
    if x2>nX
        x2=nX;
    end
    if y1<1
        y1=1;
    end
    if y2>nY
        y2=nY;l
    end
end

xcorrpic=zeros(nY,nX);
for i=x1:x2
    for j=y1:y2
        %double floating faster than single floating...
        xcorrpic(j,i)=corr2(double(squeeze(pic(round(celly),round(cellx),:))),double(squeeze(pic(j,i,:))));
    end
end

%use percentile as estimate on threshold, can be manually recalc later
temp=xcorrpic(:);
temp=temp(temp>0);
lowXcorr=prctile(temp,90); %75th
highXcorr=1;
set(handles.edit5,'String',num2str(lowXcorr));
set(handles.edit3,'String',num2str(highXcorr));
    
cla(handles.axes6);
edges=[-1:0.02:1];
n=histc(xcorrpic(:),edges);
bar(handles.axes6,edges,n,'histc');
hold(handles.axes6);
plot(handles.axes6,[lowXcorr lowXcorr],[0 1.1*max(n)],'r','LineWidth',3);
plot(handles.axes6,[highXcorr highXcorr],[0 1.1*max(n)],'r','LineWidth',3);
axis(handles.axes6,[0 1 0 1.1*max(n(round(numel(n)*3/4):end))]);
hold(handles.axes6);

bw=(xcorrpic>lowXcorr) & (xcorrpic<=highXcorr);
for i=1:1:nZ    
    cellf(i)=sum(sum(pic(:,:,i).*uint16(bw)));
end
cellf=cellf/sum(sum(bw));   %per-pixel fluorescence

zprojWhiteLev=str2num(get(handles.edit11,'String'));
roi_pic=round(zprojWhiteLev*253)*rightpanel_pic/nanmax(rightpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roi_pic(roi_pic>252)=252;
roi_pic(savedroi_pic)=254;
roi_pic(bw)=255;    %add newly selected ROI as red
cla(handles.axes1);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes1,graymap);
image(roi_pic,'Parent',handles.axes1);
axis(handles.axes1, 'off', 'equal');

set(handles.text3,'String',[num2str(sum(sum(bw))) ' pixels. F(t=0)=' num2str(cellf(1)) '.']);

cla(handles.axes3);
plot(handles.axes3,[1:1:nZ],(cellf-median(cellf))/median(cellf),'k');
axis(handles.axes3, 'tight');

leftWhiteLev=str2num(get(handles.edit10,'String'));
roiPic_overlay=round(leftWhiteLev*253)*leftpanel_pic/nanmax(leftpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roiPic_overlay(roiPic_overlay>252)=252;
roiPic_overlay(savedroi_pic)=254;
roiPic_overlay(bw)=255;    %add newly selected ROI as red
cla(handles.axes4);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes4,graymap);
image(roiPic_overlay,'Parent',handles.axes4);
axis(handles.axes4, 'off', 'equal');

%if selecting xcorr, then turn off the auto-select circle/poly function
set(handles.checkbox6, 'Value',0);
set(handles.checkbox7, 'Value',0);

% --- Executes on button press in RecalcXcorr.
function RecalcXcorr_Callback(hObject, eventdata, handles)
% hObject    handle to RecalcXcorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global pic rightpanel_pic savedroi_pic
global leftpanel_pic bw
global xcorrpic
global nX nY nZ
global cellf
global pathname

lowXcorr=str2num(get(handles.edit5,'String'));
highXcorr=str2num(get(handles.edit3,'String'));

cla(handles.axes6);
edges=[-1:0.02:1];
n=histc(xcorrpic(:),edges);
bar(handles.axes6,edges,n,'histc'); 
hold(handles.axes6);
plot(handles.axes6,[lowXcorr lowXcorr],[0 1.1*max(n)],'r','LineWidth',3);
plot(handles.axes6,[highXcorr highXcorr],[0 1.1*max(n)],'r','LineWidth',3);
axis(handles.axes6,[0 1 0 1.1*max(n(round(numel(n)*3/4):end))]);
hold(handles.axes6);

bw=(xcorrpic>lowXcorr) & (xcorrpic<=highXcorr);
for i=1:1:nZ    
    cellf(i)=sum(sum(pic(:,:,i).*uint16(bw)));
end
cellf=cellf/sum(sum(bw));   %per-pixel fluorescence

zprojWhiteLev=str2num(get(handles.edit11,'String'));
roi_pic=round(zprojWhiteLev*253)*rightpanel_pic/nanmax(rightpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roi_pic(roi_pic>252)=252;
roi_pic(savedroi_pic)=254;
roi_pic(bw)=255;    %add newly selected ROI as red
cla(handles.axes1);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes1,graymap);
image(roi_pic,'Parent',handles.axes1);
axis(handles.axes1, 'off', 'equal');

set(handles.text3,'String',[num2str(sum(sum(bw))) ' pixels. F(t=0)=' num2str(cellf(1)) '.']);

cla(handles.axes3);
plot(handles.axes3,[1:1:nZ],(cellf-median(cellf))/median(cellf),'k');
axis(handles.axes3, 'tight');

% --- Executes on button press in redrawVarPic.
function redrawVarPic_Callback(hObject, eventdata, handles)
% hObject    handle to redrawVarPic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global red_pic shift_pic
global leftpanel_pic savedroi_pic

%if we are drawing the red pic
if get(handles.radiobutton13,'Value')==1
    shift = [-str2num(get(handles.edit16,'String')) str2num(get(handles.edit15,'String'))];
    shift_pic = circshift(red_pic, shift);
    
    leftpanel_pic = shift_pic;
end

leftWhiteLev=str2num(get(handles.edit10,'String'));
temp_pic=round(leftWhiteLev*253)*leftpanel_pic/nanmax(leftpanel_pic(:));
temp_pic(temp_pic>=252)=252;
temp_pic(savedroi_pic)=254;
cla(handles.axes4);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
image(temp_pic,'Parent',handles.axes4);
colormap(handles.axes4,graymap);
axis(handles.axes4, 'off', 'equal');

% --- Executes on button press in RedrawZProjPic.
function RedrawZProjPic_Callback(hObject, eventdata, handles)
% hObject    handle to RedrawZProjPic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rightpanel_pic savedroi_pic

zprojWhiteLev=str2num(get(handles.edit11,'String'));
temp_pic=round(zprojWhiteLev*253)*rightpanel_pic/nanmax(rightpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
temp_pic(temp_pic>252)=252;
temp_pic(savedroi_pic)=254;
cla(handles.axes1);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';

image(temp_pic,'Parent',handles.axes1);
colormap(handles.axes1,graymap);
axis(handles.axes1, 'off', 'equal');


% --- Executes on button press in loadoldROI.
function loadoldROI_Callback(hObject, eventdata, handles)
% hObject    handle to loadoldROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global pic rightpanel_pic savedroi_pic shift_pic
global leftpanel_pic bw cellmask
global xcorrpic
global nX nY nZ
global cellf
global pathname

numOldROIs=str2num(get(handles.edit12,'String'));
savedroi_pic=false(nY,nX);  %mask of all the saved roi
cellmask=false(nY,nX);

nonCell_masks = false;
for i=1:numOldROIs
    filenameheader=get(handles.edit8,'String');
    try
        load(strcat(pathname,filenameheader,'_cell',sprintf('%03d',i),'.mat'))
    catch
        if nonCell_masks == false
           temp = i; %get idx of first non-existent ROI ref (to allow load w/ non-cell masks ref'ed with arbitrarily large idx)
           warning('Some ROIs referenced do not exist.');
           nonCell_masks = true;
        end
    end
    savedroi_pic=savedroi_pic | bw;
    cellmask(:,:,i)=bw;
end    

if nonCell_masks == true
    set(handles.edit2,'String',num2str(temp)); %modified to continue after last consec ROI (in case of non-cell masks ref'ed with arbitrarily large idx)
else
    set(handles.edit2,'String',num2str(numOldROIs+1));
end

zprojWhiteLev=str2num(get(handles.edit11,'String'));
roi_pic=round(zprojWhiteLev*253)*rightpanel_pic/nanmax(rightpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roi_pic(roi_pic>252)=252;
roi_pic(savedroi_pic)=254;
cla(handles.axes1);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes1,graymap);
image(roi_pic,'Parent',handles.axes1);
axis(handles.axes1, 'off', 'equal');

leftWhiteLev=str2num(get(handles.edit10,'String'));
roileftpanel_pic=round(leftWhiteLev*253)*leftpanel_pic/nanmax(leftpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roileftpanel_pic(roileftpanel_pic>252)=252;
roileftpanel_pic(savedroi_pic)=254;
cla(handles.axes4);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes4,graymap);
image(roileftpanel_pic,'Parent',handles.axes4);
axis(handles.axes4, 'off', 'equal');

% --- Executes on button press in calcNeuropil.
function calcNeuropil_Callback(hObject, eventdata, handles)
% hObject    handle to calcNeuropil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global pic rightpanel_pic savedroi_pic
global leftpanel_pic bw cellmask
global xcorrpic
global nX nY nZ
global cellf
global pathname

disp('Getting neuropil data...')

for i=1:size(cellmask,3)
    [column,row]=find(cellmask(:,:,i)==1);
    centroidx(i)=mean(row);     %cell location centroid, x
    centroidy(i)=mean(column);  %cell location centroid, y
end
totalcellmask=(sum(cellmask,3)>0);  %spatial mask for all cells
            
for i=1:size(cellmask,3)
    %for each cell, find the associated neighboring neuropil area, an annulus defined by r1 and r2
    r1=2*sqrt(sum(sum(cellmask(:,:,i)))/pi);  %inner ring, 2x cell radius
    r2=3*sqrt(sum(sum(cellmask(:,:,i)))/pi);  %outer ring, 3x cell radius
    neuropilmask(:,:,i)=zeros(size(cellmask(:,:,i)));
    for x=1:nX
        for y=1:nY
            if (((x-centroidx(i))^2 + (y-centroidy(i))^2) > r1^2) && (((x-centroidx(i))^2 + (y-centroidy(i))^2) < r2^2)
                if x>=1 && x<=nX && y>=1 && y<=nY
                    neuropilmask(y,x,i)=1;      %the spatial mask for neuropil
                end
            end
        end
    end
    %remove from neuropil spatial mask those pixels that are part of ROIs of neurons
    neuropilmask(:,:,i)=neuropilmask(:,:,i) & ~(totalcellmask);
    clear r1 r2;
    
    %find the fluorescence of the neuropil area
    for t=1:size(pic,3)
        temp=pic(:,:,t);
        neuropilf(t)=sum(temp(neuropilmask(:,:,i)==1));
    end

    filenameheader=get(handles.edit8,'String');
    subtractmask=neuropilmask(:,:,i);
    neuropilf=neuropilf/sum(sum(subtractmask)); %per-pixel fluorescence
    temp = sprintf('%03d',i);
    save(strcat(pathname,filenameheader,'_cell',temp,'.mat'),'neuropilf','subtractmask','-append');
    
    disp(['... neuropil for cell ' int2str(i)]);
end
totalneuropilmask=(sum(neuropilmask,3)>0);  %spatial mask for all neuropil

zprojWhiteLev=str2num(get(handles.edit11,'String'));
roi_pic=round(zprojWhiteLev*253)*rightpanel_pic/nanmax(rightpanel_pic(:));   %because 254 set to yellow (for old ROI); 255 is set to red (for new ROI)
roi_pic(roi_pic>252)=252;
roi_pic(savedroi_pic)=254;
roi_pic(totalneuropilmask)=255;
cla(handles.axes1);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes1,graymap);
image(roi_pic,'Parent',handles.axes1);
axis(handles.axes1, 'off', 'equal');

% --- Executes on button press in loadRed.
function loadRed_Callback(hObject, eventdata, handles)
% hObject    handle to loadRed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global pic savedroi_pic red_pic shift_pic
global leftpanel_pic

[c2_fname, c2_pathname]=uigetfile('*.tif','Load red channel z-projection');
red_pic=double(imread([c2_pathname c2_fname],'tif'));  %load overlay z-stack
shift = [-str2num(get(handles.edit16,'String')) str2num(get(handles.edit15,'String'))];
shift_pic = circshift(red_pic, shift);

set(handles.radiobutton13,'Value',1);

leftpanel_pic=shift_pic;
leftWhiteLev=str2num(get(handles.edit10,'String'));
temp_pic=uint16(round(leftWhiteLev*253)*leftpanel_pic/nanmax(leftpanel_pic(:))); 
temp_pic(temp_pic>252)=252;
temp_pic(savedroi_pic)=254;
cla(handles.axes4);
graymap=[linspace(0,1,253) 1 1; linspace(0,1,253) 1 0; linspace(0,1,253) 0 0]';
colormap(handles.axes4,graymap);
image(temp_pic,'Parent',handles.axes4);
axis(handles.axes4, 'off', 'equal');

%------------------------------------------------------
function edit1_Callback(hObject, eventdata, handles)

function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox1_Callback(hObject, eventdata, handles)

function edit2_Callback(hObject, eventdata, handles)

function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)

function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_Callback(hObject, eventdata, handles)

function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)

function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)

function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit7_Callback(hObject, eventdata, handles)

function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox3_Callback(hObject, eventdata, handles)

function edit8_Callback(hObject, eventdata, handles)

function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox4_Callback(hObject, eventdata, handles)

function checkbox5_Callback(hObject, eventdata, handles)






function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function AuxProjection_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AuxProjection_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8




