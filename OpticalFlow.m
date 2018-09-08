function [ X,Y,Vx,Vy,mov,varargout ] = OpticalFlow( MovieName,BinaryMask,BoxSize,BlurSTD,ArrowSize,scale,dt,varargin )

%AdvectiveImageVelocimetry.m Uses the Advective Image Velocimetry (AIV)
%method developed by Dhruv K. Vig, Alex E. Hamby and Charles W. Wolgemuth
%to compute the velocity field (Vx,Vy) from a time sequence of images stored in
%the TIFF or AVI file MovieName. The velocity field is output on a coarser grid
%than the original images.  This grid is stored in the vectors (X,Y). The user
%can also define a region of interest, which is input as a BinaryMask.

%  License: This code is free to use, provided
%           any publications resulting from the use of this code
%           reference the original code/authors.
%   
%%  Copyright 2015 Arizona Board of Regents on behalf of The University of Arizona

%%  Author:  Dhruv Kumar Vig (dvig89@gmail.com) and Charles Wolgemuth(wolg@email.arizona.edu)
%%  Date:    1/2015
%
%  Please notify the author of any bugs, and contribute any
%  modifications or bug fixes back to the original author.

warning off

%% Identify the format of the Image Sequences as either AVI or TIFF

    num = strfind(MovieName,'.tif');
    % num = strfind(MovieName,'.');
    if strcmpi(MovieName(num+1:end),'avi') == true
    
        Mov = VideoReader(MovieName); 
        NumFrames = Mov.NumberOfFrames;
  
        W = Mov.Height;
        L = Mov.Width;
     
        FullFrame = read(Mov,1);
    
    elseif strcmpi(MovieName(num+1:end),'tif') == true
        
        Info = imfinfo(MovieName);
        NumFrames = numel(Info);  
   
        W = Info(1).Height;
        L = Info(1).Width;
            
        FullFrame = imread(MovieName,1);
        
    else
        
        error('Incorrect image sequence file. Should be .avi or .tif')

    end
    
%% Define Additional Filter for Rotational AIV

   narginchk(7,8)
   
   if isempty(varargin) == true || strcmpi(varargin{1},'none') == true
       
   elseif strcmpi(varargin{1},'Rotation') == true
       
         MaxR = floor(BoxSize./2);
         Range = (-MaxR:MaxR);
         [Xfilt,Yfilt] = meshgrid(Range,Range);
        
   elseif strcmpi(varargin{1},'React') == true
       
   else
       
         error('Incorrect final input. Should be either [], ''none'', ''React'' or ''Rotation''')    
       
   end
   
%% Identify whether a BinaryMask is neccessary
    
   if isempty(BinaryMask)
       
       Mask = true(W,L);
       
   else
       
      Binarynum = strfind(BinaryMask,'.');
      Mask = false(W,L);
      
   end

%% AIV Filters 

  BlurSize = ceil(3.5.*BlurSTD) + mod(ceil(3.5.*BlurSTD),2) + 1;
  h = fspecial('gaussian',[BlurSize BlurSize],BlurSTD);         
  h2 = ones(BoxSize,BoxSize);
   
%% Subtract the background by fitting the image to a cubic function. See Supplemental Notes __

  FrameOld = double(FullFrame(:,:,1));
  
  [X,Y] = meshgrid((1:L),(1:W));
  [Mat] = BackgroundMatrix(X,Y,W,L);
  
  [FrameOld] = BackgroundSubtract(FrameOld,Mat,X,Y);
  
  FrameOld = ImageFilter(FrameOld,h,'symmetric'); 
    
%% Define SubMask otherwise the velocites will be calculated at every pixel in the image.

  SubMask1 = false(W,L);
  SubMask = SubMask1;

  SubMask1(11:ArrowSize:W-10,:) = true;
  SubMask(:,11:ArrowSize:L-10) = true;

  SubMask = SubMask & SubMask1;
  
  Xsub = X(SubMask);
  Ysub = Y(SubMask);

%% Iterate through Image Pairs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
for t = 1:NumFrames-1    
    
%% read in Frame at time t+1 and subtract background as above.
   
    if strcmpi(MovieName(num+1:end),'avi') == true
   
        FullFrame = read(Mov,t);
    
    elseif strcmpi(MovieName(num+1:end),'tif') == true
                 
        FullFrame = imread(MovieName,t);

    end

    FrameNew = double(FullFrame(:,:,1));

    [FrameNew] = BackgroundSubtract(FrameNew,Mat,X,Y);
    
    FrameNew = ImageFilter(FrameNew,h,'symmetric');
    
%% Compute image difference and average

    dIdt = FrameNew - FrameOld;
    Iavg = 0.5.*(FrameNew+FrameOld);

%% Compute image gradients

    Ix = 0.5.*( cshift2(Iavg,[0 -1]) - cshift2(Iavg,[0 1]) );
    Iy = 0.5.*( cshift2(Iavg,[-1 0]) - cshift2(Iavg,[1 0]) );
  
%% Load Binary ROI Stack

    if strcmpi(BinaryMask(Binarynum+1:end),'avi') == true
        
        Mov = VideoReader(BinaryMask); 
        Mask = read(Mov,t);
    
    elseif strcmpi(BinaryMask(Binarynum+1:end),'tif') == true
          
        Mask = imread(BinaryMask,t);
    end
 
      SumMask = Mask & SubMask;

%% Compute the velocity

   if isempty(varargin) == true || strcmpi(varargin{1},'none') == true
       
    % standard AIV

       [Vx,Vy] = AIVVelocity(Ix,Iy,dIdt,h2,Mask,W,L);
    
   elseif strcmpi(varargin{1},'Rotation') == true
       
    % AIV with Rotation       
     
       [Vx,Vy,Om] = AIVVelocitywRotation(Ix,Iy,dIdt,h2,Mask,Xfilt,Yfilt,W,L);
    
        Omout(:,t) = Om(SubMask);
         
        opengl('software');
    
        f = figure;
        set(f,'Renderer','zbuffer');
        set(f,'color','w');
    
        surf(X,Y,Om-1.05.*max(Om(:)));shading interp;view([0 90]);axis([0 length(FullFrame) 0 length(FullFrame)]); axis off; grid off;
        hold on
        quiver(X(SumMask),Y(SumMask),7.*Vx(SumMask),7.*Vy(SumMask),0,'k','LineWidth',1);
    
        OmMov(:,t) = getframe;
        close all 
       
   elseif strcmpi(varargin{1},'React') == true
       
    % AIV with a Reaction term

       [Vx,Vy,Gamma] = AIVVelocitywReaction(Ix,Iy,dIdt,h2,Mask,W,L);
    
          Gammaout(:,t) = Gamma(SubMask);
    
   end

% Overlay the velocity field on the image. 

    clf;
    subplot(1,2,1); % add by Wanyu
    imshow(FullFrame,[min(FullFrame(:)),quantile(FullFrame(:),0.99)]); % original imshow(FullFrame,[])
    hold on
    quiver(X(SumMask),Y(SumMask),3.*Vx(SumMask),3.*Vy(SumMask),0,'g','LineWidth',1)  

    Vx(~Mask) = NaN;
    Vy(~Mask) = NaN;

    VxSub(:,t) = (Vx(SubMask).*scale)./dt;
    VySub(:,t) = (Vy(SubMask).*scale)./dt;
    
    subplot(1,2,2); % add by Wanyu
    imagesc(FullFrame); colorbar;  % add by Wanyu
    
    mov(:,t) = getframe;
    
    FrameOld = FrameNew;
    
end

Vx = VxSub;
Vy = VySub;

X = X(SubMask);
Y = Y(SubMask);

   if isempty(varargin) == true || strcmpi(varargin{1},'none') == true
       
             varargout{1} = [];
       
   elseif strcmpi(varargin{1},'Rotation') == true
       
      varargout{1} = Omout;
      varargout{2} = OmMov;
       
   elseif strcmpi(varargin{1},'React') == true
       
       varargout{1} = Gammaout;
       
   end
    
end

%% Subroutines %%

function [Mat] = BackgroundMatrix(X,Y,W,L)

    Xsum = sum(X(:));
    Ysum = sum(Y(:));
    X2 = sum(X(:).^2);
    XY = sum(X(:).*Y(:));
    Y2 = sum(Y(:).^2);
    X3 = sum(X(:).^3);
    X2Y = sum(X(:).^2.*Y(:));
    XY2 = sum(X(:).*Y(:).^2);
    Y3 = sum(Y(:).^3);
    X4 = sum(X(:).^4);
    X3Y = sum(X(:).^3.*Y(:));
    X2Y2 = sum(X(:).^2.*Y(:).^2);
    XY3 = sum(X(:).*Y(:).^3);
    Y4 = sum(Y(:).^4);
    X5 = sum(X(:).^5);
    X4Y = sum(X(:).^4.*Y(:));
    X3Y2 = sum(X(:).^3.*Y(:).^2);
    X2Y3 = sum(X(:).^2.*Y(:).^3);
    XY4 = sum(X(:).*Y(:).^4);
    Y5 = sum(Y(:).^5);
    X6 = sum(X(:).^6);
    X5Y = sum(X(:).^5.*Y(:));
    X4Y2 = sum(X(:).^4.*Y(:).^2);
    X3Y3 = sum(X(:).^3.*Y(:).^3);
    X2Y4 = sum(X(:).^2.*Y(:).^4);
    XY5 = sum(X(:).*Y(:).^5);
    Y6 = sum(Y(:).^6);

    Mat = [ W.*L  Xsum  Ysum  X2     XY     Y2     X3      X2Y     XY2      Y3;
            Xsum           X2    XY    X3     X2Y    XY2    X4      X3Y     X2Y2    XY3;
            Ysum           XY    Y2    X2Y    XY2    Y3     X3Y     X2Y2    XY3      Y4;
            X2             X3    X2Y   X4     X3Y    X2Y2   X5      X4Y     X3Y2   X2Y3;
            XY             X2Y   XY2   X3Y    X2Y2   XY3    X4Y     X3Y2    X2Y3    XY4;
            Y2             XY2   Y3    X2Y2   XY3    Y4     X3Y2    X2Y3    XY4      Y5;
            X3             X4    X3Y   X5     X4Y    X3Y2   X6      X5Y     X4Y2   X3Y3;
            X2Y            X3Y   X2Y2  X4Y    X3Y2   X2Y3   X5Y     X4Y2    X3Y3   X2Y4;
            XY2            X2Y2  XY3   X3Y2   X2Y3   XY4    X4Y2    X3Y3    X2Y4    XY5;
            Y3             XY3   Y4    X2Y3   XY4    Y5     X3Y3    X2Y4    XY5      Y6];

end

function [Frame] = BackgroundSubtract(Frame,Mat,X,Y)

 SourceBG = [ sum(Frame(:));
                 sum(X(:).*Frame(:));
                 sum(Y(:).*Frame(:));
                 sum(X(:).^2.*Frame(:));
                 sum(X(:).*Y(:).*Frame(:));
                 sum(Y(:).^2.*Frame(:));
                 sum(X(:).^3.*Frame(:));
                 sum(X(:).^2.*Y(:).*Frame(:));
                 sum(X(:).*Y(:).^2.*Frame(:));
                 sum(Y(:).^3.*Frame(:))           
                 ];

    Answer = Mat\SourceBG;

    BG =   Answer(1) + Answer(2).*X + Answer(3).*Y + Answer(4).*X.^2 ...
          + Answer(5).*X.*Y + Answer(6).*Y.^2 + Answer(7).*X.^3 ...
          + Answer(8).*X.^2.*Y + Answer(9).*X.*Y.^2 + Answer(10).*Y.^3;

    Frame = Frame - BG;
    
end

function [Vx,Vy] = AIVVelocity(Ix,Iy,dIdt,h2,Mask,W,L)

%% initialize velocity field

    Vx = zeros(W,L);
    Vy = Vx;

%% compute coefficients for LSM
    
    A = ImageFilter(Ix.*Ix,h2,'symmetric');
    B = ImageFilter(Ix.*Iy,h2,'symmetric');
    C = ImageFilter(Iy.*Iy,h2,'symmetric');
    
    ItIx = ImageFilter(dIdt.*Ix,h2,'symmetric');
    ItIy = ImageFilter(dIdt.*Iy,h2,'symmetric');

%% determine velocity
    
    Vx(Mask) = (-C(Mask).*ItIx(Mask) + B(Mask).*ItIy(Mask))./(A(Mask).*C(Mask)-B(Mask).^2);
    Vy(Mask) = (-A(Mask).*ItIy(Mask) + B(Mask).*ItIx(Mask))./(A(Mask).*C(Mask)-B(Mask).^2);
            
end

function [Vx,Vy,Om] = AIVVelocitywRotation(Ix,Iy,dIdt,h2,Mask,Xfilt,Yfilt,W,L)

%% initialize velocity field

    Vx = zeros(W,L);
    Vy = Vx;
    Om = Vx;

%% compute coefficients for LSM
    
    A = ImageFilter(Ix.*Ix,h2,'symmetric');
    B = ImageFilter(Ix.*Iy,h2,'symmetric');
    C = ImageFilter(Ix.*Iy,Xfilt,'symmetric') - ImageFilter(Ix.^2,Yfilt,'symmetric');
    D = ImageFilter(Iy.*Iy,h2,'symmetric');
    E = ImageFilter(Iy.*Iy,Xfilt,'symmetric') - ImageFilter(Ix.*Iy,Yfilt,'symmetric');
    F = ImageFilter(Iy.^2,Xfilt.^2,'symmetric') + ImageFilter(Ix.^2,Yfilt.^2,'symmetric') - 2.*ImageFilter(Ix.*Iy,Xfilt.*Yfilt,'symmetric');
    
    ItIx = ImageFilter(dIdt.*Ix,h2,'symmetric');
    ItIy = ImageFilter(dIdt.*Iy,h2,'symmetric');  
    ItRot = ImageFilter(dIdt.*Iy,Xfilt,'symmetric') - ImageFilter(dIdt.*Ix,Yfilt,'symmetric');

%% determine velocity
    
    Denom = -(F(Mask).*B(Mask).^2 - 2.*B(Mask).*C(Mask).*E(Mask) + D(Mask).*C(Mask).^2 + A(Mask).*E(Mask).^2 - A(Mask).*D(Mask).*F(Mask));

    Vx(Mask) = ((E(Mask).^2 - D(Mask).*F(Mask)).*ItIx(Mask) + (B(Mask).*F(Mask)-C(Mask).*E(Mask)).*ItIy(Mask) + (C(Mask).*D(Mask) - B(Mask).*E(Mask)).*ItRot(Mask)) ./Denom;
    Vy(Mask) = ((B(Mask).*F(Mask) - C(Mask).*E(Mask)).*ItIx(Mask) + (C(Mask).^2 - A(Mask).*F(Mask)).*ItIy(Mask) + (A(Mask).*E(Mask) - B(Mask).*C(Mask)).*ItRot(Mask)) ./Denom;
    Om(Mask) = ((C(Mask).*D(Mask) - B(Mask).*E(Mask)).*ItIx(Mask) + (A(Mask).*E(Mask) - B(Mask).*C(Mask)).*ItIy(Mask) + (B(Mask).^2 - A(Mask).*D(Mask)).*ItRot(Mask)) ./Denom;
     
end

function [Vx,Vy,Gamma] = AIVVelocitywReaction(Ix,Iy,dIdt,h2,Mask,W,L)

%% initialize velocity field

    Vx = zeros(W,L);
    Vy = Vx;
    Gamma = Vx;

%% compute coefficients for LSM
    
    A = ImageFilter(Ix.*Ix,h2,'symmetric');
    B = ImageFilter(Ix.*Iy,h2,'symmetric');
    C = ImageFilter(Ix,h2,'symmetric');
    D = ImageFilter(Iy.*Iy,h2,'symmetric');
    E = ImageFilter(Iy,h2,'symmetric');
    F = ImageFilter(ones(size(Ix)),h2,'symmetric');
    
    ItIx = ImageFilter(dIdt.*Ix,h2,'symmetric');
    ItIy = ImageFilter(dIdt.*Iy,h2,'symmetric');
    It = ImageFilter(dIdt,h2,'symmetric');

%% determine velocity
    
    Denom = -(F(Mask).*B(Mask).^2 - 2.*B(Mask).*C(Mask).*E(Mask) + D(Mask).*C(Mask).^2 + A(Mask).*E(Mask).^2 - A(Mask).*D(Mask).*F(Mask));    

    Vx(Mask) = ((E(Mask).^2 - D(Mask).*F(Mask)).*ItIx(Mask) + (B(Mask).*F(Mask)-C(Mask).*E(Mask)).*ItIy(Mask) + (C(Mask).*D(Mask) - B(Mask).*E(Mask)).*It(Mask))./Denom;
    Vy(Mask) = ((B(Mask).*F(Mask) - C(Mask).*E(Mask)).*ItIx(Mask) + (C(Mask).^2 - A(Mask).*F(Mask)).*ItIy(Mask) + (A(Mask).*E(Mask) - B(Mask).*C(Mask)).*It(Mask))./Denom;
    Gamma(Mask) = ((C(Mask).*D(Mask) - B(Mask).*E(Mask)).*ItIx(Mask)+ (A(Mask).*E(Mask) - B(Mask).*C(Mask)).*ItIy(Mask) + (B(Mask).^2 - A(Mask).*D(Mask)).*It(Mask))./Denom;
    
end

function [B] = cshift2(A,compass)

[W,L] = size(A);

if compass(1)<0
    
    Iy = [-compass(1)+1:W 1:-compass(1)];
    
elseif compass(1)>0
    
    Iy = [W-compass(1)+1:W 1:W-compass(1)];
    
else
    
    Iy = (1:W);
    
end

if compass(2)<0
    
    Ix = [-compass(2)+1:L 1:-compass(2)];
    
elseif compass(2)>0
    
    Ix = [L-compass(2)+1:L 1:L-compass(2)];
    
else
    
    Ix = (1:L);
    
end

B = A(Iy,Ix);

end

function [ ImageOut ] = ImageFilter( ImageIn, kernel, varargin )
%ImageFilter Filters the matrix ImageIn using the kernel.  The user can
%also define how to pad the edges.  The choices are 'symmetric' (default),
%'circular', or a constant value.

% determine kernel and image sizes

[ Lk, Wk ] = size(kernel);
[ L, W ] = size(ImageIn);

% test input arguments

Num = nargin;

if Num == 2;
    varargin = 'symmetric';
end

if Num > 3
    error('myApp:argChk','Too many input arguments')
end

% pad Image matrix using user-defined boundary conditions

if varargin{1} == 'symmetric'
    
    ImagePad = [ ImageIn(Lk+1:-1:2,Wk+1:-1:2)      ImageIn(Lk+1:-1:2,:)      ImageIn(Lk+1:-1:2,W-1:-1:W-1-Wk);
                 ImageIn(:,Wk+1:-1:2)              ImageIn                   ImageIn(:,W-1:-1:W-1-Wk);
                 ImageIn(L-1:-1:L-1-Lk,Wk+1:-1:2)  ImageIn(L-1:-1:L-1-Lk,:)  ImageIn(L-1:-1:L-1-Lk,W-1:-1:W-1-Wk) ];
             
elseif varargin{1} == 'circular'
    
    ImagePad = [ ImageIn(L-Lk:L,W-Wk:W)            ImageIn(L-Lk:L,:)        ImageIn(L-Lk:L,1:Wk);
                 ImageIn(:,W-Wk:W)                 ImageIn                  ImageIn(:,1:Wk);
                 ImageIn(1:Lk,W-Wk:W)              ImageIn(1:Lk,:)          ImageIn(1:Lk,1:Wk) ];
             
elseif isfloat(varargin{1}) == true
    
    c = varargin{1};
    
    ImagePad = [ c.*ones(Lk,Wk)       c.*ones(Lk,W)    c.*ones(Lk,Wk);
                 c.*ones(L,Wk)        ImageIn          c.*ones(L,Wk);
                 c.*ones(Lk,Wk)       c.*ones(Lk,W)    c.*ones(Lk,Wk) ];
             
else
    error('myApp:argChk','Incorrect format for third input argument')
end

%% test kernel for separability

[S, Hcol, Hrow] = isfilterseparable(kernel);

% filter Image

if S == true
    
    ImageOut = filter2(Hrow,ImagePad,'same');
    ImageOut = filter2(Hcol,ImageOut,'same');
    
else

    ImageOut = filter2(kernel,ImagePad,'same');

end

ImageOut = ImageOut(Lk+1:L+Lk,Wk+1:W+Wk);

end

function [S, HCOL, HROW]  = isfilterseparable(H)
% ISFILTERSEPARABLE  Check filter separability. 
% S = ISFILTERSEPARABLE(H) takes in the filter kernel H and returns true
% if filter is separable, otherwise it returns false.
%
% [S, HCOL, HROW]  = ISFILTERSEPARABLE(H) optionally returns the vertical
% coefficients HCOL  and horizontal coefficients HROW, if the filter is
% separable, otherwise HCOL and HROW are empty. 
%
% Class Support
% -------------
% H can be logical or numeric, 2-D, and nonsparse. 
% S is logical, HCOL and HROW are the same class as H if H is float 
% otherwise they are double. 

%   Copyright 2003-2005 The MathWorks, Inc. 

S = false;
if (~isa(H,'float')),  H = double(H); end
if all(isfinite(H(:)))
  % Check rank (separability) of H
  [u,s,v] = svd(H);
  s = diag(s);
  tol = length(H) * eps(max(s));
  rank = sum(s > tol);   
  S = (rank ==1);
end
HCOL = [];
HROW = [];
if S
    HCOL = u(:,1) * sqrt(s(1));
    HROW = conj(v(:,1)) * sqrt(s(1));
end
HROW = HROW.'; 
end