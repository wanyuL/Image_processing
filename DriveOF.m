function DriveOF (MovieName,ROIName,varargin)
%% This program is a driver for AdvectiveImageVelocimetry.m which computes the local velocites filed that advects the image intensity from one frame to the next. 
%% Specifically we solve dI/dt = -v*grad(I) for v using a least squares method. 

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

%% Input Parameters 

% ROIName is the identifying name for the Binary Image Stacks that can be used to
% identify regoins-of-interests (ROIs) within an image. The file type can
% be either AVI or TIF.  The format for the
% name is FileName.avi or FileName.tif.  For example, if the user has
% three image stacks that they want to process, and the ROIs are stored in
% files named ROIMask1.tif, ROIMask2.tif, and ROIMask3.tif, then the
% ROIName should be ROIMask.tif.  The identifying label for the files,
% i.e., the integers 1, 2 and 3 in the example, can have any form.
%
% MovieName is the name of the image sequence file that will be read in.
% The file type can be either AVI or TIF, and the name follows the same
% format as was described above for the ROIName.  The identifying label
% needs to be consistent betwen the MovieNames and the ROINames.

% ResultsName is the name of the outputed workspace file that will be
% saved.  It includes the X and Y positions of the velocity vectors, as well as
% the x-component and y-components of the velocities. Vx and Vy are
% outputed as matrices, where the rows correspond with spatial location
% and the columns correspond to the frame number. 

%% AIV Parameters

% BoxSize needs to be larger than the size of a single cell. 

% BlurSTD sets the standard deviation for the Guassian blur. Should be approximately 
%         half the maximum velocity between two frames in pixels. 

% SmoothSize sets the size of a Guassian Kernel. 

% ArrowSize sets the size of the SubMask region, otherwise the velocities
% will be calculated at every pixel. 

%% Driver Parameters 

TotalMovies = 1;             %total number of movies the user wishes to analyze. This value has to be at least 1.

scale = 1;              %1 pixel = (blank) microns
dt = 1;                      %time between frames

%% AIV Parameters

BoxSize = 21;               %should be set large enough that each box contains at least one feature

BlurSTD = 4;

ArrowSize = 21;

%% Find File Format for Movies in Directory

num1 = strfind(MovieName,'.');

if strcmp(MovieName(num1+1:end),'avi') == true
    
    ImageSequenceFormat = '.avi';
    
elseif strcmp(MovieName(num1+1:end),'tif') == true
    
    ImageSequenceFormat = '.tif';
    
else
    
    error('Oops! You have made an error. Please check file format and try again.')
    
end

%% Find File Format for Binary ROI Masks in Directory

num2 = strfind(ROIName,'.');

if strcmp(ROIName(num2+1:end),'avi') == true
    
    BinaryFormat = '.avi';
    
elseif strcmp(ROIName(num2+1:end),'tif') == true
    
    BinaryFormat = '.tif';
    
elseif isempty(ROIName) == true
    
    BinaryMask = [];
        
else 
    
    error('Oops! You have made an error. Please check file format and try again.')
        
end

%% Find Prefix

prefix = MovieName(1:num1-1);

ImageFiles = dir(strcat(prefix,'*',ImageSequenceFormat));

for m = 1:length(ImageFiles)
    
    num = strfind(ImageFiles(m).name,'.');
    
    if num > 7
        q(m) = ~strcmp(ImageFiles(m).name(num-7:num-1),'Tracked');
    else
        q(m) = true;       
    end
    
end

ImageFiles = ImageFiles(q);

if isempty(ROIName) == true
        
    for m = 1:length(ImageFiles)
    ROIFiles(m).name = [];
    end
    
else
    
    ROIFiles = dir(strcat(ROIName(1:num2-1),'*',BinaryFormat));

end

%% Run AIV and Loop over movies

for i = 1:TotalMovies
            
    MovName = ImageFiles(i).name;
    BinaryMask = ROIFiles(i).name;
    
    num = strfind(MovName,'.');
    prefix = MovName(1:num-1);
    
    ResultsName = strcat(prefix,'Results','.mat');
    TrackedMovie = strcat(prefix,'Tracked','.avi'); 
    
if isempty(varargin) == true || strcmpi(varargin{1},'none') == true
    
    [ X,Y,Vx,Vy,mov] = OpticalFlow( MovieName,BinaryMask,BoxSize,BlurSTD,ArrowSize,scale,dt,'none' );
     
    save(ResultsName,'X','Y','Vx','Vy','mov')
    movie2avi(mov,TrackedMovie,'compression','none');
    
    elseif strcmp(varargin{1},'Rotation') == true
    
    TrackedMovieOM = strcat(prefix,'OmTracked','.avi');   
    
   [ X,Y,Vx,Vy,mov,Om,OmMov ] = OpticalFlow( MovieName,BinaryMask,BoxSize,BlurSTD,ArrowSize,scale,dt,'Rotation' );
       
    save(ResultsName,'X','Y','Vx','Vy','Om','OmMov','mov')
    movie2avi(mov,TrackedMovie,'compression','none');
    movie2avi(OmMov,TrackedMovieOM,'compression','none');
    
    elseif strcmp(varargin{1},'React')  == true
    
    [ X,Y,Vx,Vy,mov,Gamma ] = OpticalFlow( MovieName,BinaryMask,BoxSize,BlurSTD,ArrowSize,scale,dt,'React' );
   
    save(ResultsName,'X','Y','Vx','Vy','Gamma','mov');
    movie2avi(mov,TrackedMovie,'compression','none');
    
    else
    
     error('Oops! You have made an error. Please check your inputs.')
        
end

end
