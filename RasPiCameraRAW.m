% This function saves image from Raspberry Pi RGB camera and outputs separate 
% 4 color components as [1232x1640x4] 3D matrix, 
% imOut(:,:,1) is red component, imOut(:,:,2) green 1, imOut(:,:,3) green 2,
% and imOut(:,:,4) blue, see how Bayer filter works.
% Saves two temporary image files RAW and 16bit TIFF.
% Uses external "raspistill" software to capture RAW and custom "dcraw" to 
% convert raw into TIFF, finally TIFF is loaded into Octave.
% Assumes that custom version of "dcraw" which can handle peculiar Raspberry Pi 
% camera RAWs is installes, for more info see Raspberry Pi forum.
%
% input:
% fileName: filename for temporary image files,
% ISO value: detector sensitivity, 100 is slowest and good default,
% expTime: exposure time in micro seconds
%
% Henri Partanen 21.11.2018  

function [imOut]=RasPiCameraRAW(fileName,ISOvalue,expTime)


RAWheader='.jpg';
IMAGEheader='.tiff';

saturatedValue=4603; % experimentally found out that saturated values of
                     % Raspberry Pi camera V2 have byte value of 4603

% command line parameters to run dcraw to convert RAW images,
% -D = no Bayer filter removal, -6 = 16-bit output image, 
% -W = no automatic brightening, -T = output as tiff
dcrawString='dcraw -D -6 -W -T ';

% command line parameters to capture image with Raspbery Pi camera,
% -f = full scree preview, -ISO = detector sensitivity, 
% -ss = exposure time (in micro seconds?), 
% -awb off = automatic white balance off,
% -t = preview time, -o = output filename
% As and ad hoc raspistill saves RAW data into metadata of jpg file,
% for more information see Raspberry Pi forums.
captureString=['raspistill -f -ISO ',num2str(ISOvalue),' -ss ',num2str(expTime),' -awb off -awbg 1.5,1.5 -t 5 -r -o '];

% save image file with external raspistill software
system([captureString,fileName,RAWheader]);
disp('picture taken')

%% convert just saved RAW file into tiff using external dcraw software
%% needs a modified version of dcraw that can handle these strange RAW jpgs.
system([dcrawString,fileName,RAWheader]);
%
%% open just converted tiff image file
im=double(imread([fileName,IMAGEheader]));
% break it into color components
imB =im(1:2:end, 1:2:end);  % blue
imG1=im(1:2:end, 2:2:end);  % green 1 (two green values because Bayer filter)
imG2=im(2:2:end, 1:2:end);  % green 2
imR =im(2:2:end, 2:2:end);  % red

% insert color components into 3D matrix
imOut(:,:,1)=imR;
imOut(:,:,2)=imG1;
imOut(:,:,3)=imG2;
imOut(:,:,4)=imB;

% display how many saturated pixels image had
maxValue=max(im(:));
maxPercents=maxValue/saturatedValue*100;
Nsaturated=sum(im(:)>=saturatedValue);
disp(['Image has ',num2str(Nsaturated),' saturated pixels']);
disp(['max value is ',num2str(maxPercents),' percents of saturated']);