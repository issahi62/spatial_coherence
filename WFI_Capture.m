% This function captures 4 images from Raspberry Pi camera and saves the into
% Matlab compatible binary data file.
% Needs "RasPiCameraRAW" sub function.
%
% Input:
% saveName: filename as text string without file extension
% expTime: expposure time in micro seconds
% AOI: area of interes, cropped area of image as 4 componen vectorize
%      [x minimum, x maximum, y minimum, y maximum]
%      maximum image size is 1232 x 1640 pixels [y times x].
% whichColor: color component of RGB color camera saved 
%             1 = red, 2 = green1, 3 = gree2, 4 = blue
% plotImages: 1 = plots monochrome image of saves data (slower), 0 = doesn't plot
%
% Saved data file contains variables: 'ICut','I1Cut','I2Cut','IdCut','expTime'
%
% For example "WFI_Capture('testFile', 50e3, [1,1640,400,950], 2, 1,0)"
% saves data into "testfile.mat" file, exposure time 50 milliseconds, 
% image area [1,1640,400,950], saves green color component, and plots the image.
%
% Henri Partanen 21.11.2018


function WFI_Capture(saveName, expTime, AOI, whichColor, plotImages, testmode)


close all

more off

%plotImages=1;

fileName='temppi'
ISOvalue=100;

%whichColor=2; %1=R, 2=G,  4=B

xArea=[AOI(1):AOI(2)];
yArea=[AOI(3):AOI(4)];

saturatedValue=4603;


  disp('Open both arms')
  pause
    imOut=RasPiCameraRAW(fileName,ISOvalue,expTime);
    I=imOut(:,:,whichColor);

if(testmode==0)
  disp('Block arm 1')
  pause
    imOut=RasPiCameraRAW(fileName,ISOvalue,expTime);
    I1=imOut(:,:,whichColor);

  disp('Block arm 2')
  pause
    imOut=RasPiCameraRAW(fileName,ISOvalue,expTime);
    I2=imOut(:,:,whichColor);

  disp('Block everything')
  pause
    imOut=RasPiCameraRAW(fileName,ISOvalue,expTime);
    Id=imOut(:,:,whichColor);
end

  ICut=I(yArea,xArea); 

if(testmode==0) 
  I2Cut=I2(yArea,xArea);
  I1Cut=I1(yArea,xArea);
  
  IdCut=Id(yArea,xArea);
  
    save([saveName,'.mat'],'ICut','I1Cut','I2Cut','IdCut','expTime','-mat-binary')
end 

  ICutPlot=ICut;
  ICutPlot(ICutPlot==saturatedValue)=0;
 
  


if(plotImages==1)
    figure
    imagesc(ICutPlot)
    caxis([1,saturatedValue])
    
    colorbar
    colormap hot
end
