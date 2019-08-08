%%
% Spatial_Coherence_Analysis 
%*******************************************************
% Developed by Ibrahim Issah & Srijoyee Datta & Mustafa Aboulsaad5 
%****************************************************************
%%


close all;
%%
% section for loading the files 

 


myFolder = '/Users/kobbyTilly/Desktop/Spatial_coherence_IMS'; 
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
%%
% Used constants in the codes
% radius to draw the circle on each fftimage (Assumend value was around 36.8) 
r=40;

% identified coefficients with getpts command
px = 6.808024193548388e+02;
py = 6.285981996726675e+02;

% radius for removing unwanted section of the fftimage
kobbys_circle_radius = 30; 

% incrementors in the circshift command : still trying shiftdim.

%% incrementors in the xc ( x coordinates)
ShifterX=0;
AnotherSrijoyee = [];
AnotherSrijoyee2 = []; 
% incrementors in the y-axis or yc
ShifterY=+0;
%%
vect = 1; 
filePattern = fullfile(myFolder, 'He_Ne*.mat');
matFiles = dir(filePattern);
[~, ind] = sort([matFiles.datenum]);
matFiles = matFiles(ind);

for k = 1:length(matFiles)
  matFilename = fullfile(myFolder, matFiles(k).name);
  matData = load(matFilename);
  
  hasField = isfield(matData, 'expTime');
  if ~hasField
    % Alert user in popup window.
    warningMessage = sprintf('Data is not in %s\n', matFilename);
    uiwait(warndlg(warningMessage));
    % If the file is not there then skip to the next iteration
    continue; % Go to the next iteration.
  end
  
  I1Cut = matData.I1Cut; 
  I2Cut = matData.I2Cut; 
  IdCut = matData.IdCut; 
  ICut = matData.ICut; 

    %Normaliztion section for each file 
  newI1Cut = I1Cut - IdCut; 
  newI2Cut = I2Cut - IdCut; 
  newICut  =  ICut - IdCut; 
  Snorm = real((newICut - newI1Cut - newI2Cut)./(2.*sqrt(newI1Cut.*newI2Cut)));
  Snorm(newI1Cut==0| newI2Cut==0) = 0; 
  % Setting normalized values in the range of -1 to 1 using the find command 

  %% *************************************
% CenterSnorm(find(Snorm<-1)) = -1;
% CenterSnorm(find(isnan(Snorm))) = 0;
% CenterSnorm(find(Snorm>1)) = 1;
%***************************************

%%
 % setting the normalized .mat values to range of -1 to 1

 Snorm(Snorm<-1)=-1;
 Snorm(isnan(Snorm))=0;
 Snorm(Snorm>1)=1;
%***********************************************
% Fourier analysis of the normalized image 
%***********************************************

forSnorm = fftshift(fft2(Snorm)); 

%*********************************
 % setting the size of the image basically within 1232, 1623
%***********************************
SizeofSnormal=size(Snorm);
FI = SizeofSnormal; 
%%
%***********************************************************************
% envelope of the image of the Snorm values could be used to identify hasField
% in other to identify the envelope of the image: It has been identified
% that hilbert and envelope command in matlab could be utilized to do that
% [up, lu] = envelope(Snorm);
% figure; imagesc(up); 
%Degree of coherence based on the different arms data and the enveloped
%mutualPic = (I1Cut.*I2Cut.*hilbert(up).^2)/(I1Cut.^2);
% %figure; imagesc(mutualPic);  
%coherence value
% disp('The degree of Coherence : '); 
% u = sum(I1Cut.*I2Cut.*up.^2)/sum(I1Cut.^2);
% disp(u); 
%*************************************************************************
%%

%*************************************************
% cropping section of the image using a circle 
% using a meshgrid command to form the circle to block sections of the fftimage
%*************************************************

[X, Y] = meshgrid(1 : SizeofSnormal(2), 1 : SizeofSnormal(1));

%********************************************
%first Radius for deleting the center image 
%*******************************************

Z =sqrt((X-SizeofSnormal(2)/2).^2+(Y-SizeofSnormal(1)/2).^2);

%***********************************
% Centering of the selected image
%***********************************
CenterSnorm = forSnorm;
if(vect==1)
 absX=abs(X-FI(2)/2);
 absY=abs(Y-FI(1)/2);
 CenterSnorm(absX<=1)=0;
 CenterSnorm(absY<=1)=0;
end

%***********************************
%setting zero frequency of the image to 0 
%****************************************

CenterSnorm(Z<kobbys_circle_radius/2)=0;

% Deleting the second section of the two frequencies 

[maxValue,indMaximum]=max(CenterSnorm(:));
[ny,nx]=ind2sub(SizeofSnormal, indMaximum);
Rshift=sqrt((X-nx).^2+(Y-ny).^2);
%%
%set values outside side maximum zero
Fcut=CenterSnorm;
Fcut(Rshift>kobbys_circle_radius)=0;

%circshift or shiftdim could be used for the shifting perpose
% shift circled area to the center
Fcenter=circshift(Fcut,[round(FI(1)/2 - ny)+ShifterY, round(FI(2)/2 - nx)+ShifterX]);
FFtfinal = ifftshift(Fcenter); 
% Actual image (Complete)
InverseFourier = 2.*(ifft2(FFtfinal));
 
selector = 546; 
InverseFourier_Collector =  InverseFourier(selector:selector+4, :);
AnotherSrijoyee = [(AnotherSrijoyee); abs(sum(InverseFourier_Collector)/5)];
AnotherSrijoyee2 = [(AnotherSrijoyee2); circshift(sort(abs((sum(InverseFourier_Collector)/5))), [0, 850])]; 
phaseofSpatial = angle(InverseFourier);


% %%
% % %Plotting of figures
% 
% 
% h=figure; 
% 
% %image of one arm opened. 
% subplot(4,2,1)
% imagesc(I1Cut)
% axis equal; axis tight;
% title('(a) $S_{1}(x,y)$', 'Interpreter', 'latex')
% colorbar;
% 
% %image of one arm opened 
% subplot(4,2,2)
% imagesc(I2Cut)
% axis equal; axis tight;
% title('(b) $S_{2}(x,y)$', 'Interpreter', 'latex')
% colorbar;
% 
% % image of both arms opened
% subplot(4,2,3)
% imagesc(ICut)
% axis equal; axis tight;
% title(' (d) $S(x,y)$', 'Interpreter','latex')
% colorbar;
% 
% %Subplot for imaging the normalized image(Snorm)
% subplot(4,2,4)
% hold on
% imagesc(Snorm)
% axis equal; axis tight;
% title(' (d) $S_{norm}(x,y)$', 'Interpreter','latex')
% colorbar;
% 
% %Subplot for imaging the FourierTransform of the Snorm
% subplot(4,2,5)
% hold on
% I = imagesc(abs(forSnorm));
% [nx,ny,d] = size(I) ;
% [X,Y] = meshgrid(1:ny,1:nx) ;
% hold on
% th = linspace(0,2*pi) ;
% xc = round(px)+r*cos(th); 
% yc = round(py)+r*sin(th); 
% plot(xc,yc,'r') ;
% axis equal; axis tight;
% title(' (e) $FT[S_{norm}(x,y)]$', 'Interpreter','latex')
% colorbar;
% 
% %subplot for imaging the FTcropped image 
% subplot(4,2,6)
% hold on
% imagesc(abs(Fcenter))
% axis equal; axis tight;
% title(' (f) $FT[S_{norm}(x,y) cropped]$', 'Interpreter','latex')
% colorbar;
% 
% %subplot for imaging the mcoherence 
% subplot(4,2,7)
% hold on
% imagesc(abs(InverseFourier))
% axis equal; axis tight;
% title(' (g) $|\mu(x,y)|$', 'Interpreter','latex')
% colorbar;
% 
% % subplot for imaging the argumentSnorm 
% subplot(4,2,8)
% hold on
% imagesc(abs(phaseofSpatial))
% axis equal; axis tight;
% title(' (h) $arg[\mu(x,y)]$', 'Interpreter','latex')
% colorbar;
% % %%
% %saving images for each loop 
% fname = '/Users/kobbyTilly/Desktop/Msc(UEF)_Biophotonics/Third_Period/LAB2/NewLab2/Spatial_coherence_IMS/Images';
% saveas(h,fullfile(fname, sprintf('FIG%d.png',k)));
% saveas(h,fullfile(fname, sprintf('FIG%d.eps',k)));
% %%
% % Section for the mouse Clicking
% 
% % disp('The first Analysis is done: Do you want to continue to the Next one: ');
% %         keydown = waitforbuttonpress;
% %         if (keydown == 0)
% %             disp('Mouse button was pressed');
% %         else
% %             disp('Key was pressed');
% %         end
% %         close(h);
% %
%  close all; 
end

 
f = msgbox('WFI Complete_Analysis: Developed By Ibrahim Issah','Success');
tilly  = AnotherSrijoyee; 
[staff, ahmed] = size(tilly); 
Kplus = zeros(staff, ahmed); 
for n = 1:length(matFiles) 
    Kplus(n, :) = [tilly(n,((staff+1)-n)*20+ShifterX:end), tilly(n, 1:((staff+1)-n)*20-1+ShifterY)];
end
finalimage =  [Kplus(:, 1001:end), Kplus(:, 1:1000)];
 mustafa = figure(100); 
 imagesc(finalimage)
title('$|\mu (\delta{x}, \bar{x})|$', 'Interpreter','latex')
xlabel('$ \delta{x} [arb.u.]$', 'Interpreter', 'latex'); 
ylabel('$ \bar{x} [arb.u.]$', 'Interpreter', 'latex'); 
colorbar; caxis([0 ,1]); 
saveas(mustafa,fullfile(fname, sprintf('FIG_fimage.png')));
saveas(mustafa, fullfile(fname, sprintf('FIG_fimage.eps')));
srijoyee = figure(99); 
imagesc(AnotherSrijoyee)
title('$|\mu (\delta{x}, \bar{x})|$', 'Interpreter','latex')
xlabel('$ \delta{x} [arb.u.]$', 'Interpreter', 'latex'); 
ylabel('$ \bar{x} [arb.u.]$', 'Interpreter', 'latex'); 
saveas(srijoyee,fullfile(fname, sprintf('FIG_tilted.png')));
saveas(srijoyee, fullfile(fname, sprintf('FIG_tilted.eps')));


%%
% *Done with this code* 