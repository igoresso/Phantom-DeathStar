% Spatial resolution measurement
%
% Copyright (c) 2021
% Youheng (Peter) Zeng
% Melbourne Brain Centre Imaging Unit (MBCIU)
% University of Melbourne

clear variables
close all
clc

% Read DICOM
[fileName, filePath]    = uigetfile('*.dcm', 'Select a DICOM');
info                    = dicominfo(fullfile(filePath,fileName));
Y                       = dicomread(info);
pixel_spacing           = info.PixelSpacing(1);

figure
imshow(Y,[])
title('Imported Image')
camzoom(5)
[cx,cy,c]  = improfile(50);
% [c_im,rec] = imcrop;
[pks,locs] = findpeaks(c);

figure
x = info.PixelSpacing(1)*(1:length(c));
plot(x,c);
hold on
x1 = info.PixelSpacing(1)*locs;
plot(x1,pks,'o');
xlabel('Distance, mm');
ylabel('Intensity, a.u.');
title('Line Profile');

figure
imshow(Y,[])
h = images.roi.Line(gca,'Position',[cx(1) cy(1);cx(end) cy(end)],'Label',...
'Inspection line','LabelAlpha',0,'LineWidth',1,'MarkerSize',3,'LabelTextColor','red');
camzoom(5);