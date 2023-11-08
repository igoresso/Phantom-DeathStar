% Slice profile (double wedge) measurement with interpolation
%
% Copyright (c) 2023
% Youheng (Peter) Zeng
% Melbourne Brain Centre Imaging Unit (MBCIU)
% University of Melbourne

clear variables
close all
clc

% Parameters
slope   = 15;   % Wedge angle
w       = 10;   % Shift aug data
d       = 3;    % Central difference factor
f       = 16;   % Interpolation factor

% Read DICOM
[fileName, filePath]    = uigetfile('*.dcm', 'Select a DICOM');
info                    = dicominfo(fullfile(filePath,fileName));
Y                       = dicomread(info);
pixel_spacing           = info.PixelSpacing(1);

figure
imshow(Y,[])
title('Imported Image')

% Extract the wedge
[c_im_1, rec_1] = imcrop;
rec_2           = [ rec_1(1), ...
                    rec_1(2) - (rec_1(4) + w), ...
                    rec_1(3), ...
                    rec_1(4) ]; % Create rec_2
c_im_2          = imcrop(Y,rec_2);

% QA step
figure
plot(c_im_1')
hold on
plot(c_im_2')

% Average across y direction to get Edge Response function ERF - ROI1
ERF_1 = average_x(c_im_1);

% Differentiate the ERF to get the Projected Slice Profile
l_1     = length(ERF_1);
l_a_1   = (1:l_1) * pixel_spacing * tand(slope);
dy_1    = central_dist(ERF_1,d);
dy_1    = abs(dy_1);
x_1     = 1:f*l_1;
x_1     = pixel_spacing / f * tand(slope) * x_1; % Normalise X-axis
dy_1    = interp1(l_a_1,dy_1,x_1);

% Find two points to get FWHM
[y_max_1,I_1]   = max(dy_1);
n_dy_1          = dy_1/y_max_1;
Y_1             = abs(n_dy_1 - 0.5);
[~,lp_1]        = min(Y_1(1:I_1));
[~,rp_1]        = min(Y_1(I_1:end));
rp_1            = rp_1 + I_1 - 1 ;

% Average across the Y direction to get Edge Response function ERF - ROI2
ERF_2   = average_x(c_im_2);

% Differentiate the ERF to get the Projected Slice Profile
l_2     = length(ERF_2);
l_a_2   = (1:l_2) * pixel_spacing * tand(slope);
dy_2    = central_dist(ERF_2,d);
dy_2    = abs(dy_2);
x_2     = 1:f*l_2;
x_2     = pixel_spacing / f * tand(slope) * x_2; % Normalise X-axis
dy_2    = interp1 (l_a_2,dy_2,x_2);

% Find two points to get FWHM
[y_max_2,I_2]   = max(dy_2);
n_dy_2          = dy_2/y_max_2;
Y_2             = abs(n_dy_2 - 0.5);
[~,lp_2]        = min(Y_2(1:I_2));
[~,rp_2]        = min(Y_2(I_2:end));
rp_2            = rp_2 + I_2 - 1;

% Plot(Ensure the position of ROI is stored in rec_1)
figure
imshow(Y_1,[])
title('Phantom Image')
r_1 = rectangle('Position',rec_1);
r_1.EdgeColor = 'b';
r_1.LineWidth = 1;
hold on
r_2 = rectangle('Position',rec1);
r_2.EdgeColor = 'r';
r_2.LineWidth = 1;
hold off

figure
subplot(1,2,1);
plot(l_a_1,ERF_1);
hold on
plot(l_a_2,erf_2);
title('Edge Response Function');
xlabel('X, mm');
ylabel('Signal, a.u.');

subplot(1,2,2);
plot(x_1,n_dy_1);
hold on
plot(x_2,n_dy_2);
yline(0.5,'r');
hold on
plot((lp_1) * pixel_spacing/f * tand(slope),n_dy_1(lp_1),'b*');
plot((rp_1) * pixel_spacing/f * tand(slope),n_dy_1(rp_1),'b*');
plot((lp_2) * pixel_spacing/f * tand(slope),n_dy_2(lp_2),'r*');
plot((rp_2) * pixel_spacing/f * tand(slope),n_dy_2(rp_2),'r*');
title('Slice Profile');
xlabel('Z, mm');
ylabel('Normalised Signal');

% Standardize lp and lr
lp_1    = lp_1 * pixel_spacing/f ;
rp_1    = rp_1 * pixel_spacing/f ;
lp_2    = lp_2 * pixel_spacing/f ;
rp_2    = rp_2 * pixel_spacing/f ;
w_1     = rp_1 - lp_1;
w_2     = rp_2 - lp_2;

% Rotation correction
theta   = asin(abs(w_1-w_2) * sin(2*slope*2*pi/360) / (w_1 + w_2))/2;
theta_d = rad2deg(theta);
FWHM    = (w_1 * tand(slope) + w_2 * tand(slope)) / 2;

if w_1 > w_2
    c_FWHM = w_2 * tand(slope + theta_d);
else
    c_FWHM = w_1 * tand(slope + theta_d);
end

function dy = central_dist(y,d)
%   [B1_filtered, r] = filterMaps(B1_raw,x,y,z)
%       calculated central difference.
%    
% Inputs:
%   y       : sample points, 
%   d       : distance used in the calculation of differential values 
%             (2d is the distance between two points)
%
% Output:
%   dy      : central difference

n               = length(y);
dy              = zeros(n,1);
pixel_spacing   = 0.2604;
% if d==3
% A = [-1 9 -45 45 -9 1];%parameter matrix position:[-3:3]
% for i = (d+1) : (n-d)
%
% dy(i) = (A(1)*y(i-3)+A(2)*y(i-2)+A(3)*y(i-1) +A(4)*y(i+1)+...
% A(5)*y(i+2)+A(6)*y(i+3))/(60*0.2604);
%
% end
% end
for i = (d + 1) : (n - d)
    y_l = sum(y(i-d:i-1))/d/pixel_spacing;
    y_r = sum(y(i+1:i+d))/d/pixel_spacing;
    dy(i) = (y_r-y_l)/d+1/pixel_spacing;
end

% Normalisation
dy = dy./max(dy);

end

% Function of average
% input is the image of interest im, output is ERF :a single line profile across
% x direction
function ERF = average_x(im)
%   ERF = average_x(im)
%       computes the averaged edge response function
%    
% Input:
%   im      : image of interest
%
% Output:
%   ERF     : a single line profile across the X direction

[numRows, numCols] = size(im);
if numRows < numCols
    im = im.';
    [numRows, numCols] = size(im);
end

ERF = zeros(numRows,1);

for cols = 1:numRows
    r_d = im(cols,:);
    ERF(cols) = mean(r_d);
end

end