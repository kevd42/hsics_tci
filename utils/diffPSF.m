function H = diffPSF(ms,ns,wls,rad_max)
%% DIFFPSF
% Analytical Fraunhofer diffraction kernel of square aperture
%
% Usage:
% H = DIFFPSF(ms,ns,wls,rad_max)
%
% Input:
% ms      - resolution of FPA in a single direction
% ns      - resolution of scene in a single direction
% wls - wavelengths vector 
%
% Output:
% H      - 3-D array containing one kernel per wavelength
%
% see README.txt
%
% Authors: K. Degraux, V. Cambareri
% Date: 2017
%
% (c) UCLouvain 2018

if numel(ns)==1
    ns = [ns,ns];
end
if numel(ms)==1
    ms = [ms,ms];
end
s = ms + ns - 1;


% rad_max = max(wls)*F / (wfpa*wslm);

% Set the system parameters

wslm  = 80e-6; % width of slm pixel in m

Dlens = sqrt(2)*wslm*(max(s)); % min diameter of lens
F     = 0.706005*Dlens; % focal length in m (cannot be smaller than the diameter of the lense?

wfpa   = max(wls)*F / (rad_max*wslm);  %10*5.5e-6; % width of fpa pixel in m


% Lfpa   = ms * wfpa;
% px_bin = wfpa/5.5e-6;


% Dlens = sqrt(2)*wslm*(max(s)); % min diameter of lens



% fprintf('should be more than 0.7: %f\n',F/Dlens);


% Set empty output
nf = numel(wls);
H  = zeros([s,nf]);

% Plane wave input
% hi =  @(i, j, f, delta, l)  (sinc(delta*j/(l*f)).*sinc(delta*i/(l*f))).^2;
hi =  @(i, f, delta, l)  sinc(delta*i/(l*f)).^2;

% Sampling grid
step = 0.2;
i = (-(max(s)-1)/2:step:(max(s)-1)/2)*wfpa;
% j = (-(s(2)-1)/2:step:(s(2)-1)/2)*wfpa;
% [jj, ii] = meshgrid(j, i);
hwin = fspecial('disk', (max(s)-1)/2);

for band = 1:numel(wls)
	%% Compute monochromatic PSF
%     hh = hi(ii, jj, F, wslm, wls(band));
%     hint = diff(diff(cumtrapz(i, cumtrapz(j,hh,2)),1,1),1,2)/(step^2);
% 	H(:,:,band) = imresize(hint,s,'box'); % approximate the integral average
%     H(:,:,band) = H(:,:,band)/sum(sum(abs(H(:,:,band))));
    
    hh = hi(i, F, wslm, wls(band));
    hint = diff(cumtrapz(i, hh))/(step);
    % approximate the integral average
	HH = hwin.*imresize(hint(:)* hint(:).',[max(s),max(s)],'box');
    
    % In case s(1)~=s(2), crop H
    % Normalization of the kernel
    H(:,:,band) = HH(ceil(max(s)/2)+floor(-s(1)/2+1:s(1)/2),...
                     ceil(max(s)/2)+floor(-s(2)/2+1:s(2)/2))/sum(abs(HH(:)));
    
end

% Example:
% H = diffPSF(128, 512, linspace(400,700,16)*1e-9)
% for i = 1:16
%     mesh(H{i})
%     pause(1)
% end