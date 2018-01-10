

function y = importMosaicExample(gt_name,roi,idL,m)
% IMPORTMOSAICEXAMPLE
% Import the experimental data from the mosaic sensor.
% Usage:
% y = IMPORTMOSAICEXAMPLE(gt_name,roi,idL,m)
% Inputs:
% gt_name: 'leaves' or 'RedYellow'
% roi: region of interest
% idL: mapping of the FPA to band indices
% m: size of the FPA
% Author: K. Degraux
%  (c) UCLouvain 2018
if isempty(roi)
    roi = [1,256,145,256];
end

switch gt_name
    case 'leaves';    cubename='images/Large/MosaicExample/VIS_leaves_FD1_FN0_original.hdr';
    case 'RedYellow'; cubename='images/Large/MosaicExample/RedYellowSet_3_FD1_FN0.hdr';
end
[ cubename,lines,samples,bands,defaultBands, interleave, wavelengths ] = parsecube( cubename ); %#ok<NASGU,ASGLU>

% disp([ num2str(lines) 'x' num2str(samples) 'x' num2str(bands) ' cube'])
cube = multibandread(cubename,[lines,samples,bands],'uint16',0,interleave,'ieee-le');
cube = cube(roi(1):roi(1)+roi(2)-1,roi(3):roi(3)+roi(4)-1,:);
cube = cube/max(cube(:));
Y = zeros(m(1:2));
for band = 1:bands
    Y(idL==band) = vec(cube(:,:,band));
end
y =Y(:); 