function x = importCaveGroundtruth(gt_name,gt_center,n)
%% IMPORTCAVEGROUNDTRUTH(gt_name,n)
% Imports the groundtruth from the Cave multispectral dataset [1].
% Returns a 3-D array with the specified dimensions n.
% List of gt_name accepted:
% balloons, beads, cd, chart_and_stuffed_toy, clay, cloth, egyptian_statue,
% face, fake_and_real_{beers,food,lemon_slices,lemons,peppers,strawberries,
% sushi,tomatoes}, feathers, flowers, glass_tiles, hairs, jelly_beans, 
% oil_painting, paints, photo_and_face, pompoms, real_and_fake_{peppers,
% apples}, sponges, stuffed_toys, superballs, thread_spools, watercolors.
%
% [1] "Generalized Assorted Pixel Camera: Post-Capture Control of 
%      Resolution, Dynamic Range and Spectrum,"
%      F. Yasuma, T. Mitsunaga, D. Iso and S.K. Nayar,
%      IEEE Transactions on Image Processing,
%      Vol. 99, Mar. 2010.
%
% Author: K. Degraux
% Date:   7 Apr 2017
%
%  (c) UCLouvain 2018

if nargin == 1
    n = [512,512,16];
end

folder = ['images/complete_ms_data/',gt_name,'_ms/',gt_name,'_ms/'];

X = zeros(512,512,16);
for l = 1:16
    X(:,:,l) = imread([folder,gt_name,'_ms_',sprintf('%2.2d',l+7),'.png']);
end

% Record extreme values of raw data
maxval = max(X(:));
minval = min(X(:));
    
% Resize in spatial domain (bicubic interpolation)
if isempty(gt_center)

    X = imresize(X,n(1:2));
    % Clamp the resized image to initial range
    X = min(max(X,minval),maxval);
else
    X = X(gt_center(1)+(-n(1)/2:n(1)/2-1),gt_center(2)+(-n(2)/2:n(2)/2-1),:);
end

% Resize in spectral domain (cubic interpolation) % Already 16 bands :-)
%X = permute(imresize(permute(X,[3,1,2]),[n(3),n(1)]),[2,3,1]);

% Bring back the range between 0 and 1
X = (X-minval)/(maxval-minval);

x=X(:);

end