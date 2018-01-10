function y0 = interpolateFPA(y,ms,nl,idL)
% INTERPOLATEFPA - Use scattered linear interpolation to infer the data
% that was not acquired at the detector.
%
% Usage:
%  y0 = INTERPOLATEFPA(y,ms,nl,idL)
%
% Inputs:
%  y   is the vector of measurements
%  ms  = [mi,mj] is the spatial size of the focal plane array
%  nl  is the number of bands
%  idL (of size ms) is the index map of the spectral filters.
%
% Output:
%  y0  is the vectorized array of size mi x mj x nl containing the
% interpolation.
%
% See also: scatteredInterpolant
% 
% (c) K. Degraux, 2017
%  (c) UCLouvain 2018

mi = ms(1);
mj = ms(2);
[J,I]   = meshgrid(1:mj,1:mi);

method = 'linear'; 
bandspacing = 1; % Arbitrary virtual spacing (in pixels units) between 
                   % bands in the cube
interpolant = scatteredInterpolant(J(:),I(:),idL(:)*bandspacing,y(:),method);
y0 = interpolant(repmat(J,1,1,nl),repmat(I,1,1,nl),...
                 repmat(permute(1:nl,[1,3,2]),mi,mj,1) );
             

y0 = y0(:);

end