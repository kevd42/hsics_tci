
function colorout(cmin,cmax)
% COLOROUT
% Thin colorbar outsize the axes.
% Usage:
% COLOROUT(cmin,cmax)
% Author: K. Degraux
%  (c) UCLouvain 2018

subplot('position',[0.98 0.1 0.01 0.8])
imagesc([],[cmin,cmax])
axis off
colorbar('location','WestOutside');

end