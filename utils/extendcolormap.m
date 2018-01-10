

function newmap = extendcolormap(map,ncol)
% EXTENDCOLORMAP
% newmap = EXTENDCOLORMAP(map,ncol)
% Each color of map is extended from black to its value.
%
% Example:
% imagesc(randn(10))
% colormap(extendcolormap(jet(8),16))
% colorbar
%
% Author : K. Degraux
% Date : Oct 2013
%  (c) UCLouvain 2018

for i = 1:size(map,1)
    newmap((i-1)*ncol+1:ncol*i,:) = linspace(0,1,ncol)'*map(i,:);
end

end
