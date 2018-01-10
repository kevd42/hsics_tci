function displayTargetVolume(x,n,h)
% DISPLAYTARGETVOLUME
% Helper for in-loop display of the target volume.
% Usage DISPLAYTARGETVOLUME(x,n,h)
% Author: K. Degraux
% (c) UCLouvain 2018
Xten = reshape(x,n);
if n(3)>1
    HSinspector(Xten,h,[0 1]);
else
    softfig(h);
    imagesc(Xten,[0 1]);
    colormap gray;
end
set(h,'name','Target volume');