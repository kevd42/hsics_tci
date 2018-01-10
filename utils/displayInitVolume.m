function displayInitVolume(x,n,h)
% DISPLAYINITVOLUME
% Helper for in-loop display of the initial volume.
% Usage DISPLAYINITVOLUME(x,n,h)
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
set(h,'name','Initial point');

% Yten = reshape(y,m);
% if m(3)>1
%     HSinspector(Yten,h*10,[],'jet');
% else
%     softfig(h*10);
%     imagesc(Yten);
%     colormap jet;
%     colorout;
% end
% set(h*10,'name','Measurement Vector');