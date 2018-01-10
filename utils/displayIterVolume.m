function displayIterVolume(Psi0,u,n,h)
% DISPLAYITERVOLUME
% Helper for in-loop display of the current volume.
% Usage DISPLAYITERVOLUME(x,n,h)
% Author: K. Degraux
% (c) UCLouvain 2018
Xten = reshape(Psi0*u,n);
if n(3)>1
    HSinspector(Xten,h);
%     HSinspector(Xten,h);
else
    softfig(h);
    imagesc(Xten);
    colormap gray;
end
set(h,'name','Current iterate');

% Yten = reshape(Phi*u,m);
% if m(3)>1
%     HSinspector(Yten,h*10,[],'jet');
% else
%     softfig(h*10);
%     imagesc(Yten);
%     colormap jet;
%     colorout;
% end
% set(h*10,'name','Reobserved iterate');