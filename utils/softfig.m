function h = softfig(h)
% SOFTFIG
% Create figure or set current figure without poping it up.
% Author : K. Degraux
%  (c) UCLouvain 2018
if nargin<1
    h = figure;
else
    if ishandle(h);
        set(0,'CurrentFigure',h);
    else
        figure(h);
    end
end
end