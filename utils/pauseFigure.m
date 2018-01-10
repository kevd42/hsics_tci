function hFig = pauseFigure(hFig,init_str)
%% PAUSEFIGURE Create a figure window that allows to pause execution and 
%  enter debug mode
%
% Placed before a loop,
% hFig = PAUSEFIGURE; 
% hFig = PAUSEFIGURE('init');
% hFig = PAUSEFIGURE(hFig,'init');
% creates a figure (with an optionally specified number hFig)
% and returns it's handle. The figure will be used to pause the execution.
% 
% Inside the loop, 
% hFig = PAUSEFIGURE(hFig);
% Checks if the figure still exists. If it has been closed by the user, the
% iteration pauses and MATLAB enters in DEBUG MODE.
% 
% see also: keyboard, return, dbquit
%
% Author : K. Degraux
% Date   : Aug 2014
%  (c) UCLouvain 2018
init_bool = 0;

% initialization, hFig number is specified
if nargin == 2
    if strcmp(init_str,'init')
        init_bool = 1;
        if ishandle(hFig)
            close(hFig);
        end
    end
end
    
% initialization, hFig number is not specified
if nargin<1 || strcmp(hFig,'init')
    init_bool = 1;
    hFig = figure;
    close(hFig);
end

if ~ishandle(hFig)
    hFig = figure(hFig); 
    set(hFig,'name','Close this window to pause');
    pos = get(hFig,'OuterPosition');
    pos(3:4) = pos(3:4)./[2,4];
    set(hFig,'OuterPosition',pos);
    imagesc([]);
    axis off
    text(0.5,0.5,'\fontsize{20} \bf Close this window to pause','HorizontalAlignment','Center');
    if ~init_bool
        keyboard;
    end
end
drawnow;
end