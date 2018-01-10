function f = HSinspector(Vol,f,clim,cmode,wlim)
% HSINSPECTOR create an interactive visualization tool for 3-D arrays (e.g.
% hyperspectral data cubes). 
% Usage:
%   HSinspector(Vol)
%   HSinspector(Vol,f)
%   HSinspector(Vol,f,clim)
%   HSinspector(Vol,f,clim,cmode)
%   HSinspector(Vol,f,clim,cmode,wlim)
%   f = HSinspector(...)
%   All input parameters are optional. An empty array will select the
%   default value.
%   
%   Vol   : 3-D array to display.
%   f     : int to select the figure number in which we create the tool.
%   clim  : [cmin cmax] minimum and maximum value of the color space
%   cmode : 'hs' (default) or a colormap. choosing 'hs' will set a different
%           color for every band, approx. corresponding to reality.
%   wlim  : sets the min and max value of the wavelength (in nm) so that if 
%           cmode is 'hs', the colors correspond (approx.) to reality.
%
% Author : K. Degraux, ISPGroup (UCLouvain)
% 2014
%
%  (c) UCLouvain 2018


if nargin<1
    Vol = randn(512,512,16);
end

if nargin<2||isempty(f)
    f = figure;
    newfig = true;
else
    if ishandle(f)
        set(0,'CurrentFigure',f);
        clf(f);
        newfig = false;
    else
        figure(f);
        newfig = true;
    end
end

if nargin<3 || isempty(clim)
    minVol = min(Vol(Vol>-Inf));
    maxVol = max(Vol(Vol<Inf));
    clim   = [minVol,maxVol];
else
    minVol = clim(1);
    maxVol = clim(2);
    Vol(Vol<minVol) = minVol;
    Vol(Vol>maxVol) = maxVol;
end

if nargin<4 || isempty(cmode)
    cmode = 'hs';
end

n = size(Vol);

if nargin<5 || isempty(wlim)
    wlim = [473,632];
end

if ~strcmp(cmode,'hs')
    cmap1 = colormap(cmode);
    cmaplength = length(cmap1);
else
    cmaplength = 256;
end

cmap2 = hot(cmaplength);


if numel(n)<3
    imagesc(Vol,clim);
    axis equal tight;
    return;
end


set(0,'Units','pixels');
ScreenSize = get(0,'ScreenSize');

Spacer = 0.05;

RelFigW = (n(2) + max(n(3),(n(1)/2))) * (1+3*Spacer);
RelFigH = (n(1) + max(n(3),(n(1)/2))) * (1+3*Spacer);


FigRatio = RelFigW/RelFigH;
ScreenRatio = ScreenSize(3)/ScreenSize(4);

if FigRatio>ScreenRatio
    FigSize(3) = ScreenSize(3);
    FigSize(4) = FigSize(3) / FigRatio;
else
    FigSize(4) = ScreenSize(4);
    FigSize(3) = FigSize(4) * FigRatio;
end

FigSize(1) = round(ScreenSize(3)/2 - FigSize(3)/2);
FigSize(2) = round(ScreenSize(4)/2 - FigSize(4)/2);



%  Create the GUI
set(f,'WindowButtonUpFcn',@figure_WindowButtonUpFcn);
if newfig
    set(f,'Position',FigSize);
end




cclim = [minVol,2*maxVol-minVol+2/cmaplength];




curr_i = round(n(1)/2);
curr_j = round(n(2)/2);
curr_l = round(n(3)/2);
prev_i = curr_i;
prev_j = curr_j;
prev_l = curr_l;
hscmode = strcmp(cmode,'hs');
if hscmode
    w = linspace(min(wlim),max(wlim),n(3));
    cmap1 = colorMapGen(w(curr_l), cmaplength);
end
cmap = [cmap1;cmap1(end,:);cmap2(1,:);cmap2];
colormap(cmap);

im_ij = cell(n(3),1);
im_lj = cell(n(1),1);
im_il = cell(n(2),1);
ha1 = axes('Units','Normalized',...
    'Position',[Spacer,Spacer,n(2)/RelFigW,n(1)/RelFigH],...
    'ButtonDownFcn',@axes1_ButtonDownFcn);
hold on;
im_ij{curr_l} = imagesc(Vol(:,:,curr_l),cclim);
set(im_ij{curr_l},'hittest','off');
axis tight;
axis ij;
set(ha1,'YaxisLocation','right','XaxisLocation','top');
hcmap = colorbar('location','west');
set(hcmap,'Units','Normalized','position',[0.5*Spacer,Spacer,Spacer/5,n(1)/RelFigH],'Ylim',[minVol,maxVol]);


ha2 = axes('Units','Normalized',...
    'Position',[Spacer,2*Spacer+n(1)/RelFigH,n(2)/RelFigW,max(n(3),n(1)/2)/RelFigH],...
    'ButtonDownFcn',@axes2_ButtonDownFcn);

hold on;
im_lj{curr_i} = imagesc(permute(Vol(curr_i,:,:)+maxVol+2/64-minVol,[3,2,1]),cclim);
set(im_lj{curr_i},'hittest','off');
axis tight;
axis xy;

ha3 = axes('Units','Normalized',...
    'Position',[2*Spacer+n(2)/RelFigW,Spacer,max(n(3),n(1)/2)/RelFigW,n(1)/RelFigH],...
    'ButtonDownFcn',@axes3_ButtonDownFcn);

hold on;
im_il{curr_j} = imagesc(permute(Vol(:,curr_j,:)+maxVol+2/64-minVol,[1,3,2]),cclim);
set(im_il{curr_j},'hittest','off');
axis tight;
axis ij

ha4 = axes('Units','Normalized','Position',[2*Spacer+n(2)/RelFigW,2*Spacer+n(1)/RelFigH,max(n(3),n(1)/2)/RelFigW,max(n(3),n(1)/2)/RelFigH]); 
    
hp1 = cell(1,4);
ht1 = cell(1,2);
hp2 = cell(1,2);
ht2 = cell(1,1);
hp3 = cell(1,2);
ht3 = cell(1,1);
click_state = 0;

axes1_update
axes2_update
axes3_update
axes4_update

    function axes1_update
        set(f,'CurrentAxes',ha1);
        axesLimits = [get(ha1,'XLim').' get(ha1,'YLim').'];
        if isempty(im_ij{curr_l})
            hold on;
            im_ij{curr_l} = imagesc(Vol(:,:,curr_l),cclim);
            set(im_ij{curr_l},'hittest','off');
        else
            set(im_ij{curr_l},'Visible','on');
        end
        xlim(axesLimits(:,1));
        ylim(axesLimits(:,2));

        
        if prev_l~=curr_l
            set(im_ij{prev_l},'Visible','off');
            if hscmode
                cmap1 = colorMapGen(w(curr_l), cmaplength);
                cmap = [cmap1;cmap1(end,:);cmap2(1,:);cmap2];
                colormap(cmap);
            end
            set(hcmap,'ylim',[minVol,maxVol]);
        end
        
        if ~isempty(hp1{1})
            for i = 1:4
                delete(hp1{i});
            end
        end
        hp1{1} = plot([curr_j,curr_j],[0,n(1)]+0.5,'r:');
        hp1{2} = plot([0,n(2)]+0.5,[curr_i,curr_i],'r:');
        hp1{3} = plot(curr_j,curr_i,'or');
        hp1{4} = plot(curr_j,curr_i,'+r');
        for i = 1:4
            set(hp1{i},'hittest','off');
        end
        if ~isempty(ht1{1})
            for i = 1:2
                delete(ht1{i});
            end
        end
        ht1{1} = text(axesLimits(2,1)+0.055*diff(axesLimits(:,1)),curr_i,num2str(curr_i),'horizontalalignment','center');
        ht1{2} = text(curr_j,axesLimits(1,2)-0.05*diff(axesLimits(:,2)),num2str(curr_j),'horizontalalignment','center');
    end
    function axes2_update
        set(f,'CurrentAxes',ha2);
        axesLimits = [get(ha2,'XLim').' get(ha2,'YLim').'];
        if isempty(im_lj{curr_i})
            hold on;
            im_lj{curr_i} = imagesc(permute(Vol(curr_i,:,:)+maxVol+2/64-minVol,[3,2,1]),cclim);
            set(im_lj{curr_i},'hittest','off');
        else
            set(im_lj{curr_i},'Visible','on');
        end
        xlim(axesLimits(:,1));
        ylim(axesLimits(:,2));
        set(ha2,'xtick',[]);
        
        if prev_i~=curr_i
            set(im_lj{prev_i},'Visible','off');
        end
        
        
        if ~isempty(hp2{1})
            for i = 1:2
                delete(hp2{i});
            end
        end
        hp2{1} = plot([0,n(2)]+0.5,[curr_l,curr_l],'g.:','markersize',15);
        hp2{2} = plot([curr_j,curr_j],[0,n(3)]+0.5,'b.:','markersize',15);
        for i = 1:2
            set(hp2{i},'hittest','off');
        end
        if ~isempty(ht2{1})
            delete(ht2{1});
        end
        ht2{1} = text(axesLimits(1,1)-0.06*diff(axesLimits(:,1)),curr_l,num2str(curr_l),'horizontalalignment','center');
    end
    function axes3_update
        set(f,'CurrentAxes',ha3);
        axesLimits = [get(ha3,'XLim').' get(ha3,'YLim').'];
        if isempty(im_il{curr_j})
            hold on;
            im_il{curr_j} = imagesc(permute(Vol(:,curr_j,:)+maxVol+2/64-minVol,[1,3,2]),cclim);
            set(im_il{curr_j},'hittest','off');
        else
            set(im_il{curr_j},'Visible','on');
        end
        xlim(axesLimits(:,1));
        ylim(axesLimits(:,2));
        set(ha3,'ytick',[]);
        hold on;
        
        if prev_j~=curr_j
            set(im_il{prev_j},'Visible','off');
        end
        
        if ~isempty(hp3{1})
            for i = 1:2
                delete(hp3{i});
            end
        end
        hp3{1} = plot([curr_l,curr_l],[0,n(1)]+0.5,'g.:','markersize',15);
        hp3{2} = plot([0,n(3)]+0.5,[curr_i,curr_i],'b.:','markersize',15);
        for i = 1:2
            set(hp3{i},'hittest','off');
        end
        if ~isempty(ht3{1})
            delete(ht3{1});
        end
        ht3{1} = text(curr_l,axesLimits(2,2)+0.05*diff(axesLimits(:,2)),num2str(curr_l),'horizontalalignment','center');
    end
    function axes4_update
        set(f,'CurrentAxes',ha4);
        plot(permute(Vol(curr_i,curr_j,:),[3,1,2]),'linewidth',2);hold on;
        plot(curr_l,Vol(curr_i,curr_j,curr_l),'go','linewidth',2,'markersize',10); hold off;
        xlim([1,n(3)]);
        if maxVol>minVol
            ylim([minVol,maxVol]);
        end
        drawnow;
    end
   
    function axes1_ButtonDownFcn(src,eventdata) %#ok<INUSD>
        
        click_state = 1;
        while click_state

            click_ij=round(get_coords(ha1));

            prev_i = curr_i;
            prev_j = curr_j;

            axesLimits = ceil(get(ha1,'YLim'));
            axesLimits(2) = axesLimits(2)-1;
            curr_i = max(min(axesLimits(2)-click_ij(1,2)+axesLimits(1),axesLimits(2)),1);
            curr_j = max(min(click_ij(1,1),n(2)),1);
            axes1_update;
            axes2_update;
            axes3_update;
            axes4_update;
            drawnow;
            pause(1/24);
        end

    end
    function axes2_ButtonDownFcn(src,eventdata) %#ok<INUSD>
        
        click_state = 1;
        while click_state

            click_ij=round(get_coords(ha2));

            prev_l = curr_l;
            prev_j = curr_j;

            curr_l = max(min(click_ij(1,2),n(3)),1);
            curr_j = max(min(click_ij(1,1),n(2)),1);
            axes1_update;
            axes2_update;
            axes3_update;
            axes4_update;
            drawnow;
            pause(1/24);
        end

    end
    function axes3_ButtonDownFcn(src,eventdata) %#ok<INUSD>
        
        click_state = 1;
        while click_state

            click_ij=round(get_coords(ha3));

            prev_i = curr_i;
            prev_l = curr_l;
            
            axesLimits = ceil(get(ha1,'YLim'));
            axesLimits(2) = axesLimits(2)-1;
            
            curr_i = max(min(axesLimits(2)-click_ij(1,2)+axesLimits(1),axesLimits(2)),1);
            curr_l = max(min(click_ij(1,1),n(3)),1);
            axes1_update;
            axes2_update;
            axes3_update;
            axes4_update;
            drawnow;
            pause(1/24);
        end

    end

    function figure_ResizeFcn(src,eventdata) %#ok<DEFNU,INUSD>
    end
    
    function figure_WindowButtonUpFcn(src,eventdata) %#ok<INUSD>
        click_state = 0;
    end


% % Make the GUI visible.
% set(f,'Visible','on');
   
   
end

function value = get_in_units(hObject,propName,unitType)

  oldUnits = get(hObject,'Units');  %# Get the current units for hObject
  set(hObject,'Units',unitType);    %# Set the units to unitType
  value = get(hObject,propName);    %# Get the propName property of hObject
  set(hObject,'Units',oldUnits);    %# Restore the previous units

end
function coords = get_coords(hAxes)

  %# Get the screen coordinates:
  coords = get_in_units(0,'PointerLocation','pixels');

  %# Get the figure position, axes position, and axes limits:
  hFigure = get(hAxes,'Parent');
  figurePos = get_in_units(hFigure,'Position','pixels');
  axesPos = get_in_units(hAxes,'Position','pixels');
  axesLimits = [get(hAxes,'XLim').' get(hAxes,'YLim').'];

  %# Compute an offset and scaling for coords:
  offset = figurePos(1:2)+axesPos(1:2);
  axesScale = diff(axesLimits)./axesPos(3:4);

  %# Apply the offsets and scaling:
  coords = (coords-offset).*axesScale+axesLimits(1,:);

end