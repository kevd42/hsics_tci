function colmap = colorMapGen(w, n)
% Color conversion taken from codingmess.blogspot.com
% Table
R = 0.0;
G = 0.0;
B = 0.0;
if w >= 380 && w < 440
    R = -(w - 440.) / (440. - 350.);
    G = 0.0;
    B = 1.0;
end;
if w >= 440 && w < 490
    R = 0.0;
    G = (w - 440.) / (490. - 440.);
    B = 1.0;
end
if w >= 490 && w < 510
    R = 0.0;
    G = 1.0;
    B = -(w - 510.) / (510. - 490.);
end
if w >= 510 && w < 580
    R = (w - 510.) / (580. - 510.);
    G = 1.0;
    B = 0.0;
end
if w >= 580 && w < 645
    R = 1.0;
    G = -(w - 645.) / (645. - 580.);
    B = 0.0;
end
if w >= 645 && w <= 780
    R = 1.0;
    G = 0.0;
    B = 0.0;
end
% grey level colormap for INFRARED and ULTRAVIOLET 
if w>780 || w<380
    R = 1;
    G = 1;
    B = 1;
end

% Intensity correction
% SSS = 0.0
% if w >= 380 && w < 420
%     SSS = 0.3 + 0.7*(w - 350) / (420 - 350);
% end
% if w >= 420 && w <= 700
%     SSS = 1.0;
% end
% if w > 700 && w <= 780
%     SSS = 0.3 + 0.7*(780 - w) / (780 - 700);
% end
% SSS = SSS*255;

maxCol = [R G B];
maxColProp=0.5;
colmap = [linspace(0, maxCol(1), floor(n*maxColProp))', ...
          linspace(0, maxCol(2), floor(n*maxColProp))', ...
          linspace(0, maxCol(3), floor(n*maxColProp))'];
colmap = [colmap(1:end-1,:);
          linspace(maxCol(1), 1, ceil(n*(1-maxColProp))+1)', ...
          linspace(maxCol(2), 1, ceil(n*(1-maxColProp))+1)', ...
          linspace(maxCol(3), 1, ceil(n*(1-maxColProp))+1)'];