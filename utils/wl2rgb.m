function M = wl2rgb(w)
% Color conversion taken from codingmess.blogspot.com
[R, G, B] = arrayfun(@wl2rgb_map, w);
M = [R; G; B]';
end

function [R, G, B] = wl2rgb_map(w)
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
end