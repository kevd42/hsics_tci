
function v = proxBl2(x,y,epsilon)
% v = PROXBL2(x,y,epsilon)
% prox of the indicator function of (ie projection on) the convex set
%
% norm(x-y) <= epsilon
%
% Author : K. Degraux
% Date : Oct 2013
%  (c) UCLouvain 2018

if epsilon > 0
    nxy = norm(y-x);
    if nxy>epsilon
        v = x + (y-x)*max(1-epsilon/nxy,0);
    else
        v = x;
    end
else
    v = y;
end

end
