
function a = proxIC(x,Phi,L,U,y)
% a = PROXIC(x,Phi,L,U,y)
% Author : K. Degraux
% Date : Oct 2013
%
%  (c) UCLouvain 2018

if condest(L)>1e6 ||condest(U)>1e6
    L = [];
    U = [];
end

if (isempty(L)||isempty(U))
    a = x + Phi'*((Phi*Phi')\(y-Phi*x));
else
    a = x + Phi'*(U\(L\(y-Phi*x)));
end

end
