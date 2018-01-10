
function [y,news] = proxNucl(x,lambda,i,s,n)
% [y,news] = PROXNUCL(x,lambda,i,s,n)
% Proximal operator of the nuclear norm with memory for accelerated 
% computations.
% Author : K. Degraux
% Date : Oct 2013
%
%  (c) UCLouvain 2018

X = tensor(x,n);
Xi = tenmat(X,i);

[newXi,news] = stsvds(Xi.data,lambda,s);

Y = tensor(newXi,size(X));
y = Y(:);

end
