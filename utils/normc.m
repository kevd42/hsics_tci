function nX = normc(X)
% NORMC
% Normalize the columns of X
% Usage: nX = normc(X)
% Author: K. Degraux
%  (c) UCLouvain 2018
norms = sqrt(sum(X.^2, 1));
nX = X*diag(1./norms);
end