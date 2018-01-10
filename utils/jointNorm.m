
function out = jointNorm(X,p1,p2)
% JOINTNORM(X,p1,p2)
% X is a matrix
% p1 is the exponent of the norm computed for each column (default 2)
% p2 is the exponent of the norm computed for the p1-norms vector (default
% 1)
%
% Author : K. Degraux
% Date : April 2012
% see also : jointSparsity, jointHs, norm
%  (c) UCLouvain 2018

    if nargin < 2
        p1=2;
        p2=1;
    end
    out = norm(sum(abs(X).^p1,1).^(1/p1),p2);
    
    
end
