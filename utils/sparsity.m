%SPARSITY sparsity of the signal
%   usage:
%   K = sparsity(x)
%   number of nonzero elements in x (equivalent to nnz(x) )
%   K = sparsity(x,epsilon)
%   number of elements in abs(x)/norm(x) that are > epsilon 
%   epsilon must be >= 0
%
%   Author : K. Degraux
%   Date : Oct 2013
%  (c) UCLouvain 2018

function K = sparsity(x,epsilon)



    if nargin == 1
        epsilon = 0;
    elseif isempty(epsilon)
        epsilon = 0;
    elseif epsilon < 0
        error('epsilon must be non negative');
    end
    K = nnz(abs(x)/norm(x)>epsilon);

end