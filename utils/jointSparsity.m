%JOINTSPARSITY row or column joint sparsity of the signal with tolerance
%
%   usage:
%   K = JOINTSPARSITY(X)
%   number of nonzero rows in the matrix X
%   K = JOINTSPARSITY(X,dim)
%   number of nonzeros rows (dim == 1, default) or columns (dim == 2) in
%   matrix X
%   K = JOINTSPARSITY(X,dim,epsilon)
%   number of rows (dim == 1, defaults) or columns (dim == 2) in X such
%   that any normalized element is >epsilon in absolute value
%
%   epsilon must be >= 0
%
%   Author : K. Degraux
%   Date : Oct 2013
%  (c) UCLouvain 2018
function K = jointSparsity(X,dim,epsilon)

if nargin == 1
    dim = 1;
    epsilon = 0;
end
if nargin == 2
    epsilon = 0;
end
if isempty(dim)
    dim = 1;
elseif dim>2 || dim <1
    dim = 1;
end

if isempty(epsilon)
    epsilon = 0;
elseif epsilon < 0;
    error('epsilon must be non negative');
end
% normX = sqrt(sum(X.^2,dim));
% if dim == 1
%     normX = repmat(normX,[size(X,1),1]);
% else
%     normX = repmat(normX,[1,size(X,2)]);
% end
% K = nnz(sum(abs(X)./normX>epsilon,3-dim));

normX = sqrt(sum(X.^2,3-dim));

K = nnz(abs(normX)./norm(normX)>epsilon);


end