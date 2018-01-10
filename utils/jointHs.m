

function [X,suppX] = jointHs(A,s)
% JOINTHS(A,s)
% A and X are matrices
% s is the number of non zero columns in X
% non linear operator taht sets all but the largest (in L2 norm) s
% columns of A to zero. If there is no unique such set, a set can be
% selected either randomly or based on a predefined ordering of the 
% elements. 
% Author : K. Degraux
% Date : April 2012
% see also : Hs, jointSparsity, jointNorm
%  (c) UCLouvain 2018

    X = A;
    
    [~,ix] = sort(sum(A.^2,1),'descend');
    
    X(:,ix(s+1:end)) = 0;
    if (nargout == 2)
        suppX = ix(1:s);
    end
end
