% Author : K. Degraux
% Date : Oct 2013

function V = proxl12(X,gamma)
% V = PROXL12(X,gamma)
% Prox of the matrix function ||X||12, i.e.
%
% sum_i norm(X(:,i))
% Author : K. Degraux
% Date : Oct 2013
%
%  (c) UCLouvain 2018


% /!\ joint-sparse in columns

nxi = sqrt(sum(X.^2,1));
a = max(nxi-gamma,0);
a(a>0) = a(a>0)./nxi(a>0);
V = repmat(a,size(X,1),1).*X;


end
