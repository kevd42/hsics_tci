
function [x,suppx] = Hs(a,s)
% Hs(a,s)
% a and x are vectors
% s is the number of non zero elements in x
% non linear operator taht sets all but the largest (in magnitude) s
% elements of a to zero. If there is no unique such set, a set can be
% selected either randomly or based on a predefined ordering of the 
% elements. 
% Author : K. Degraux
% Date : April 2012
%  (c) UCLouvain 2018
    x = a;
    [~,ix] = sort(abs(a),'descend');
    x(ix(s+1:end)) = 0;
    if (nargout == 2)
        suppx = ix(1:s);
    end
end

% V2 : hard thresholding
% lambda0p5
% function x = Hlambda(a,lambda0p5)
%     x = a;
%     x(abs(x)<lambda0p5) = 0;
% end
