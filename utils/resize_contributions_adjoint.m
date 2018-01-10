function [weights_adj, indices_adj] = resize_contributions_adjoint(in_length, out_length, ...
                                            weights, indices)
% RESIZE_CONTRIBUTIONS_ADJOINT                                  
% Construct the sparse matrix associated to the adjoint of the resize
% operator.
% Usage:
% [weights_adj, indices_adj] = RESIZE_CONTRIBUTIONS_ADJOINT(in_length,...
%                                                           out_length,... 
%                                                           weights,...
%                                                           indices) 
% see also resize_contribtuions
% 
% Author : K. Degraux
%  (c) UCLouvain 2018

ind_width = ceil(size(indices,2)*out_length/in_length);
                                        

[i,~,v]=find(indices);

[i_sort,ord] = sort(i);

v_sort = v(ord);
weights_sort = weights(ord);

indices_adj = ones(in_length,ind_width);
weights_adj = zeros(in_length,ind_width);

for k = 1:in_length
    tmp_w = weights_sort(v_sort==k);
    tmp_i = i_sort(v_sort==k);
    i = 1;
    while length(tmp_i)>ind_width
        if tmp_i(i) == tmp_i(i+1)
            tmp_i(i) =[];
            tmp_w(i+1) = tmp_w(i)+tmp_w(i+1);
            tmp_w(i) = [];
        else
            i = i+1;
        end
    end
    compl = ind_width - length(tmp_i);
    indices_adj(k,:) = [tmp_i;ones(compl,1)].';
    weights_adj(k,:) = [tmp_w;zeros(compl,1)].';
end

end