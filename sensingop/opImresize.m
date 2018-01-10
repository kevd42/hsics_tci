function out = opImresize(in,mode,init)
% OPIMRESIZE Spot operator using imresize.
%
% Usage:
% out = opImresize(in,mode,init)
%   
% see also: opImresize_init
% Author: Kevin Degraux
%  (c) UCLouvain 2018

if issparse(in)
    in = full(in);
end



if mode == 1
    tmp = reshape(in,[init.in_length,init.channels]);
    for dim = 1:2
        tmp = resize_along_dim(tmp, dim, init.weights{dim}, init.indices{dim});
    end
else
    tmp = reshape(in,[init.out_length,init.channels]);
    for dim = 1:2
        tmp = resize_along_dim(tmp, dim, init.weights_adj{dim}, init.indices_adj{dim});
    end
end
out = tmp(:);

end