function init = opImresize_init(n,p,kernel,kernel_width,antialiasing)
% OPIMRESIZE_INIT initialize the imresize operator compatible with the spot
% toolbox.
%
% Usage:
% init = OPIMRESIZE_INIT(n,p,kernel,kernel_width,antialiasing)
%   
% see also opImresize
% Author: Kevin Degraux
%  (c) UCLouvain 2018

if numel(n) > 2
    init.channels = n(3);
else
    init.channels = 1;
end

for dim = 1:2
    init.in_length(dim) = n(dim) ;
    init.out_length(dim) = n(dim)*p(dim) ;
    scale = p(dim);
    
    [init.weights{dim}, init.indices{dim}] = ...
        resize_contributions(init.in_length(dim), init.out_length(dim), ...
                             scale, kernel, ...
                             kernel_width, antialiasing);

    [init.weights_adj{dim}, init.indices_adj{dim}] = ...
        resize_contributions_adjoint(init.in_length(dim), init.out_length(dim), ...
                                     init.weights{dim}, init.indices{dim});
                                        
end