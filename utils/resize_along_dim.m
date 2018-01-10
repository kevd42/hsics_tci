function out = resize_along_dim(in, dim, weights, indices)
% RESIZE_ALONG_DIM Resize along a specified dimension
% Usage:
% out = RESIZE_ALONG_DIM(in, dim, weights, indices)
% 
% in           - input array to be resized
% dim          - dimension along which to resize
% weights      - weight matrix; row k is weights for k-th output pixel
% indices      - indices matrix; row k is indices for k-th output pixel
% 
% See also IMRESIZE IMRESIZEMEX RESIZE_CONTRIBUTIONS
% 
% Copied from the image processing toolbox imresize function.
% Renamed to avoid conflicts.
% Original name : resizeAlongDim
% Copyright 1992-2013 The MathWorks, Inc.
%
% 

% if isPureNearestNeighborComputation(weights)
%     out = resizeAlongDimUsingNearestNeighbor(in, dim, indices);
%     return
% end

out_length = size(weights, 1);

size_in = size(in);
size_in((end + 1) : dim) = 1;

if (ndims(in) > 3)
    % Reshape in to be a three-dimensional array.  The size of this
    % three-dimensional array is the variable pseudo_size_in below.
    %
    % Final output will be consistent with the original input.
    pseudo_size_in = [size_in(1:2) prod(size_in(3:end))];
    in = reshape(in, pseudo_size_in);
end

% The 'out' will be uint8 if 'in' is logical 
% Otherwise 'out' datatype will be same as 'in' datatype
out = imresizemex(in, weights', indices', dim);

if ( (length(size_in) > 3) && (size_in(end) > 1) )
    % Restoring final output to expected size
    size_out = size_in;
    size_out(dim) = out_length;
    out = reshape(out, size_out);
end
%---------------------------------------------------------------------
