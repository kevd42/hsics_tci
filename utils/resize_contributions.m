function [weights, indices] = resize_contributions(in_length, out_length, ...
                                            scale, kernel, ...
                                            kernel_width, antialiasing)
% RESIZE_CONTRIBUTIONS Computes weights and indices of an input pixels to
% resize an image.
%
% Usage :
% [weights, indices] = resize_contributions(in_length, out_length, ...
%                                             scale, kernel, ...
%                                             kernel_width, antialiasing)
% Outputs :
% weights      - weight matrix; row k is weights for k-th output pixel
% indices      - indices matrix; row k is indices for k-th output pixel
% 
% Inputs :
% in_length    - length of input array along specific dimension
% out_length   - length of output array along specific dimension
% scale        - ratio out_length/in_length
% kernel       - function handle returning a 1-D interpolation kernel
% kernel_width - width of the kernel (in pixels)
% antialiasing - logical. If true and scale<1, modify the kernel to
%                simulatneously interpolate and antialias
%
% See also IMRESIZE IMRESIZEMEX RESIZE_ALONG_DIM
% 
% Copied from the image processing toolbox imresize function.
% Renamed to avoid conflicts.
% Original name : contributions
% Help written by K. Degraux
% Copyright 1992-2013 The MathWorks, Inc.
%
%

if (scale < 1) && (antialiasing)
    % Use a modified kernel to simultaneously interpolate and
    % antialias.
    h = @(x) scale * kernel(scale * x);
    kernel_width = kernel_width / scale;
else
    % No antialiasing; use unmodified kernel.
    h = kernel;
end

% Output-space coordinates.
x = (1:out_length)';

% Input-space coordinates. Calculate the inverse mapping such that 0.5
% in output space maps to 0.5 in input space, and 0.5+scale in output
% space maps to 1.5 in input space.
u = x/scale + 0.5 * (1 - 1/scale);

% What is the left-most pixel that can be involved in the computation?
left = floor(u - kernel_width/2);

% What is the maximum number of pixels that can be involved in the
% computation?  Note: it's OK to use an extra pixel here; if the
% corresponding weights are all zero, it will be eliminated at the end
% of this function.
P = ceil(kernel_width) + 2;

% The indices of the input pixels involved in computing the k-th output
% pixel are in row k of the indices matrix.
indices = bsxfun(@plus, left, 0:P-1);

% The weights used to compute the k-th output pixel are in row k of the
% weights matrix.
weights = h(bsxfun(@minus, u, indices));

% Normalize the weights matrix so that each row sums to 1.
weights = bsxfun(@rdivide, weights, sum(weights, 2));

% Clamp out-of-range indices; has the effect of replicating end-points.
indices = min(max(1, indices), in_length);

% If a column in weights is all zero, get rid of it.
kill = find(~any(weights, 1));
if ~isempty(kill)
    weights(:,kill) = [];
    indices(:,kill) = [];
end
