%% IMRESIZEMEX Unofficial help of built-in function IMRESIZEMEX
% Usage :
% out = IMRESIZEMEX(in, weights, indices, dim);
%
% is equivalent to
% 
% if dim ==2 
%   in = in.';
% end
%
% out_size = size(indices,2);
% in_size  = size(in,2);
% 
%  
% out = zeros(out_size,in_size);
% for i=1:out_size
%   out(i,:) = (in(indices(:,i),:)*weights(:,i)).' ;
% end
%
% if dim ==2 
%   out = out.';
% end
% 
% with limited precision but fast.
% 
% Columns of indices contain the indices of the pixel rows (col if dim==2) 
% of in that participate to the corresponding row (col if dim==2) of out.
% Columns of weights contain the weights of the pixel rows (col if dim==2) 
% of in that participate to the corresponding row (col if dim==2) of out.
%
% See also RESIZE_ALONG_DIM IMRESIZE RESIZE_CONTRIBUTIONS
% 
% Author of the Help M-file : K. Degraux
% 
% Mex function Copyright 1992-2013 The MathWorks, Inc.
%
%