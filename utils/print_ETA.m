function print_ETA(niter_time,rem_iter,n)
% PRINT_ETA
% Print the remaining time, given the time taken by n iterations and the
% number of iterations remaining.
% Usage:
%  print_ETA(niter_time,rem_iter,n)
% Inputs:
% niter_time = time taken by n iterations
% rem_iter = number of iterations remaining
% n = number of iterations between tic and toc
% Author: K.Degraux (ISPGroup)
% Date: 21 Aug 2015
%  (c) UCLouvain 2018
rem_sec = ceil(rem_iter * niter_time/n);
rem_min = floor(rem_sec/60);
rem_h   = floor(rem_min/60);
rem_d   = floor(rem_h/24);

rem_sec = rem_sec-60*rem_min;
rem_min = rem_min-60*rem_h;
rem_h   = rem_h-24*rem_d;


fprintf('ETA : %i days %2i:%2i:%2i \n',rem_d,rem_h,rem_min,rem_sec);
end