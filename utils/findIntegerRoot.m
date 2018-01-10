function [P,Q] = findIntegerRoot(N)
% FINDINTEGERROOT
% Find integer factorization of N closest to square root.
% [P,Q] = FINDINTEGERROOT(N)
% So that P*Q==N
% Author: K. Degraux
%  (c) UCLouvain 2018
    P = floor(sqrt(N));
    while(mod(N,P)~=0)
        P=P-1;
    end
    Q = N/P;
end