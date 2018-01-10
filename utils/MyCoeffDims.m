function [pext, qext, levels] = MyCoeffDims(p, q, levels)
% MYCOEFFDIMS
% Compute the size of the UDWT decomposition.
% Extracted from the opWavelet2 function from the spot toolbox.
% Usage:
% [pext, qext, levels] = MYCOEFFDIMS(p, q, levels)
%
% See also: opWavelet2
%
% Author: K. Degraux
%  (c) UCLouvain 2018
         if p >= 2^levels
            plevels = levels;
            if q >= 2^levels
               qext = ceil(q/(2^levels))*2^levels;
            elseif q > 1
               qlevels = floor(log2(q));
               levels = min(plevels,qlevels);
               qext = ceil(q/(2^levels))*2^levels;
            else
               qext = q;
            end
            pext = ceil(p/(2^levels))*2^levels;
         elseif p > 1
            plevels = floor(log2(p));
            if q >= 2^levels
               levels = min(levels,plevels);
               qext = ceil(q/(2^levels))*2^levels;
            elseif q > 1
               qlevels = floor(log2(q));
               levels = min(plevels,qlevels);
               qext = ceil(q/(2^levels))*2^levels;
            else
               levels = min(levels,plevels);
               qext = q;
            end
            pext = ceil(p/(2^levels))*2^levels;
         else
            pext = p;
            if q >= 2^levels
               qext = ceil(q/(2^levels))*2^levels;
            elseif q > 1
               levels = floor(log2(q));
               qext = ceil(q/(2^levels))*2^levels;
            else
               qext = q;
            end
         end
      end