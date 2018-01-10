function [Phi,info] = constructInFocus(m,n,idLy,SR)
% Phi = CONSTRUCTINFOCUS(slm,psf,idLy,SR,m,n,s)
% Pre-Compute useful arrays and return an opFunction
%  INPUTS:
%   idLy : cell of logical arrays of size ms that contain, for each band,
%          the location of the corresponding filters on the FPA.
%   SR : Spectral response matrix to be applied on every pixel spectrum.
%        Size : nl x nf (nl = n(3) = number of wl, nf = number of filters)
%        If empty, defaults as 1.
%   m : output size [mi, mj, mS] 
%       with ms = [mi mj] spatial dimensions of the FPA.
%       and  mS = number of snapshots.
%       NB : mS should always be 1
%   n : input size [ni, nj, nl]
%       with ns = [ni, nj] spatial dimensions of the target image.
%       and  nl = spectral dimension (number of wl).
%  OUTPUT:
%   Phi : spot opFunction that models the InFocus device behavior.
%
%
%  NOTE : This function should work in the Panchromatic case as well, 
%         setting n(3)=1
%
% Author: Kevin Degraux
%  (c) UCLouvain 2018


ns    = n(1:2);
ms    = m(1:2);
mS    = m(3);

Ns = prod(ns);
Ms = prod(ms);

if isscalar(SR)
    nf = n(3);
else
    error('SR not yet supported');
    nf = size(SR,2);
end

if isequal(ms,ns)
    Up = opEye(Ns);
    info.invPhi_T_Phi_rhoId_h = @(rho) opBlockDiag(nf,(1/(1+rho))*opEye(Ns),0);
else
    init = opImresize_init(ns,ms./ns,@(x) lanczos(x,3),6,false);
    Up = opFunction(Ms,Ns,@(in,mode) opImresize(in,mode,init));
    % dUp = double(Up); % sparse matrix version of the Up operator;
    partUp = cell(2,1);
    for dim = 1:2
        partUp{dim} = spalloc(ms(dim),ns(dim),6*ms(dim));
        for i=1:6
            idx = sub2ind([ms(dim),ns(dim)],(1:ms(dim))',init.indices{dim}(:,i));
            partUp{dim}(idx) = partUp{dim}(idx) + init.weights{dim}(:,i);
        end
    end
    UpTUp = opKron(partUp{1}'*partUp{1},partUp{2}'*partUp{2});
    %dUpTUp = kron(partUp{1}'*partUp{1},partUp{2}'*partUp{2});
    % To improve speed, implement custom opSpot with internal memory and
    % hotstart pcg. Maybe later (already reasonably fast).
    info.invPhi_T_Phi_rhoId_h = @(rho) opBlockDiag(nf,...
             opFunction(Ns,Ns,@(in,mode) pcg_quiet( UpTUp+rho*opEye(Ns), in)), 0); 
end

FilterMaps = cellfun(@(i) opMask(i)*Up,idLy,'UniformOutput', false);
FilterMap = [FilterMaps{:}];

Phi  = FilterMap;
fpa_idx = reshape(1:(Ms*nf),[ms,nf]);
info.Rm_idx = zeros(ms);
for band = 1:nf
    fpa_idx_band = fpa_idx(:,:,band);
    info.Rm_idx(idLy{band}) = fpa_idx_band(idLy{band});
end
info.Rm_idx = info.Rm_idx(:);
info.Rn_idx = ':'; % identity

info.FTDF = opBlockDiag(nf,Up,0); %

end

% function invdUpTUp_rhoId = precomp_invdUpTUp_rhoId(rho,Ns,UpTUp)
%     % Precompute the matrix to invert (depending on rho) and the
%     % preconditionner to speed up the conjugate gradient (cfr. pcg).
%     % Use opSpot implementation of UpTUp.
%     % Observation : Using preconditionning slows down the process
%     
% %     A = dUpTUp + rho*speye(Ns); % Matrix to invert (symmetric PSD).
% %     L = ichol(A); % Incomplete Cholesky factorization for preconditionning
%     invdUpTUp_rhoId = @(in,mode) pcg( UpTUp+rho*opEye(Ns), in, [],[]);
% end

function [x,flag,relres,iter,resvec] = pcg_quiet(varargin)
% Simple wrapper to suppress the output message of pcg.
[x,flag,relres,iter,resvec] = pcg(varargin{:});
end
