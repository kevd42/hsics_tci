function [u0,u0_info] = computeInitFP(info,y,epsilon,idL,optical_path,m,n)
% COMPUTEINITFP Compute the linear interpolation of the known samples on
% the FP spectral sensor (see (8) in the paper) then compute the tikhonov 
% regularization and return the initial volume.
% 
%
% Usage:
% [u0,u0_info] = COMPUTEINITFP(info,y,epsilon,idL,optical_path,m,n)
%
% 
% Author: Kevin Degraux
% Date: 2017
%
%
%  (c) UCLouvain 2018



ms = m(1:2);
mi = m(1);
mj = m(2);
mS = m(3);
ns = n(1:2);
ni = n(1);
nj = n(2);
nl = n(3);

if strcmp(optical_path,'OffFocus')
    % TODO: put this code snippet where it should be...
    s = ms+ns-1;
    si = s(1);
    sj = s(2);
    mm = prod(s)*mS*nl;
    nn = prod(s)*nl;
    % Valid part of the convolution
    idx_valid = {(ni : si);
                 (nj : sj)}; 
    % Central part of the valid part
    offset    = floor((s-ns+1)/2 - ms/2 );
    fpa_crop  = {(1:mi) + offset(1);
                 (1:mj) + offset(2)}; 
    % Combine the cropping regions
    fpa_idx = {idx_valid{1}(fpa_crop{1});
               idx_valid{2}(fpa_crop{2})};
    y0    = zeros([s,mS,nl]);
elseif strcmp(optical_path,'InFocus')
    mm = prod(ms)*nl;
%     nn = prod(ms)*nl;
%     s = ms;
    fpa_idx = {1:mi;
               1:mj};
    y0    = zeros([m,nl]);
end

y_fpa = reshape(y,m);

for snap = 1:mS
    y_snap = y_fpa(:,:,snap);
    y0(fpa_idx{:},snap,:)  = reshape(interpolateFPA(y_snap(:),ms,nl,idL),[ms,1,nl]);
end

Rm_idxc = 1:mm;
Rm_idxc(info.Rm_idx) = [];
u0_info.yc0 = y0(Rm_idxc);

% Tikhonov regularization with unmasked unpadded convolution
% v0 = info.invPhi_T_Phi_rhoId_h(epsilon^2) * info.FTDF'*y0(:); 
% 
% u0 = v0(info.Rn_idx);


if strcmp(optical_path,'OffFocus')
    idx = reshape(1:mm,[s,mS,nl]);
    Rm_test_idx = [];
    for snap = 1:mS
        for band = 1:nl
            Rm_test_snap_band = idx(fpa_idx{:},snap,band);
            Rm_test_idx = [Rm_test_idx;Rm_test_snap_band(:)];
        end
    end
    Rm_test = opRestriction(mm,Rm_test_idx);
    Phi_test = Rm_test*info.FTDF*opRestriction(nn,info.Rn_idx)';
    tic;
    u0 = pcg((Phi_test'*Phi_test + epsilon^2 * opEye(prod(ns)*nl)),Phi_test'*y0(Rm_test_idx),1e-3,10);
    toc;
elseif strcmp(optical_path,'InFocus')
    u0 = info.invPhi_T_Phi_rhoId_h(epsilon^2) * info.FTDF'*y0(:);
end



u0 = max(min(u0,1),0);
u0 = u0-min(u0);
u0 = u0/max(u0);

end