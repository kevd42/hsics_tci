function [u,alg_info] = launch_ABPDN_box_ADMM(settings,rho)
% LAUNCH_ABPDN_BOX_ADMM launch the analysis-BPDN with box
% contraint ADMM algorithm.
%
% Usage:
% [u,alg_info] = LAUNCH_ABPDN_BOX_ADMM(settings,rho)
%
% Author: Kevin Degraux
%
% see README.txt
%
% see also setenv, hsics_setup, hsics_jobid, hsics_settings,
% hsics_reconstruction, hsics_gather_psnr, hsics_psnr_figure,
% hsics_patches_figure
%
%  (c) UCLouvain 2018

n     = settings.n;
m     = settings.m;
Phi   = settings.sensing.Phi;

% ns = n(1:2);
% nf = n(3);
% ms = m(1:2);
% mS = m(3);

if isfield(settings.sensing,'Phi_info')
    Phi_info = settings.sensing.Phi_info;
else
    Phi_info = struct();
end

% if isfield(Phi_info,'s')
%     s = Phi_info.s;
% else
%     s = n(1:2)+m(1:2)-1;
% end

if isfield(Phi_info,'FTDF')
	FTDF = Phi_info.FTDF;
else
    error('please define FTDF');
end
if isfield(Phi_info,'Rm_idx')
	Rm_idx = Phi_info.Rm_idx;
else
    error('please define Rm_idx');
end
if isfield(Phi_info,'Rn_idx')
	Rn_idx = Phi_info.Rn_idx;
else
    error('please define Rn_idx');
end
if isfield(Phi_info,'normPhi2')
	normPhi2 = Phi_info.normPhi2;
else
    normPhi2 = normest(FTDF,1e-3)^2;
end
if isfield(Phi_info,'invPhi_T_Phi_rhoId_h')
	invPhi_T_Phi_rhoId_h = Phi_info.invPhi_T_Phi_rhoId_h;
else
    invPhi_T_Phi_rhoId_h = @(rho) inv( FTDF'*FTDF + rho * opEye(nn) );
end
    
    
y         = settings.measurements.y;
epsilon   = settings.measurements.epsilon;
A         = settings.prior.A;
if isfield(settings.prior,'A_info') && isfield(settings.prior.A_info,'om')
    %alpha = 0;% Original L2-normed wavelet
    %alpha = 1;% tight frame normed
    %alpha = 2;% Besov with s=1 ??
    alpha = 2;
    om = 1./(settings.prior.A_info.om.^(alpha-1));
    om = om/mean(om); % For numerical stability. Does not change solution.
else
    om = 1;
end

u0 = settings.measurements.u0;
u0_info = settings.measurements.u0_info;

yc0 = u0_info.yc0;

param.upperbnd = settings.upperbnd;

param.lowerbnd = settings.lowerbnd;

param.verbosity = 2;

if ~isdeployed
    param.display_u = @(u) displayIterVolume(1,u,n,4);
    graphics = true;
else
    graphics = false;
end
    
param.checkpoint = ['matfiles/checkpoint/',settings.identifier,...
                    '/cp_',settings.filename,'.mat'];


param.maxiter  = 2000;
param.tol      = 5e-5;
param.iter_div = 10;

norm_y     = norm(y);
if settings.groundtruth.available
    x          = settings.groundtruth.x;
    if isfield(settings.groundtruth,'dyn')
        dyn = settings.groundtruth.dyn;
    else
        dyn = 1;
    end
    itmax = param.maxiter/param.iter_div;
    param.qfun = @(u_est,info,it) printPSNRs(u_est,info,it,itmax,...
               param.iter_div,y,norm_y,Phi,1,A,om,x,dyn,true,graphics,123);
else
    itmax = param.maxiter/param.iter_div;
    param.qfun = @(u_est,info,it) printPSNRs(u_est,info,it,itmax,...
              param.iter_div,y,norm_y,Phi,[],A,om,[],[],true,graphics,123);
end


param.rho = rho;


% TODO : Infer the missing information for a better initialization of ADMM
%yc0 = zeros([...]);


[u,alg_info] = solve_ABPDN_box_ADMM(FTDF, Rm_idx, Rn_idx, y, yc0,...
                                    epsilon, invPhi_T_Phi_rhoId_h, ...
                                    A, om, u0, normPhi2, param);

end