function [u,alg_info] = launch_lys_syn_pos_PIHT(settings,K)
% LAUNCH_LYS_SYN_POS_PIHT launch the analysis-synthesis with positivity
% contraint PIHT algorithm
%
% Usage:
% [u,alg_info] = launch_lys_syn_pos_PIHT(settings,K)
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

Phi       = settings.sensing.Phi;
if isfield(settings.sensing,'Phi_info') && isfield(settings.sensing.Phi_info,'PhiTPhi')
    Phi_T_Phi = settings.sensing.Phi_info.PhiTPhi;
else
    Phi_T_Phi = Phi'*Phi;
end
y         = settings.measurements.y;
Phi_T_y   = Phi'* y;
Psi       = settings.prior.Psi;
Psi0      = settings.prior.Psi0;
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


param.upperbnd = settings.upperbnd;

param.lowerbnd = settings.lowerbnd;

param.verbosity = 2;

if ~isdeployed
    param.display_u = @(u) displayIterVolume(Psi0,u,n,4);
    graphics = true;
else
    graphics = false;
end
    
param.checkpoint = ['matfiles/checkpoint/',settings.identifier,'/cp_',settings.filename,'.mat'];


param.maxiter  = 1000;
param.tol      = 1e-6;
param.iter_div = 10;

norm_y     = norm(y);
if settings.groundtruth.available
    x          = settings.groundtruth.x;
    norm_x     = norm(x);
    itmax = param.maxiter/param.iter_div;
    param.qfun = @(u_est,info,it) printPSNRs(u_est,info,it,itmax,param.iter_div,y,norm_y,Phi,Psi0,A,om,x,norm_x,true,graphics,123);
else
    itmax = param.maxiter/param.iter_div;
    param.qfun = @(u_est,info,it) printPSNRs(u_est,info,it,itmax,param.iter_div,y,norm_y,Phi,[],A,om,[],[],true,graphics,123);
end


param.K = K*size(A,1)/100;

[u,alg_info] = solve_lys_syn_pos_PIHT(Phi,Phi_T_Phi,Phi_T_y,Psi,A,om,u0,param);


end