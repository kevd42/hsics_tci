function [u,info] = solve_lys_syn_pos_PIHT(Phi,Phi_T_Phi,Phi_T_y,Psi,A,om,u0,param)
%% solve_lys_syn_pos_PIHT
% Solve Analysis-Synthesis with positivity constraint with the
% PIHT
% min_{u}  1/2 * norm(Phi*u - y)^2  s.t.  ||A*u||_0 <= K_end
%                                        0 <= Psi*u <= upperbnd
%
% USAGE
% solve_lys_syn_pos_PIHT(Phi_T_Phi,Phi_T_y,Psi,A,om,u0,mu,normPhi2,param)
%
% INPUTS : 
%
% Phi_T_Phi : linear operator computing Phi'*Phi*u for input u
% Phi_T_y : Data pre-multiplied by Phi'
% Psi : linear operator for the positivity constraint (synthesis op). Must
%       be orthonormal.
% A : linear operator for the analysis prior. Must be a tight frame such
%     that A'*A = Id; (NB: not mathematically necessary but implementation
%     requires it)
% om : optional weights for the rows of A applied before hard thresholding.
% u0 : initial guess for u (al0 is computed as al0 = A*u0)
% mu : Regularization parameter (cfr above)
% normPhi2 : norm of Phi_T_Phi or a good upper bound (Lipschitz constant of
%            the gradient of F)
% optional fields of param:
% verbosity : level of verbosity of the function. (0, 1 or 2)
% display_u : function handle to display the current iterate 
%                        every 10 iterations (must accept u as argument).
% maxiter
% tol
% iter_div
% upperbnd
% checkpoint
%
%
% Author: K. Degraux (ISPGroup)
% Date: 1 Sep 2015
%  (c) UCLouvain 2018


%% 
if ~isfield(param,'verbosity'); verbosity = 0;     else; verbosity = param.verbosity; end
if ~isfield(param,'maxiter');   maxiter   = 10000; else; maxiter   = param.maxiter;   end
if ~isfield(param,'tol');       tol       = 1e-9;  else; tol       = param.tol;       end
if ~isfield(param,'iter_div');  iter_div  = 50;    else; iter_div  = param.iter_div;  end
if ~isfield(param,'upperbnd');  upperbnd  = inf;   else; upperbnd  = param.upperbnd;  end
if ~isfield(param,'lowerbnd');  lowerbnd  = 0;     else; lowerbnd  = param.lowerbnd;  end
if ~isfield(param, 'K'); param.K = round([size(A,1)/(2*maxiter),size(A,1)/2]); end

if isfield(param,'checkpoint')
    checkpoint = param.checkpoint;
    if ~ischar(checkpoint)
        checkpoint = 'matfiles/checkpoint/default_checkpoint.mat';
    end
end
    

K_chunks = length(param.K)-1;
Kxx = [1,floor(maxiter/K_chunks)*(1:K_chunks-1),maxiter];
K = round(linterp(Kxx,param.K,1:maxiter));

%% Algorithm
rel_err = zeros(maxiter/iter_div,1);
info    = struct;
t_init  = 1;
if isfield(param,'qfun')
    info.qfun_out = [];
end

u = u0;
om = abs(om);

if isfield(param,'checkpoint') && exist(checkpoint,'file')
    load(checkpoint);
    fprintf('Checkpoint loaded, restarting at iteration %i\n',t_init);
end

if verbosity>1
    timer = tic;
end
if isfield(param,'display_u')
    hFig = pauseFigure(42,'init');
end
for t = t_init:maxiter
    
        
    uprev = u;
        
    tmp = Phi_T_y - Phi_T_Phi*u;
    
    norm_tmp_sq = norm(tmp)^2;
    if norm_tmp_sq>0
        tau = norm_tmp_sq / norm(Phi*tmp)^2;
    else
        tau = 0;
    end
    
    tmp = u + tau*tmp;

    b = A*tmp;
    
    idx = om~=0;
    b(idx) = Hs(b(idx)./om(idx),K(t)).*om(idx);
    
    

    u = A'*b;
    
    u = Psi'*min(max(Psi*u,lowerbnd),upperbnd);
    

   
    if verbosity>2
        fprintf('%5i ',t);
    end
    
    if isfield(param,'display_u')
        hFig = pauseFigure(hFig);
        drawnow;
    end
    if mod(t,iter_div)==0
        rel_err(t/iter_div) = norm(u-uprev)/norm(u);
        if verbosity>0
            fprintf('\n %5i  : \t relative error = %g \n',t,rel_err(t/iter_div));
        end
        if isfield(param,'display_u') && isa(param.display_u,'function_handle')
            param.display_u(u);
            drawnow;
        end
        
        if isfield(param,'qfun')
            info.qfun_out = param.qfun(u,info,t/iter_div);
        end
       
        
        if rel_err(t/iter_div) < tol
            break;
        end
        if verbosity>1
            teniter_time = toc(timer);
            timer=tic;
            print_ETA(teniter_time,maxiter-t,iter_div);
        end
        
        t_init = t+1; %#ok<NASGU>
        if isfield(param,'checkpoint')
            if exist(checkpoint,'file')
                copyfile(checkpoint,[checkpoint,'.bak']);
            end
            save(checkpoint,'t_init','u','rel_err','info');
            
        end
        
    end
end

info.rel_err  = rel_err;
info.niter    = t;

