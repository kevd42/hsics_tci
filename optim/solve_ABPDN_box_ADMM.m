function [u,info] = solve_ABPDN_box_ADMM(Phi, Rm_idx, Rn_idx, y, yc0,...
                                         epsilon, invPhi_T_Phi_rhoId_h, ...
                                         A, om, u0, normPhi2, param)
%% solve_ABPDN_box_ADMM
% Solve ABPDN with box constraint with ADMM
% min_{u} ||om*A*u||_1 s.t.  norm(Rm*Phi*Rn'*u - y) <= epsilon
%                            lowerbnd <= u <= upperbnd
%
% which, with v = Rn'*u, is equivalent to
%
% min_{v} ||[om;ones].*[A*Rn ; Rn_c]*v||_1  s.t.
%                                        norm(Rm*Phi*v - y) <= epsilon
%                                        lowerbnd <= Rn*v <= upperbnd
%                                        Rn_c*v == 0
%
% Where 
%        u is the original optimization variable
%        v is the zero-padded optimization variable
%        om is a vector of weights
%        Rn is a restriction operator of size n x nn
%        Rn_c is the complementary restriction (of size nn-n x nn)
%        Rm is a restriction operator of size m x mm
%        Phi is a sensing operator of size mm x nn
%        y is a data vector of size m
%        epsilon, lowerbnd and upperbnd are scalar bounds (applied
%        element-wise)
%
% USAGE
% [u,info] = solve_ABPDN_box_ADMM(Phi, Rm, Rn, y, epsilon, ...
%                                invPhi_T_Phi_rhoId,A,om,u0,normPhi2,param)
%
% INPUTS : 
% Phi : linear operator. Must support matvec and conjugate transposition.
% Rm_idx  : Indices of the restriction operator Rm (see above).
%           Can include permutations.
% Rn_idx  : Indices of the restriction operator Rn (see above).
%           Can include permutations.
% y   : data vector
% yc0 : complementary (interpolated) data vector (unobserved data). 
%       If empty, the default is 0.
% epsilon : scalar bound for the ABPDN constraint
% invPhi_T_Phi_rhoId_h : function handle with input rho that returns a
%   linear operator (input in) which computes the solution (out) to
%   (Phi'*Phi + rho*Id) * out = in
% A   : linear operator for the analysis prior. Must be a tight frame such
%       that A'*A = Id;
% om  : weights in the l1 norm (cfr above)
% u0  : initial guess for u
% normPhi2 : norm of Phi_T_Phi or a good upper bound (Lipschitz constant of
%            the gradient of F)
% optional fields of param:
% verbosity : level of verbosity of the function. (0, 1 or 2)
% display_u : function handle to display the current iterate 
%             every iter_div iterations (must accept u as argument).
% maxiter  : maximum number of iterations.
% tol      : tolerance on relative increment between two iterations.
% iter_div : number of iterations between two computations of convergence
% lowerbnd, upperbnd : scalar bounds for the box constraint
% rho : parameter controlling the convergence speed of ADMM. 
% checkpoint : path to the checkpoint file
%
%
% Author: K. Degraux (ISPGroup)
% Date: 13 Apr 2017
%  (c) UCLouvain 2018


%% 
if ~isfield(param,'verbosity'); verbosity = 0;     else; verbosity = param.verbosity; end
if ~isfield(param,'maxiter');   maxiter   = 10000; else; maxiter   = param.maxiter;   end
if ~isfield(param,'tol');       tol       = 1e-9;  else; tol       = param.tol;       end
if ~isfield(param,'iter_div');  iter_div  = 50;    else; iter_div  = param.iter_div;  end
if ~isfield(param,'upperbnd');  upperbnd  = inf;   else; upperbnd  = param.upperbnd;  end
if ~isfield(param,'lowerbnd');  lowerbnd  = 0;     else; lowerbnd  = param.lowerbnd;  end
if ~isfield(param,'rho');       rho       = 20;    else; rho       = param.rho;  end

if isfield(param,'checkpoint')
    checkpoint = param.checkpoint;
    if ~ischar(checkpoint)
        checkpoint = 'matfiles/checkpoint/default_checkpoint.mat';
    end
end
%% Algorithm
rel_err = zeros(maxiter/iter_div,1);
info    = struct;
t_init  = 1;
if isfield(param,'qfun')
    info.qfun_out = [];
end

% ADMM constants (must be >0 for convergence)
mu1 = rho*50/normPhi2;
mu2 = rho;
mu3 = rho;

nu       = 2;
tau_incr = 1.5;
tau_decr = 1.5;

alpha = 1.8;

% Construction of the inverse operator
invPhi_T_Phi_rhoId = (1/mu1) * invPhi_T_Phi_rhoId_h( (mu2+mu3)/mu1 );

% Sizes of the operators
[mm,nn] = size(Phi);
if Rn_idx == ':'
    n = nn;
else
    n = length(Rn_idx);
end
if Rm_idx == ':'
    m = mm;
else
    m = length(Rm_idx);
end
    

p = size(A,1);

% Complementary index maps
Rn_idxc = 1:nn;
Rn_idxc(Rn_idx) = [];

Rm_idxc = 1:mm;
Rm_idxc(Rm_idx) = [];

% Initialization
nu1 = zeros(mm,1);
nu1(Rm_idx)   = y;
if ~isempty(yc0); nu1(Rm_idxc)  = yc0; end
d1 = zeros(mm,1);


% nu2(p+1:end) is always zero so we redefine nu2 = nu2(1:p)
nu2 = A*u0;
d2 = zeros(p+nn-n,1);

% nu3(Rn_idxc) is always zero so we redefine nu3 = nu3(Rn_idx)
nu3 = u0;

d3 = zeros(nn,1);

v = mu1 * Phi'*nu1;
v(Rn_idx) = v(Rn_idx) + mu2 * A'*nu2 + mu3 * nu3;

v = invPhi_T_Phi_rhoId * v;



% Residuals (commented out because we only need the norm)
% r1 defined inside the loop
% r2 = zeros(p+nn-n,1);
% r3 = zeros(nn,1);
% 
% r1 = zeros(maxiter,1);
% r2 = zeros(maxiter,1);
% r3 = zeros(maxiter,1);
% s1 = zeros(maxiter,1);
% s2 = zeros(maxiter,1);
% s3 = zeros(maxiter,1);



if isfield(param,'checkpoint') && exist(checkpoint,'file')
    load(checkpoint);
    fprintf('Checkpoint loaded, restarting at iteration %i\n',t_init);
    
end

Hv1 = Phi*v;
Hv2 = A*v(Rn_idx);

if verbosity>1
    timer = tic;
end
if isfield(param,'display_u')
    hFig = pauseFigure(42,'init');
end
for t = t_init:maxiter
    
    % ==== Data fidelity constraint ====

    ze1 = alpha*Hv1 + (1-alpha)*nu1 - d1;
    %ze1 = Hv1 - nu1 + ze1;
    
    %nu1_prev = nu1;
    % Projection of the selected entries on epsilon*B2(y)
    res = ze1(Rm_idx) - y;
    norm_res = norm(res);
    if norm_res > epsilon
        nu1(Rm_idx) = y + res*(epsilon/norm_res)  ;
    else
        nu1(Rm_idx) = ze1(Rm_idx);
    end
    nu1(Rm_idxc) = ze1(Rm_idxc);
    
    d1 = nu1 - ze1; % = d1 + nu1 - Hv1;
    
    % ==== Constrained weighted soft thresholding ====
    % instead of minimizing the L1-norm over the entire vector, we 
    % explicitly add the "zero-padding" constraint, which is also
    % equivalent to putting "infinite" weights on the border entries.
    % NB : We have the right to do that since the L1 norm is separable.
    
    ze2     = alpha*Hv2 + (1-alpha)*nu2 - d2(1:p);
    %ze2(p+1:end) = v2(Rn_idxc) - d2(p+1:end);
    
%     ze2(1:p)     = Hv2         + ze2(1:p)     - nu2;
%     ze2(p+1:end) = v(Rn_idxc)  + ze2(p+1:end) ;
    
    %nu2_prev = nu2;
    nu2 = om .* wthresh((1./om) .* ze2,'s',1/mu2);
    
    d2(1:p)     = nu2 - ze2;
    if ~isempty(Rn_idxc)
        d2(p+1:end) = d2(p+1:end) - (alpha * v(Rn_idxc)); % = - ze2(p+1:end);
    end
    % ==== Box constraint and "zero-padding" constraint again ====
    
    ze3 = alpha*v(Rn_idx)+(1-alpha)*nu3 - d3(Rn_idx);
    
    %ze3 = ze3 + v; 
    %ze3(Rn_idx) = ze3(Rn_idx) - nu3;
    
    %nu3_prev = nu3;
    nu3 = min(max( ze3 , lowerbnd ), upperbnd);
    
    d3(Rn_idx)  = nu3 - ze3;
    if ~isempty(Rn_idxc)
        d3(Rn_idxc) = d3(Rn_idxc) - alpha*v(Rn_idxc);    % - ze3(Rn_idxc);
    end
    % ==== Update step ====
    v_prev    = v;
    v         =               mu1 * Phi'*(d1         + nu1);
    v(Rn_idx) = v(Rn_idx)   + mu2 * A'  *(d2(1:p)    + nu2) ...
                            + mu3 *      (d3(Rn_idx) + nu3);
    if ~isempty(Rn_idxc)
        v(Rn_idxc) = v(Rn_idxc) + mu2 *       d2(p+1:end) ...
                                + mu3 *       d3(Rn_idxc);
    end

    v = invPhi_T_Phi_rhoId * v;
    

    
    % ==== Primal Residuals ====
    Hv1 = Phi * v ;
%     %r1 = nu1 - Phi * v;
%     r1(t) = norm(nu1 - Hv1);
%     
    Hv2 = A*v(Rn_idx);
%     %r2(1:p)     = nu2(1:p) - A * v(Rn_idx);
%     %r2(p+1:end) =          -     v(Rn_idxc);
%     r2(t) = sqrt(sum((nu2(1:p) - Hv2).^2)+sum(v(Rn_idxc).^2));
%     
%     %r3(Rn_idx)  = nu3 - v(Rn_idx);
%     %r3(Rn_idxc) =     - v(Rn_idxc);
%     r3(t) = sqrt(sum((nu3 - v(Rn_idx)).^2) + sum((- v(Rn_idxc)).^2));
%     
%     
%     %r = [r1;r2;r3];
%     
%     if mod(t,iter_div)==0
%         % ==== Dual Residuals ====
%         tic;
%         s1(t) = norm(mu1 * Phi' * (nu1_prev - nu1));
%         s2(t) = norm(mu2 * A'   * (nu2_prev - nu2));
%         s3(t) = norm(mu3 *        (nu3_prev - nu3));
%         toc;
% 
%         %s = s1+s2+s3;
% 
%         flag_mu = false;
%         mu1_old = mu1;
%         if r1(t)/2 + r2(t)/4 + r3(t)/4 > s1(t)*nu
%             mu1 = mu1 * tau_incr;
%             flag_mu = true;
%         elseif r1(t)/2 + r2(t)/4 + r3(t)/4 < s1(t)/nu
%             mu1 = mu1 / tau_decr;
%             flag_mu = true;
%         end
% 
%         mu2_old = mu2;
%         if r2(t)/2 + r3(t)/4 + r1(t)/4 > s2(t)*nu
%             mu2 = mu2 * tau_incr;
%             flag_mu = true;
%         elseif r2(t)/2 + r3(t)/4 + r1(t)/4 < s2(t)/nu
%             mu2 = mu2 / tau_decr;
%             flag_mu = true;
%         end
% 
%         mu3_old = mu3;
%         if r3(t)/2 + r1(t)/4 + r2(t)/4 > s3(t)*nu
%             mu3 = mu3 * tau_incr;
%             flag_mu = true;
%         elseif r3(t)/2 + r1(t)/4 + r2(t)/4 < s3(t)/nu
%             mu3 = mu3 / tau_decr;
%             flag_mu = true;
%         end
% 
%         if flag_mu
%             % update of the inverse operator
%             invPhi_T_Phi_rhoId = (1/mu1) * invPhi_T_Phi_rhoId_h( (mu2+mu3)/mu1 );
%             d1 = mu1_old /mu1 * d1;
%             d2 = mu2_old /mu2 * d2;
%             d3 = mu3_old /mu3 * d3;
%         end
% 
%         softfig(666);
%         semilogy(r1(1:t)); hold on;
%         semilogy(r2(1:t));
%         semilogy(r3(1:t));
%         semilogy(s1(1:t));
%         semilogy(s2(1:t));
%         semilogy(s3(1:t)); 
% 
%         hold off;
% 
%         fprintf('(r1/2 + r2/4 + r3/4)/s1 =  %g \n',(r1(t)/2 + r2(t)/4 + r3(t)/4)/s1(t));
%         fprintf('(r2/2 + r3/4 + r1/4)/s2 =  %g \n',(r2(t)/2 + r3(t)/4 + r1(t)/4)/s2(t));
%         fprintf('(r3/2 + r1/4 + r2/4)/s3 =  %g \n',(r3(t)/2 + r1(t)/4 + r2(t)/4)/s3(t));
%         fprintf('(r/s =  %g \n',sqrt(r1(t).^2+r2(t).^2+r3(t).^2)/sqrt(s1(t).^2+s2(t).^2+s3(t).^2));
%         fprintf('mu1 = %g, mu2 = %g, mu3 = %g \n',mu1,mu2,mu3);
% 
%     end
    
    if verbosity>2
        fprintf('%5i ',t);
    end
    
    if isfield(param,'display_u')
        hFig = pauseFigure(hFig);
        drawnow;
    end
    if mod(t,iter_div)==0
        rel_err(t/iter_div) = norm(v-v_prev)/norm(v);
        if verbosity>0
            fprintf('\n %5i  : \t relative error = %g \n',t,rel_err(t/iter_div));
        end
        u = v(Rn_idx); % TODO : generalize framework to handle v ?
        info.Au = Hv2;
        % Alternatively choose u = nu3
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
            if isfield(info,'Au')
                info = rmfield(info,'Au');
            end
            save(checkpoint,'t_init','v','nu1','nu2','nu3','d1','d2','d3','rel_err','info');
            
        end
        
    end
end

info.rel_err  = rel_err;
info.niter    = t;

