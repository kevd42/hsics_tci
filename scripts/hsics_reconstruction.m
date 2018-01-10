function hsics_reconstruction(settings)
% HSICS_RECONSTRUCTION
%   HSICS_RECONSTRUCTION(settings)
%   Solve inverse problems of the form
%
%   x = Psi0 * argmin_u G( A * u )   +  F( y - Phi * u ) + C( Psi * u )
%
%   Where G is some "complexity function" (typically l1-norm or sparsity
%   constraint), F is a "Fidelity function" (typically, l2-tube
%   constraint or "sum of squares" terms), C is some additional constraint
%   (typically a box or positivity constraint).
%   The solution quality is evaluated after transformation of u with Psi0.
%   In general Psi needs not be equal to Psi0.
%   
%   This function is meant as a framework that takes care of saving/loading
%   files and calling the functions that do the actual job.
% 
%   settings is a nested structure containing the following fields:
%       identifier
%       jobid
%       filename
%       n         [ni,nj,nl] size of x, the HS volume of interest.
%       m         [mi,mj,mS] size of the stack of snapshots.
%       prior
%           name               String identifying the prior.
%           alias              Short name for file names.
%           @() generateA      Function handle generating the Analysis
%                              dictionary. Should return A or [A,A_info].
%                              A_info is an arbitrary structure containing
%                              additional info that will be saved in the
%                              prior matfile.
%           compNorm           Boolean that indicates if we need to
%                              estimate the operator norm of A.
%           @() generatePsi0   Function handle generating the synthesis
%                              dictionary. Should return Psi0 or 
%                              Psi0,Psi0_info].
%                              Psi0 is used only to generate x from u i.e.
%                              for PSNR computation (if GT is available) and
%                              for graphics. Note that in general, Psi0 must
%                              not necessarily be identical to Psi.
%                              Phi0_info is similar to A_info.
%           @() generatePsi    Function handle generating the synthesis
%                              dictionary. Should return Psi or 
%                              [Psi,Psi_info].
%                              Psi is used to generate a proxy of x from u 
%                              for the constraint C.
%                              Psi_info is similar to A_info.
%       groundtruth
%           name               String identifying the grountruth (GT).
%           alias              Short name for file names.
%           available          Boolean that indicates if GT is available.
%           @generateGT(A,Psi0)
%                              Function handle that generates the
%                              groundtruth. Output [u], [x,u] or 
%                              [x,u,GT_info]. In principle x = Psi0 * u.
%           @graphics(x,u,A)   Function handle that generates graphics from
%                              the groundtruth. No output.
%           dyn                (optional) dynamic of the groundtruth
%                              entries (default = 1)
%       sensing
%           name               String identifying the sensing matrix.
%           alias              Short name for file names.
%           renew              Boolean that indicates if we must ignore the
%                              previously saved files.
%           @generatePhi()     Function handle generating the sensing
%                              operator mapping y = Phi * u. Be careful to
%                              INCLUDE THE SYNTHESIS (Psi) in the sensing
%                              operator. Output Phi or [Phi,Phi_info]. 
%                              Phi must support multiplication, 
%                              adjoint and left-division (e.g. Spot 
%                              operators do the job). As usual, Phi_info 
%                              is general purpose arbitrary saved info.
%           compMeanCol        Boolean that indicates if we need to
%                              estimate the average column norm.
%           compNorm           Boolean that indicates if we need to
%                              estimate the operator norm.
%       measurements
%           synthetic          Boolean that indicates if the measurements
%                              are synthetic or come from real data.
%           SNRin
%               true           True PSNR of the noise added to synthetic
%                              measurements
%               target         Target PSNR if needed by the algorithm.
%           @loadMes()         Function handle importing the measurements.
%                              Output y.
%           @computeU0(Phi,normPhi,y,epsilon)
%                              Function handle computing the initial point
%                              for the algorithm, knowing Phi, the
%                              measurement vector y and epsilon. Output
%                              u0.
%           @graphicsU0(x0,u0,A)
%                              Function handle that generates graphics from
%                              the initial point. No output.
%                              
%       algorithm
%           @recontruct(settings)
%                              Function handle to the reconstruction
%                              algorithm. Output [u_est,alg_info].
% see README.txt
%
% Author : K. Degraux (ISPGroup, ICTEAM, UCLouvain)
% Date : 30 July 2015
%
%  (c) UCLouvain 2018



if ~isdeployed
    graphics = true;
else
    graphics = false;
end

%% Settings nested structures unfold for code readability
% When editing a child structure, don't forget to update the parent

sensing      = settings.sensing;
prior        = settings.prior;
groundtruth  = settings.groundtruth;
measurements = settings.measurements;
SNRin        = settings.measurements.SNRin;
identifier   = settings.identifier;
jobid        = settings.jobid;
filename     = [settings.filename,'_j',num2str(jobid)];
settings.filename = filename;

% Target resolution
n  = settings.n;
ni = n(1);      % vertical
nj = n(2);      % horizontal
nl = n(3);      % wavelength
ns = [ni,nj];   % spatial
n  = [ns,nl];   % HS volume resolution
N  = prod(n);   % total nb of voxels

% FPA resolution
m  = settings.m; % measurement "stack" resolution
mi = m(1); % vertical
mj = m(2); % horizontal
mS = m(3); % Number of snapshots
ms = [mi,mj];
Ms = prod(ms);
M  = prod(m); % Total nb of measurements   


fprintf('Target resolution : %ix%ix%i (N = %i)\n',ni,nj,nl,N);
fprintf('FPA resolution : %ix%i (Ms = %i = %1.2f N)\n',mi,mj,Ms,Ms/N);
fprintf('Number of snapshots : %i (M = %i = %0.2f)\n',mS,M,M/N);




%% Define priors
fprintf('Priors : %s\n',prior.name);

if ~exist('matfiles','dir')
    mkdir matfiles
end

prior.file = sprintf('matfiles/prior/%s_%ix%ix%i.mat',prior.alias,n);


if exist(prior.file,'file')
    prior = load(prior.file);
    A    = prior.A;
    Psi0 = prior.Psi0;
    Psi  = prior.Psi; %#ok<NASGU>
    fprintf('Prior file %s loaded\n',prior.file);
else
    if ~exist('matfiles/prior','dir')
        mkdir('matfiles/prior');
    end
    
    % Handle to the generating function of A

    [A,prior.A_info] = prior.generateA();         

    
    % Handle to the generating function of Psi0
    [Psi0,prior.Psi0_info] = prior.generatePsi0(); 

    
    % Handle to the generating function of Psi
    [Psi,prior.Psi_info] = prior.generatePsi(); 

    
    prior = rmfield(prior,{'generateA','generatePsi0','generatePsi'});
    
    prior.A      = A;
    prior.Psi0   = Psi0;
    prior.Psi    = Psi;

    
    % Compute norm of A
    if prior.compNorm
        prior.A_info.norm   = normest(prior.A,1e-3)+1e-3;
        fprintf('The norm of the analysis prior operator A is %1.2f\n',prior.A_info.norm);
    end
    

    % Save prior
    save(prior.file,'-struct','prior');
    fprintf('Prior saved in %s\n',prior.file);
end
settings.prior = prior;



%%            
if groundtruth.available
    %% Import and pre-process groundtruth
    groundtruth.file = sprintf('matfiles/groundtruth/%s_%ix%ix%i.mat',groundtruth.alias,n);

    if exist(groundtruth.file,'file')
        groundtruth = load(groundtruth.file);
        x = groundtruth.x;
        u = groundtruth.u;
        fprintf('Groundtruth file %s loaded\n',groundtruth.file)
    else
        if ~exist('matfiles/groundtruth','dir')
            mkdir('matfiles/groundtruth');
        end
        
        % Handle to the generating function of x and u
        [x,u,groundtruth.GT_info] = groundtruth.generateGT(A,Psi0); 
        if isempty(u)
            u = Psi0'*x;
            x = Psi0*u;
        end
        
        groundtruth.x = x;
        groundtruth.u = u;
        
        groundtruth = rmfield(groundtruth,'generateGT');
        
        save(groundtruth.file,'-struct','groundtruth')
        fprintf('Groundtruth saved in %s\n',groundtruth.file)
    end

    settings.groundtruth = groundtruth;


    %% Graphics (Groundtruth)
    if graphics
        groundtruth.graphics(x,u,A) % Handle to the graphics function for the groundtruth
        drawnow;
        pause(0.01);
    end
end

%% Construct the sensing operator

fprintf('Sensing operator : %s \n',sensing.name);

sensing.file = sprintf('matfiles/sensing/%s_%ix%ix%i_%ix%ix%i',...
    sensing.alias,...
    n,...
    m);

if sensing.renew
    sensing.file = sprintf('%s_j%i.mat',sensing.file,jobid);
else
    sensing.file = [sensing.file,'.mat'];
end

load_sensing = exist(sensing.file,'file') && ~sensing.renew;
if load_sensing
    try
        sensing = load(sensing.file);
        Phi     = sensing.Phi;
    %     if isfield(sensing,'Phi_info') && isfield(sensing.Phi_info,'normPhi') 
    %         normPhi = sensing.Phi_info.normPhi;
    %     else
    %         normPhi = 1;
    %     end
        fprintf('Sensing file %s loaded\n',sensing.file)
    catch
        warning('There was a problem during load');
        load_sensing = false;
        if exist(sensing.file,'file')
            fprintf('Removing corrupt files')
            delete(sensing.file);
        end
    end
end
if ~load_sensing
    if ~exist('matfiles/sensing','dir')
        mkdir('matfiles/sensing');
    end
    
    %sensing.rngstate = rng(jobid);
    
    % Handle to the generating function of Phi
    % Reminder: Phi = Phi0 * Psi (see help) 

    [Phi,sensing.Phi_info] = sensing.generatePhi();       

    sensing = rmfield(sensing,'generatePhi');
    sensing.Phi = Phi;
    

    % Monte Carlo estimation of the average column norm
    if sensing.compMeanCol
        % USEFULL?
        fprintf('Average column norm\n');
        sensing.Phi_info.meanCol= MCMeanColNorm(sensing.Phi);
        fprintf('Sensing operator average column norm is %1.5f\n',sensing.meanColPhi);
    end
    
    % Estimation of the operator norm
    sensing.Phi_info.normPhi = 1;
    if sensing.compNorm
        fprintf('Computing sensing operator norm...\n');
        sensing.Phi_info.normPhi  = normest(sensing.Phi ,1e-3);%+1e-3;
        normPhi = sensing.Phi_info.normPhi;
        fprintf('The norm of the sensing operator Phi is %g\n', normPhi);
    end
    try
        save(sensing.file,'-struct','sensing')
        fprintf('Sensing saved in %s\n',sensing.file)
    catch
        warning('There was a problem when saving sensing file');
        if exist(sensing.file,'file')
            fprintf('Removing corrupt files\n')
            delete(sensing.file);
        end
    end
end
settings.sensing = sensing;



%% Quality criterion function handles

SNR   = @(x,x_est) 20*log10( norm(x(:))/norm(x_est(:) - x(:)) );
NSNR  = @(x,x_est) 20*log10( 1/norm(x_est(:)/norm(x_est(:)) - x(:)/norm(x(:))) ); %#ok<NASGU>
MSE   = @(x,x_est) sum((x(:)-x_est(:)).^2)/numel(x);
if isfield(settings.groundtruth, 'dyn')
    dyn = settings.groundtruth.dyn;
else
    dyn = 1;
end
PSNR  = @(x,x_est) 20*log10(dyn) - 10*log10(MSE(x(:),x_est(:)));

%% Measurements, noise estimation and initial guess u0
settings.checkpointfile = sprintf('matfiles/checkpoint/%s/init_%s.mat',identifier,filename);
if ~exist(settings.checkpointfile,'file')
    if ~exist('matfiles/checkpoint','dir')
        mkdir matfiles/checkpoint;
    end
    
    if ~exist(['matfiles/checkpoint/',identifier],'dir')
        mkdir(['matfiles/checkpoint/',identifier]);
    end
    
    %% Measurements and noise estimation
    
    
    fprintf('Measurements and noise estimation\n');
    
    if measurements.synthetic

        y0 = Phi*u;

        varnoise = norm(y0)^2 / 10^(SNRin.true/10);
        
        noise = randn(M,1)*sqrt(varnoise)/sqrt(M);
        
        y = y0 + noise;
        
        if SNRin.target == SNRin.true
            epsilon = norm(noise)*1.01; % Noise oracle with 1% tolerence
        else
            epsilon = norm(y) / 10^(SNRin.target/20);
        end
        
    else
        
        y = measurements.loadMes();
        
        epsilon = norm(y) / 10^(SNRin.target/20);
        
    end
    
    measurements.y       = y;
    measurements.epsilon = epsilon;

    fprintf('Target max distance to y : epsilon = %1.3g (%2.2f dB)\n',epsilon,SNRin.target);
    
    %% Initial guess u0
    fprintf('Initial guess u0\n');
    [u0,u0_info] = measurements.computeU0(Phi,sensing.Phi_info,y,epsilon);

    measurements.u0      = u0;
    measurements.u0_info = u0_info;
    settings.measurements = measurements;
    
    %% CHECKPOINT HERE

    save(settings.checkpointfile,'-struct','measurements');

else
    measurements = load(settings.checkpointfile);
    fprintf('Initial checkpoint %s loaded\n',settings.checkpointfile);
    settings.measurements = measurements;
    y       = measurements.y;
    u0      = measurements.u0;
    epsilon = measurements.epsilon; %#ok<NASGU>
end

x0 = Psi0*u0;

%% Graphics : initial guess
if groundtruth.available
    fprintf('PSNR0 = %2.1f dB\n',PSNR(x,x0));
end
    
if graphics
    measurements.graphicsU0(x0,u0,y,A)
    
    drawnow;
    pause(0.01);
end

%% Reconstruction algorithm

if ~exist('matfiles/results','dir')
    mkdir('matfiles/results');
end
if ~exist(['matfiles/results/',identifier],'dir')
    mkdir(['matfiles/results/',identifier]);
end
settings.resultfile = sprintf('matfiles/results/%s/res_%s.mat',identifier,filename);

if ~exist(settings.resultfile,'file')
    tstart = tic;
    [u_est,alg_info] = settings.algorithm.recontruct(settings);  %#ok<ASGLU>
    time_algo = toc(tstart);                                    %#ok<NASGU>
    
    prior_file       = prior.file;                              %#ok<NASGU>
    if groundtruth.available
        groundtruth_file  = groundtruth.file;                         %#ok<NASGU>
    else
        groundtruth_file = []; %#ok<NASGU>
    end
    sensing_file      = sensing.file;                           %#ok<NASGU>                       %#ok<NASGU>
    measurements_file = settings.checkpointfile;                %#ok<NASGU>
    save(settings.resultfile,'u_est','alg_info','prior_file','groundtruth_file','sensing_file','measurements_file');

else
    load(settings.resultfile);
end

%% Compute Final PSNR/SNR

if groundtruth.available
    x_est = Psi0 * u_est;
    Xten_est = reshape(x_est,n);
    Xten = reshape(x,n);
    for band= 1:nl
        fprintf('band %i , PSNR = %3.2f (SNR = %3.2f) \n',...
            band,PSNR(Xten(:,:,band),Xten_est(:,:,band)),SNR(Xten(:,:,band),Xten_est(:,:,band)));
    end
end


end