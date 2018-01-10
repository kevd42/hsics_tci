function hsics_main(jobid)
% HSICS_MAIN launch the reconstruction corresponding to jobid
%
% Usage:
% HSICS_MAIN
% HSICS_MAIN(jobid)
%
% jobid is an integer between 0 and 705. Use hsics_jobid to get the jobid
% corresponding to a specific setup.
%
% This can easily be launched on several machines in parallel (eg one jobid
% per machine). A matfiles checkpoint system is used so that stopped 
% simulations can be recovered and relaunched from another machine (just
% launch the job with the matfiles in the appropriate subfolders of the
% 'matfiles' folder.
%
% Note: if jobid is not provided as an input argument, the function uses
% the environment variable 'MY_JOBID' (set for example in the current 
% OS terminal).
%
% Results are stored in a matfile format.
% WARNING: Each job may generate a significant amount of data (on the order
% of 500MB). Use appropriate mass storage to avoid write errors.
%
% Author: Kevin Degraux
% Date: 2017
%
% see README.txt
%
% see also setenv, hsics_setup, hsics_jobid, hsics_settings,
% hsics_reconstruction, hsics_gather_psnr, hsics_psnr_figure,
% hsics_patches_figure
%
%  (c) UCLouvain 2018

%% JOB ID
if nargin<1
    jobid = str2double(getenv('MY_JOBID'));
end
if isempty(jobid) || isnan(jobid) || (jobid<0)
    warning('jobid could not be determined, setting default');
    
    % chart_and_stuffed_toy 5px PSF MSRC 40dB SNRin (mosaic, 1 snapshot m/n=1/16)
    jobid = 676; 
    
end

%% Simulation parameters


n         = [256,256,16];

% Synthetic simulations
n_gt = 8;       % Number of groundtruths
n_q = 4;        % Number of FPA size factor / number of snapshots
n_snr = 2;      % Number of input SNR levels
n_optical = 11; % Number of synthetic scenarios
if jobid < n_optical*n_snr*n_gt*n_q 
    
    % FPA size factor / Number of snapshots
    q   = floor(mod(jobid,n_gt*n_q)/n_gt)+1; % {1,2,3,4}
    
    % Input noise
    if mod(jobid,n_snr*n_gt*n_q) < n_gt*n_q
        SNRin = 20;
    else
        SNRin = 40;
    end

    % Ground truth
    gt_list = {'balloons','feathers','jelly_beans','glass_tiles',...
               'chart_and_stuffed_toy','stuffed_toys','superballs','beads'};
       
    gt_centers = [  255,128;
                    256,256;
                    256,256;
                    256,256;
                    230,280;
                    256,256;
                    200,236;
                    256,256  ];
    gt_id     = mod(jobid,n_gt)+1;
    gt_name   = gt_list{gt_id};
    gt_center = gt_centers(gt_id,:);
    gt_center = [max(min(gt_center(1),512-n(1)/2),n(1)/2+1),...
                 max(min(gt_center(2),512-n(2)/2),n(2)/2+1)];

    % Optical setup
    optical_list = {'InFocus_mosaic';
                    'InFocus_random';
                    'OffFocus_mosaic_multi';
                    'OffFocus_mosaic_11px_multi';
                    'OffFocus_random_multi';
                    'OffFocus_random_11px_multi';
                    'OffFocus_mosaic';
                    'OffFocus_mosaic_11px';
                    'OffFocus_random';
                    'OffFocus_random_11px';
                    'OffFocus_mosaic_5px_multi'};

    optical_id = floor(jobid/(n_snr*n_gt*n_q))+1;
    optical_str = optical_list{optical_id};
    optical_str_sep = strsplit(optical_str,'_');
    
    % Optical path is the first part of optical str
    optical_path = optical_str_sep{1};

    % fp_layout is the second part of optical str
    fp_layout = optical_str_sep{2};

    % Check for the 'multi' flag at the end
    if strcmp(optical_str_sep{end},'multi')
        m = [256,256,q^2];
    else
        m = [256*q,256*q,1];
    end

    % Check for the psf size in 3rd position
    withPSF = false;
    psf_rad_max = 5.7514; % Default 55e-6 m fpa pixel pitch (bins of 10 by 10 pixels)
    if length(optical_str_sep)>=3
        if strcmp(optical_str_sep{3},'11px')
            withPSF     = true;
        elseif strcmp(optical_str_sep{3},'5px')
            withPSF     = true;
            psf_rad_max = 2.8757; % 110e-6 m fpa pixel pitch (bins of 20 by 20 5.5 micron pixels)
        end
    end
    
% Experimental data (2 cases)
elseif jobid < n_optical*n_snr*n_gt*n_q + 1 
    optical_path = 'InFocus';
    m            = [512,512,1];
    fp_layout    = 'mosaic';
    withPSF      = false;
    gt_name      = 'leaves';
    SNRin        = 40;
elseif jobid < 11*n_snr*n_gt*n_q + 2
    optical_path = 'InFocus';
    m            = [512,512,1];
    fp_layout    = 'mosaic';
    withPSF      = false;
    gt_name      = 'RedYellow';
    SNRin        = 40;
end


%%
% ==== PIHT ====
% K. Degraux, V. Cambareri, L. Jacques, B. Geelen, C. Blanch and 
% G. Lafruit, "Generalized inpainting method for hyperspectral image 
% acquisition," 2015 IEEE International Conference on Image Processing 
% (ICIP), Quebec City, QC, 2015, pp. 315-319.
% doi: 10.1109/ICIP.2015.7350811

% algo.name = 'PIHT'
% % 20 dB
% % K = [0.1,0.5,1,2]*2;
% % K = [0.01,0.05,1,2]*2;
% 
% % Noiseless
% algo.K = [0.1,0.5,5,10]*2; % "normal K"
% %K = [0.01,0.05,5,10]*2; % "K slow" gives much better results with InFocus
% % K = [0.01,0.05,0.1,0.5,1,5,10]*2; % same as previous.
% 
% % K must depend on SNRin


% ==== ADMM ====
algo.name = 'ADMM';

% OFF FOCUS
algo.rho  = 40; % OK for OffFocus (and InFocus)

% IN FOCUS
%algo.rho = 10; % Works OK for InFocus
% algo.rho = 5; % More oscillations but seemingly quicker convergence
% algo.rho = 8; % Good tradeoff
% algo.rho = 1; % Too much oscilations

settings = hsics_settings(gt_name,gt_center,fp_layout,...
                                optical_path,withPSF,psf_rad_max,SNRin,n,m,algo);

settings.jobid = jobid;

%% LAUNCH JOB
fprintf('\n\nStarting reconstruction...\n\n')
hsics_reconstruction(settings);

fprintf('Job %i terminated sucessfully\n',jobid);
 
end