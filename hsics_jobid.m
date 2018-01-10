function jobid_out = hsics_jobid(gt_name_in,...
                             q_in,...
                             optical_path_in,...
                             fp_layout_in,...
                             multi_snap_in,...
                             withPSF_in,...
                             psf_rad_max_in,...
                             SNRin_in)
% HSICS_JOBID display name of files in which the results and checkpoints 
% are stored and return job id corresponding to a specific simulation
%
% USAGE: 
% jobid = HSICS_JOBID(gt_name_in,...
%                     q_in,...
%                     optical_path_in,...
%                     fp_layout_in,...
%                     multi_snap_in,...
%                     withPSF_in,...
%                     psf_rad_max_in,...
%                     SNRin_in)
%
% diplays all jobs corresponding to the matching pater if an argument is
% the empty array []
% 
% Accepted values (except [])
% gt_name_in:      'balloons', 'feathers', 'jelly_beans', 'glass_tiles',
%                  'chart_and_stuffed_toy', 'stuffed_toys', 'superballs', 
%                  'beads' (synthetic),'leaves','RedYellow' (experimental)
% q_in:            1,2,3,4 (factor for the subsampling rate so that
%                  m/n = q^2/16)
% optical_path_in: 'InFocus' (MSVI), 'OffFocus' (MSRC)
% fp_layout_in:    'mosaic', 'random'
% multi_snap_in:   false, true (for OffFocus only)
% withPSF_in:      false, true (for OffFocus only)
% psf_rad_max_in:  2.8757 (5 px PSF), 5.7514 (11 px PSF) (for OffFocus 
%                  only)
% SNRin_in:        20, 40
%
% Total number of accepted combinations: 706
%
% returns the LAST acceptable jobid (between 0 and 705)
%
% see README.txt
%
% Author: Kevin Degraux
% Date: 2017
%
% see also hsics_main, hsics_setup
%  (c) UCLouvain 2018

                       
                       
for jobid = 0:705


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
        multi_snap = true;
    else
        m = [256*q,256*q,1];
        multi_snap = false;
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

n         = [256,256,16];


% ==== ADMM ====
algo.name = 'ADMM';
algo.rho  = 40; % OK for OffFocus (and InFocus)


% settings = HSI_journal_settings(gt_name,gt_center,fp_layout,...
%                                 optical_path,withPSF,psf_rad_max,SNRin,n,m,algo);
% 
% settings.jobid = jobid;

gt_name_out = gt_name_in;
if isempty(gt_name_in); gt_name_in = gt_name; end;

q_out = q_in;
if isempty(q_in)
    q_in = q; 
end

optical_path_out = optical_path_in;
if isempty(optical_path_in)
    optical_path_in = optical_path;
end

fp_layout_out = fp_layout_in;
if isempty(fp_layout_in)
    fp_layout_in = fp_layout;
end

multi_snap_out = multi_snap_in;
if isempty(multi_snap_in)
    multi_snap_in = multi_snap;
end

withPSF_out = withPSF_in;
if isempty(withPSF_in)
    withPSF_in = withPSF;
end

psf_rad_max_out = psf_rad_max_in;
if isempty(psf_rad_max_in)
    psf_rad_max_in = psf_rad_max;
end

SNRin_out = SNRin_in;
if isempty(SNRin_in)
    SNRin_in = SNRin; 
end


if strcmp(gt_name,gt_name_in) &&...
   q == q_in &&...
   strcmp(optical_path,optical_path_in) &&...
   strcmp(fp_layout,fp_layout_in) &&...
   multi_snap == multi_snap_in &&...
   withPSF == withPSF_in &&...
   psf_rad_max == psf_rad_max_in &&...
   SNRin == SNRin_in
    fprintf('jobid=%i  ',jobid);
    fprintf('q=%i, ',q); 
    fprintf('optical_path=%s, ',optical_path); 
    fprintf('fp_layout=%s, ',fp_layout);  
    fprintf('multi_snap=%i, ',double(multi_snap));   
    fprintf('withPSF=%i, ',double(withPSF));   
    fprintf('psf_rad_max=%f, ',psf_rad_max);  
    fprintf('SNRin=%f\n',SNRin);    
    
    jobid_out = jobid;
end


gt_name_in      = gt_name_out;
q_in            = q_out;
optical_path_in = optical_path_out;
fp_layout_in    = fp_layout_out;
multi_snap_in   = multi_snap_out;
withPSF_in      = withPSF_out;
psf_rad_max_in  = psf_rad_max_out;
SNRin_in        = SNRin_out;

end