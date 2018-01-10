function hsics_gather_psnr(results_folder,optical_setup)
%% HSICS_GATHER_PSNR
% Load 'res' files and gather the psnr and spsnr maps in 11x2 files (2 per
% optical configuration)
% 
% The list of optical configurations is
%
% optical_list = {'InFocus_mosaic';
%                 'InFocus_random';
%                 'OffFocus_mosaic_multi';
%                 'OffFocus_mosaic_11px_multi';
%                 'OffFocus_random_multi';
%                 'OffFocus_random_11px_multi';
%                 'OffFocus_mosaic';
%                 'OffFocus_mosaic_11px';
%                 'OffFocus_random';
%                 'OffFocus_random_11px';
%                 'OffFocus_mosaic_5px_multi'};
%
% Usage:
% HSICS_GATHER_PSNR
% HSICS_GATHER_PSNR(results_folder)
% HSICS_GATHER_PSNR(results_folder,optical_setup)
%
% If optical_setup is not specified (or empty), the script will loop
% through all of them.
%
% The results_folder is the 'matfiles' folder where the results are stored.
% The files will be saved in that folder.
%
% see README.txt
%
% Author: K. Degraux
% Date: 8/01/2018
%
%  (c) UCLouvain 2018


%results_folder   = '/dir/altersense/kdegraux/matfiles';

if nargin<1
    results_folder='matfiles';
end
if nargin<2
    optical_setup=[];
end

n = [256,256,16];

n_gt      = 8;
n_q       = 4;
n_snr     = 2;
%n_layout  = 2;
n_optical = 11;

gt_list = {'balloons','feathers','jelly_beans','glass_tiles',...
           'chart_and_stuffed_toy','stuffed_toys','superballs','beads'};

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

for jobid = 0 : n_gt*n_q*n_snr*n_optical-1
    
    if mod(jobid,n_snr*n_gt*n_q)==0

        psnr = zeros(n_gt,n_q,n_snr);

        spsnr_maps = zeros(n(1),n(2),n_gt,n_q,n_snr);

    end

    disp(jobid);
    %% Settings
    n_gt = 8;
    gt_id = mod(jobid,n_gt)+1;
    gt_name = gt_list{gt_id};

    n_q = 4;
    q   = floor(mod(jobid,n_gt*n_q)/n_gt)+1; % {1,2,3,4}

    n_snr = 2;
    if mod(jobid,n_snr*n_gt*n_q) < n_gt*n_q
        SNRin = 20;
    else
        SNRin = 40;
    end
    
    optical_id = floor(jobid/(n_snr*n_gt*n_q))+1;
    optical_str = optical_list{optical_id};
    
    % Optical path is the first part of optical str
    optical_str_sep = strsplit(optical_str,'_');
    optical_path = optical_str_sep{1};
    
    % Check for the 'multi' flag
    if strcmp(optical_str_sep{end},'multi')
        m = [256,256,q^2];
    else
        m = [256*q,256*q,1];
    end

    n         = [256,256,16];

    algo.name = 'ADMM';
    algo.rho  = 40;

    % settings.jobid = jobid; 
    n_str = sprintf('%ix%ix%i',n(1),n(2),n(3));
    m_str = sprintf('%ix%ix%i',m(1),m(2),m(3));

    settings.identifier = sprintf('HSI_%s_%s',algo.name,optical_path);

    if strcmp(algo.name,'PIHT')
        algstr = sprintf('%.0f_',algo.K);
    elseif strcmp(algo.name,'ADMM')
        algstr = sprintf('%.0f_',algo.rho);
    end

    settings.filename = sprintf('%s_%s%s_%s',settings.identifier,algstr,n_str,m_str);

    settings.groundtruth.name       = gt_name;
    settings.groundtruth.alias      = sprintf('%s_%s',gt_name,n_str);

    %% Loading res files
    % checkpoint_file  = [results_folder,'/checkpoint/',settings.identifier,...
    %                     '/cp_',settings.filename,'_j',num2str(jobid),'.mat'];
    % init_file        = [results_folder,'/checkpoint/',settings.identifier,...
    %                     '/init_',settings.filename,'_j',num2str(jobid),'.mat'];
    if isempty(optical_setup) || strcmp(optical_setup,optical_str)
        try
            res_file         = [results_folder,'/results/',settings.identifier,...
                                '/res_',settings.filename,'_j',num2str(jobid),'.mat'];

            load(res_file);

            groundtruth_file = [results_folder,'/groundtruth/',...
                                settings.groundtruth.alias,sprintf('_%dx%dx%d.mat',n)];

            load(groundtruth_file);
            snr_id = SNRin/20;

            psnr(gt_id,q,snr_id) = -10*log10( sum((u-u_est).^2)./length(u) );

            spsnr_maps(:,:,gt_id,q,snr_id) = ...
                                   -10*log10( sum(reshape(u-u_est,n).^2,3)./n(3) );

        catch
            warning('jobid %i, %s not available',jobid,res_file);
        end

        %% Saving gathered psnr and spsnr maps
        if mod(jobid,n_snr*n_gt*n_q)==n_snr*n_gt*n_q-1
            fprintf('Saving %s\n',optical_str)
            save([results_folder,'/',optical_str,'_psnr.mat'],'psnr');
            save([results_folder,'/',optical_str,'_spsnr_maps.mat'],'spsnr_maps');
        end
    end

end

end


