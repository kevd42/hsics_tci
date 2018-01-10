function settings = hsics_settings(gt_name,gt_center,fp_layout,...
                                       optical_path,withPSF,psf_rad_max,SNRin,n,m,algo)
% HSICS_SETTINGS
% Generates settings for the hsics_reconstruction
% Usage:
% HSICS_SETTINGS(gt_name,gt_center,fp_layout,...
%                optical_path,withPSF,psf_rad_max,SNRin,n,m,algo)
%
% Inputs:
% gt_name:      'balloons', 'feathers', 'jelly_beans', 'glass_tiles',
%               'chart_and_stuffed_toy', 'stuffed_toys', 'superballs', 
%               'beads' (synthetic),'leaves','RedYellow' (experimental)
% gt_center:    center of the ROI for cropping the groundtruth
% fp_layout:    'mosaic', 'random'
% optical_path: 'InFocus' (MSVI), 'OffFocus' (MSRC)
% withPSF:      false, true (for OffFocus only)
% psf_rad_max:  2.8757 (5 px PSF), 5.7514 (11 px PSF) or any float
% SNRin:        input SNR in dB. Used for generating the noise (synthetic
%               experiments only) and setting the regularization parameter.
% n:            size of the target (array of int of length 3)
% m:            size of the measurements and nb of snapshots (array of int 
%               of length 3)
% 
% Output:
% settings: struct containing the settings (operators, functions,...) for
%           the hsics_reconstruction function.
%
% see README.txt
%
% date : 06/04/2017
% author : K. Degraux
%
%  (c) UCLouvain 2018

fprintf('Initialization...\n')

ni = n(1);
nj = n(2);
ns = n(1:2);
nl = n(3);

mi = m(1);
mj = m(2);
ms = m(1:2);
mS = m(3);


n_str = sprintf('%ix%ix%i',n(1),n(2),n(3));
m_str = sprintf('%ix%ix%i',m(1),m(2),m(3));


% fprintf('Dimensions: n=%s, m=%s\n',n_str,m_str);

%% FPA 
fprintf('FPA: %s\n',fp_layout);
% Seed the random number generator
rng(sum(m)+sum(n)+double(withPSF)+sum(double(fp_layout)));
sensing.rngstate = rng;

if strcmp(fp_layout,'mosaic')
    fp_mosaic = [3   4   2   1;
                 11  12  10  9;
                 15  16  14  13;
                 7   8   6   5];
    idL  = repmat(fp_mosaic,ms./[4,4]);
else
    idL  = repmat(reshape(1:nl, 4, 4)', ms./[4,4]);
    idL  = reshape(idL(randperm(numel(idL))),size(idL));
end
idLy = arrayfun(@(j) (idL == j),1:nl,'UniformOutput',false);



%% Box constraint

settings.upperbnd = 1;
settings.lowerbnd = 0;

%% General

% fprintf('Algorithm: %s\n',algo.name);

settings.identifier = sprintf('HSI_%s_%s',algo.name,optical_path);

if strcmp(algo.name,'PIHT')
    algstr = sprintf('%.0f_',algo.K);
elseif strcmp(algo.name,'ADMM')
    algstr = sprintf('%.0f_',algo.rho);
end

settings.filename   = sprintf('%s_%s%s_%s',settings.identifier,algstr,n_str,m_str);
settings.n          = n;
settings.m          = m; 

%% Prior


prior.name         = 'UDWT_DCT2_kron_DCT'; 
prior.alias        = sprintf('UWDkD_%s',n_str);
prior.generateA    = @() constructTightUDWT_DCT2_kron_DCT(n);
prior.compNorm     = false;
prior.generatePsi0 = @() deal(1,[]); % dummy 1 input, 2 outputs.
prior.generatePsi  = @() deal(1,[]);

% fprintf('Prior: %s \n',prior.name );

%% Groundtruth
groundtruth.name       = gt_name;
groundtruth.alias      = sprintf('%s_%s',gt_name,n_str);
groundtruth.graphics   = @(x,u,A) displayTargetVolume(x,n,1);


fprintf('Groundtruth: %s \n',groundtruth.alias );

%% Sensing

fprintf('Sensing: %s \n',optical_path );

if strcmp(optical_path,'OffFocus')
    
    % === SLM ===
    s = ns + ms - 1;
    slm = sign(randn([s,mS]));
    %slm = slm/prod(s);
    if withPSF
        
        fprintf('Discretizing the PSF...');
        psf = diffPSF(ms,ns,linspace(470,620,16)*1e-9,psf_rad_max);
        for band = 1:nl
            psf(:,:,band) = psf(:,:,band)/sum(sum(abs(psf(:,:,band))));
        end
        fprintf('done\n')
        fprintf('Computing the diffracted Coded Aperture patterns...');
        ss = [s(1)+ceil((s(1)-1)/2),s(2)+ceil((s(2)-1)/2)]; % Size of the conv
        SLM = zeros(ss(1),ss(2),mS);
        for snap = 1:mS
            SLM(:,:,snap) = fft2(slm(:,:,snap),ss(1),ss(2));
        end
        PSF = zeros(ss(1),ss(2),nl);
        for band = 1:nl
            PSF(:,:,band) = fft2(psf(:,:,band),ss(1),ss(2));
        end
        
        slm_psf = zeros(ss(1),ss(2),nl);
        for snap = 1:mS
            for band = 1:nl
                slm_psf(:,:,snap,band) = ifft2(PSF(:,:,band).*SLM(:,:,snap),'symmetric');
            end
        end
        slm = slm_psf(end-s(1)+1:end,end-s(2)+1:end,:,:);
        
        psf_str = sprintf('PSF%1.4f',psf_rad_max);
        fprintf('done\n')
    else
        psf_str = 'noPSF';
    end
    
    sensing.generatePhi = @() constructOffFocus(m,n,slm,idLy,1,1);
else
    psf_str = 'noPSF';
    sensing.generatePhi = @() constructInFocus(m,n,idLy,1);
end
sensing.name        = sprintf('%s_%s_%s_%d',optical_path,psf_str,...
                                          fp_layout,sensing.rngstate.Seed);
sensing.alias       = sprintf('%s_%s_%s_%d',optical_path,psf_str,...
                                          fp_layout,sensing.rngstate.Seed);
sensing.renew       = false;

sensing.compMeanCol = false;
sensing.compNorm    = true;
sensing.synthetic   = true;

%% Measurements
if strcmp(gt_name,'leaves') || strcmp(gt_name,'RedYellow')
    
    groundtruth.available  = false;
    measurements.synthetic    = false;
    measurements.SNRin.true   = SNRin;
    measurements.SNRin.target = SNRin;
    
    switch gt_name
        case 'leaves';    roi = [1,128,195,128];
        case 'RedYellow'; roi = [35,128,201,128];
    end
    measurements.loadMes = @() importMosaicExample(gt_name,roi,idL,m);
else
    groundtruth.available  = true;
    groundtruth.generateGT = @(A,Psi0) deal(importCaveGroundtruth(gt_name,gt_center,n),[],[]);
    
    measurements.synthetic    = true;
    measurements.SNRin.true   = SNRin;
    measurements.SNRin.target = SNRin;
end


measurements.computeU0    = @(Phi,info,y,epsilon)...
                        computeInitFP(info,y,epsilon,idL,optical_path,m,n);
measurements.graphicsU0   = @(x0,u0,y,A) displayInitVolume(x0,n,2);


%% Algorithm
if strcmp(algo.name,'PIHT')
    algorithm.recontruct = @(settings) launch_lys_syn_pos_PIHT(settings,algo.K);
elseif strcmp(algo.name,'ADMM')
    algorithm.recontruct = @(settings) launch_ABPDN_box_ADMM(settings,algo.rho);
end


%%

settings.prior        = prior;
settings.groundtruth  = groundtruth;
settings.sensing      = sensing;
settings.measurements = measurements;
settings.algorithm    = algorithm; 