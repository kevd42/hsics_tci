function hsics_patches_figure(results_folder, gt_name,optical_list,q,...
                              SNRin,inits, with_gt, rois)
% HSICS_PATCHES_FIGURE 
% Create a figure with several regions of interest (rois) to compare
% several optical setups on one groundtruth/SNRin/q configuration.
% ONLY FOR SYNTHETIC EXAMPLES
% 
% Usage:
% HSICS_PATCHES_FIGURE
% HSICS_PATCHES_FIGURE(results_folder)
% HSICS_PATCHES_FIGURE(results_folder, gt_name)
% HSICS_PATCHES_FIGURE(results_folder, gt_name, optical_list)
% HSICS_PATCHES_FIGURE(results_folder, gt_name, optical_list, q)
% HSICS_PATCHES_FIGURE(results_folder, gt_name, optical_list, q, SNRin)
% HSICS_PATCHES_FIGURE(results_folder, gt_name, optical_list, q, SNRin,...
%                      inits)
% HSICS_PATCHES_FIGURE(results_folder, gt_name, optical_list, q, SNRin,...
%                      inits, with_gt)
% HSICS_PATCHES_FIGURE(results_folder, gt_name, optical_list, q, SNRin,...
%                      inits, with_gt, rois)
%
% Inputs:
% results_folder: folder containing all the matfiles (default 'matfiles')
% gt_name: name of the groundtruth from the following list:
%        gt_list = {'balloons','feathers','jelly_beans','glass_tiles',...
%            'chart_and_stuffed_toy','stuffed_toys','superballs','beads'};
%          default: 'chart_and_stuffed_toy'
% optical_list: list of all the optical setups that must appear on the
%               figure, in the order in which they appear. The list of
%               acceptable setupts is:
%                    {'InFocus_mosaic';
%                     'InFocus_random';
%                     'OffFocus_mosaic_multi';
%                     'OffFocus_mosaic_11px_multi';
%                     'OffFocus_random_multi';
%                     'OffFocus_random_11px_multi';
%                     'OffFocus_mosaic';
%                     'OffFocus_mosaic_11px';
%                     'OffFocus_random';
%                     'OffFocus_random_11px';
%                     'OffFocus_mosaic_5px_multi'};
%     default is 
%     optical_list = {'InFocus_mosaic';
%                     'InFocus_random';
%                     'InFocus_mosaic';
%                     'InFocus_random';
%                     'OffFocus_mosaic_multi';
%                     'OffFocus_mosaic_5px_multi';
%                     'OffFocus_mosaic_11px_multi'};
% q: 1,2,3,4 (factor for the subsampling rate so that m/n = q^2/16)
% SNRin: 20 or 40 (default) input SNR
% inits: array (at least) of the same size as optical_list. When
%       inits(i)==1, the ith row of the plot will display of initialization 
%       result of the optical_list{i} setup. 
%       Default: inits = [1,1,0,0,0,0,0];
% with_gt: if true (default) the groundtruth is used for the first row and
%          spectral errors are shown.
% rois: n_cols by 3 array such that [rois(:,1), rois(:,2)] are the upper 
%       left corners and rois(:,3) are the spectral band indices of each
%       region of interest to be displayed 
%
% see README.txt
%
% Author: K. Degraux
% Date: 8/01/2018
%
%  (c) UCLouvain 2018


if nargin<1
    results_folder = 'matfiles';
end
gt_list = {'balloons','feathers','jelly_beans','glass_tiles',...
           'chart_and_stuffed_toy','stuffed_toys','superballs','beads'};
       
if nargin<2
    gt_name = gt_list{5};
end
if nargin<3
    optical_list = {'InFocus_mosaic';
                    'InFocus_random';
                    'InFocus_mosaic';
                    'InFocus_random';
                    'OffFocus_mosaic_multi';
                    'OffFocus_mosaic_5px_multi';
                    'OffFocus_mosaic_11px_multi'};
end
if nargin<4
    q = 1;
end
if nargin<5
    SNRin = 40;
end
if nargin<6
    inits = [1,1,0,0,0,0,0];
end
if nargin<7
    with_gt = any(strcmp(gt_name,gt_list));
end
if nargin<8
    rois = [8,1,1;  162,60,8 ; 79,162,16]; % chart_and_stuffed_toy
%     rois = [170,25,1;  160,126,8; 135,192,16]; % balloons
end

poi       = [32,32];
roi_size  = [64,64];

% Number of patches (rows and columns) without counting groundtruth row nor
% the strectal column
n_rows = length(optical_list);
n_cols = size(rois,1);

% TODO: treat experimental cases
%         case 'leaves';    roi = [1,128,195,128];
%         case 'RedYellow'; roi = [35,128,201,128];


n = [256,256,16];
n_str = sprintf('%ix%ix%i',n(1),n(2),n(3));

% Find the gt_center corresponding to the chosen groundtruth
% gt_centers = [  255,128;
%                 256,256;
%                 256,256;
%                 256,256;
%                 230,280;
%                 256,256;
%                 200,236;
%                 256,256  ];
% 
% gt_center = gt_centers(strcmp(gt_list,gt_name),:);
% gt_center = [max(min(gt_center(1),512-n(1)/2),n(1)/2+1),...
%              max(min(gt_center(2),512-n(2)/2),n(2)/2+1)];



algo_name = 'ADMM';
rho = 40;
algstr = sprintf('%.0f_',rho);



wavelengths = linspace(470,620,16);

gt_row = double(with_gt);
U = cell(n_rows+gt_row,1);
if with_gt
    load(sprintf('matfiles/groundtruth/%s_256x256x16_256x256x16.mat',gt_name));
    U{1} = reshape(x, n);
end

for row = 1:n_rows
    % Synthetic examples
    if any(strcmp(gt_name,gt_list))
        optical_str = optical_list{row};
        optical_str_sep = strsplit(optical_str,'_');

        % Optical path is the first part of optical str
        optical_path = optical_str_sep{1};

        % fp_layout is the second part of optical str
        fp_layout = optical_str_sep{2};

        % the multi flag is at the end of optical str.
        multi_snap = strcmp(optical_str_sep{end},'multi');
        if multi_snap
            m = [256,256,q^2];
        else
            m = [256*q,256*q,1];
        end
        m_str = sprintf('%ix%ix%i',m(1),m(2),m(3));

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

        % Compute the corresponding jobid
        jobid = hsics_jobid(gt_name,q,optical_path,fp_layout,multi_snap,...
                            withPSF,psf_rad_max,SNRin);

        identifier = sprintf('HSI_%s_%s',algo_name,optical_path);

        settings_filename = sprintf('%s_%s%s_%s',identifier,algstr,n_str,m_str);

        if inits(row)
            file = [results_folder,filesep,'checkpoint',filesep,identifier,...
                    filesep,'init_',settings_filename,'_j',num2str(jobid),'.mat'];
            file_content = load(file);
            U{row+gt_row} = reshape(file_content.u0,n);
        else
            file = [results_folder,'/results/',identifier,...
                    '/res_',settings_filename,'_j',num2str(jobid),'.mat'];
            file_content = load(file);
            U{row+gt_row} = reshape(file_content.u_est,n);
        end
        
    else % Experimental data
        % Overwrite input parameters for the plot
        optical_list = {'Naive';
                        'Init';
                        'Result'};
        SNRin = 40;
        inits = [0,0,0];
        with_gt = false;
%         q = 2;
%         m = [256*q,256*q,1];
%         m_str = sprintf('%ix%ix%i',m(1),m(2),m(3));
        
        if strcmp(gt_name,'leaves') 
            % Leaves
            res_file=[results_folder,'/results/HSI_ADMM_InFocus/',...
                      'res_HSI_ADMM_InFocus_40_256x256x16_512x512x1_j704.mat'];
            res_file_content=load(res_file);
            rois = [12,12,1; 180,20,8 ; 128,112,16];
        elseif strcmp(gt_name,'RedYellow')
            % RedYellow
            res_file=[results_folder,'/results/HSI_ADMM_InFocus/',...
                      'res_HSI_ADMM_InFocus_40_256x256x16_512x512x1_j705.mat'];
            res_file_content=load(res_file);
            rois = [50, 5, 1; 150,64, 8; 40,120, 16];
        end

        mes_file_content = load(res_file_content.measurements_file);
        %sens_file_content = load(res_file_content.sensing_file);
        n = [256,256,16];
        ns = n(1:2);
        nl = n(3);
        U_est  = reshape(res_file_content.u_est,n);
        U_init = reshape(mes_file_content.u0,n);

        ms = [512,512];
        Y = reshape(mes_file_content.y,ms);
        fp_mosaic = [3   4   2   1;
                     11  12  10  9;
                     15  16  14  13;
                     7   8   6   5];
        idL  = repmat(fp_mosaic,ms./[4,4]);
        idLy = arrayfun(@(j) (idL == j),1:nl,'UniformOutput',false);

        U_naive = zeros(n);
        for band = 1:nl
            U_naive(:,:,band) = imresize(reshape(Y(idLy{band}),128,128),ns,'nearest');
        end
        U = {U_naive,U_init,U_est};
        n_rows = length(optical_list);
        n_cols = size(rois,1);
    end
end


h = figure();
set(h,'position',[210, 136, n_cols*110+120, (n_rows+gt_row)*110]);
clf;
lstyles = {'-','--',':'};

%axs = cell(n_rows+gt_row,3);
axs = tight_subplot(n_rows+gt_row,4,[.02 .02],[.03 .01],[.01 .03]);
for rowid = 1:n_rows+gt_row
    ax_spectral = axs(rowid*(n_cols+1));
    for patch = 1:n_cols
        %axs{rowid,patch} = subplot(n_rows+gt_row,4,patch+(rowid-1)*4);
        %axes(axs(patch+(rowid-1)*(n_cols+1)))
        ax = axs(patch+(rowid-1)*(n_cols+1));
        imagesc(ax,...
                U{rowid}(rois(patch,1)+(0:roi_size(1)-1),...
                         rois(patch,2)+(0:roi_size(2)-1),...
                         rois(patch,3) )); 
        hold(ax, 'on');
        axis(ax,'square', 'off');
        cmap  = colorMapGen(wavelengths(rois(patch,3)), 256);
        
        colormap(ax, cmap );
        plot(ax,poi(1)+0.5,poi(2)+0.5,'+','color',[0,0,0]);
        plot(ax,poi(1),poi(2),'+','color',[1,0.9,0.5]);
        col_letter=['(',char('A'-1+patch),')'];
        if rowid==1
            col_label = sprintf('%d nm',wavelengths(rois(patch,3)));
            text(ax,1.5,1.5,col_label,...
                 'color',[0,0,0],...
                 'verticalalignment','top',...
                 'horizontalalignment','left');
            text(ax,1,1,col_label,...
                 'color',[1,0.9,0.5],...
                 'verticalalignment','top',...
                 'horizontalalignment','left');
            
            text(ax,poi(1)+3.5,poi(2)+3.5,col_letter,...
                 'color',[0,0,0],...
                 'verticalalignment','bottom',...
                 'horizontalalignment','left',...
                 'fontsize',12);
            text(ax,poi(1)+3,poi(2)+3,col_letter,...
                 'color',[1,0.9,0.5],...
                 'verticalalignment','bottom',...
                 'horizontalalignment','left',...
                 'fontsize',12);
        end
        if patch == 1 && (rowid>1 || ~with_gt)
            % Linebreak after optical path
            row_label = sprintf(strrep(optical_list{rowid-gt_row},'Focus_','Focus\n'));
            % replace _ by space
            row_label = strrep(row_label,'_',' ');
            % add 'init' if appropriate
            if inits(rowid-gt_row)
                row_label = [row_label,' init']; %#ok<AGROW>
            end
            text(ax,1.5,roi_size(1)+0.5,row_label,...
                 'color',[0,0,0],...
                 'fontsize',12,...
                 'verticalalignment','bottom',...
                 'horizontalalignment','left');
            text(ax,1,roi_size(1),row_label,...
                 'color',[1,0.9,0.5],...
                 'fontsize',12,...
                 'verticalalignment','bottom',...
                 'horizontalalignment','left');
        end


        %subplot(n_rows+gt_row,n_cols+1,rowid*(n_cols+1));
        x = wavelengths;
        y1 = squeeze(U{rowid}(poi(1)+rois(patch,1),...
                              poi(2)+rois(patch,2),:));
        plot(ax_spectral,x,y1,...
            'k','linestyle',lstyles{mod(patch,3)+1},...
            'linewidth',2'); 
        hold(ax_spectral, 'on');
        text(ax_spectral,x(end),y1(end),col_letter,...
                 'color',[0,0,0],...
                 'verticalalignment','bottom',...
                 'horizontalalignment','left',...
                 'fontsize',10);
        
        % Plot the error
        if with_gt && rowid >1
            y2 =squeeze(U{1}(poi(1)+rois(patch,1),...
                                 poi(2)+rois(patch,2),:));
            plot(ax_spectral,x,y2,...
                'color',[0.5,0.5,0.5],...
                'linestyle',lstyles{mod(patch,3)+1},...
                'linewidth',1);
            X=[x(:).',fliplr(x(:).')];    % create continuous x value array for plotting
            Y=[y1(:).',fliplr(y2(:).')];  % create y values for out and then back
            fill(ax_spectral,X,Y,[1,0,0],'facealpha',0.1,'edgealpha',0); % plot filled area         
        end
    end
    xlim(ax_spectral,[wavelengths(1),wavelengths(end)]);
    ylim(ax_spectral,[0,1]);
    set(ax_spectral,'ytick',[0,1]);
    set(ax_spectral,'xtick',[]);
end
set(ax_spectral,'xtick',wavelengths(rois(:,3)));
end