% HSICS_SETUP
% Set the path for hsics_main
%
% see README.txt
%
% Author: K. Degraux
% Date: 2018
% see also hsics_main, hsics_jobid
%
%  (c) UCLouvain 2018

% Check for the presence of the CAVE dataset
% see http://www.cs.columbia.edu/CAVE/databases/multispectral/
if ~exist(['images',filesep,'complete_ms_data'], 'dir')
    warning(['please download the CAVE dataset at \n',...
             'http://www.cs.columbia.edu/CAVE/databases/',...
             'multispectral/zip/complete_ms_data.zip\n',...
             'and unzip it in the images folder'])
end

% Directory containing hsics_setup.m
hsics_dir = which('hsics_setup');
hsics_dir = hsics_dir(1:end-14);

addpath([hsics_dir,filesep,'utils'])
addpath([hsics_dir,filesep,'utils',filesep,'prox'])
addpath([hsics_dir,filesep,'utils',filesep,'imresizemex'])

addpath([hsics_dir,filesep,'toolboxes'])
addpath([hsics_dir,filesep,'toolboxes',filesep,'spot-master'])
addpath([hsics_dir,filesep,'toolboxes',filesep,'tight_subplot'])

addpath([hsics_dir,filesep,'sensingop'])

addpath([hsics_dir,filesep,'scripts'])

addpath([hsics_dir,filesep,'optim'])


