function hsics_delete_checkpoint(results_folder,jobid,verbose)
% HSICS_DELETE_CHECKPOINT
% Deletes the backup file extract convergence information from checkpoint 
% file if available.
% Only compatible with UNIX environment.
% Author: K. Degraux
% Date: 10/1/2018
% (c) UCLouvain 2018
if nargin<3
    verbose=false;
end
try
    cp_file   = ls([results_folder,'/checkpoint/*/cp*j',num2str(jobid),'.mat']);
    cp_file = cp_file(:,1:end-1);
    cp_file_content = load(cp_file);
    convergence_file_content = struct();
    if isfield(cp_file_content,'rel_err')
        convergence_file_content.rel_err=cp_file_content.rel_err;
    end
    if isfield(cp_file_content,'info')
        convergence_file_content.info=cp_file_content.info; %#ok<STRNU>
    end
    if verbose
        fprintf('compressed cp file %i (checkpoint deleted)\n',jobid);
    end
    save(cp_file,'-struct','convergence_file_content');
catch
    if verbose
        fprintf('cp file %i not found \n',jobid);
    end
end
try
    bak_file   = ls([results_folder,'/checkpoint/*/cp*j',num2str(jobid),'.mat.bak']);
    bak_file = bak_file(:,1:end-1);
    rm(bak_file)
    if verbose
        fprintf('removed bak file %i\n',jobid);
    end
catch
    if verbose
        fprintf('bak file %i not found \n',jobid);
    end
end

end