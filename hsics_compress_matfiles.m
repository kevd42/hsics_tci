function hsics_compress_matfiles(results_folder,jobid,verbose)
% HSICS_COMPRESS_MATFILES
% Delete superfluous fields in results and init matfiles
% Only compatible with UNIX environment.
% Author: K. Degraux
% Date: 10/1/2018
% (c) UCLouvain 2018
if nargin<3
    verbose=false;
end
try
    res_file  = ls([results_folder,'/results/*/*j',num2str(jobid),'.mat']);
    res_file  = res_file(1,1:end-1);
    res_file_content = load(res_file);
    if isfield(res_file_content,'alg_info')
        res_file_content = rmfield(res_file_content,'alg_info'); %#ok<NASGU>
    end
    if verbose
        fprintf('compressed res file %i\n',jobid);
    end
    save(res_file,'-struct','res_file_content');
catch
    if verbose
        fprintf('res file %i not found \n',jobid);
    end
end
try
    init_file = ls([results_folder,'/checkpoint/*/init*j',num2str(jobid),'.mat']);
    init_file = init_file(:,1:end-1);
    init_file_content = load(init_file);
    if isfield(init_file_content,'loadMes')
        init_file_content = rmfield(init_file_content,'loadMes'); 
    end
    if isfield(init_file_content,'computeU0')
        init_file_content = rmfield(init_file_content,'computeU0');
    end
    if isfield(init_file_content,'graphicsU0')
        init_file_content = rmfield(init_file_content,'graphicsU0');
    end
    if isfield(init_file_content,'u0_info')
        init_file_content = rmfield(init_file_content,'u0_info'); %#ok<NASGU>
    end
    if verbose
        fprintf('compressed init file %i\n',jobid);
    end
    save(init_file,'-struct','init_file_content');
catch
    if verbose
        fprintf('init file %i not found \n',jobid);
    end
end


end