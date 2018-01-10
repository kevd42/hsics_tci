function out = printPSNRs(u_est,info,it,itmax,iter_div,y,norm_y,...
                          Phi,Psi0,A,om,x,dyn,verbose,graphics,h)
%% PRINTPSNRS 
% Print and plot the PSNR, the SNRin and the analysis L1 norm of the
% current u_est iterate.
% Author : K. Degraux
%  (c) UCLouvain 2018

out = info.qfun_out;

if isempty(out)
    out.SNRin      = zeros(itmax,1);
    out.regularity = zeros(itmax,1);
    
    if ~isempty(x)
        out.SNR    = zeros(itmax,1);
        out.regutarget = norm(om.*(A*(Psi0'*x)),1);
    end
end

out.SNRin(it)      = 20*log10( norm_y / norm(y - Phi  * u_est) );
if isfield(info,'Au')
    Au = info.Au;
else
    Au = A*u_est;
end
out.regularity(it) = norm(om.*Au,1);
if ~isempty(x)
    MSE = sum((x - Psi0 * u_est).^2)/length(x);
    out.PSNR(it)    = 20*log10(dyn) - 10*log10( MSE );
end
if verbose
    fprintf('SNRin = %3.2fdB\n',out.SNRin(it));
    
    if ~isempty(x)
        fprintf('PSNR   = %3.2fdB\n',out.PSNR(it));
        fprintf('Regularity = %g (target %g)\n',out.regularity(it),out.regutarget);
    else
        fprintf('Regularity = %g\n',out.regularity(it));
    end
    
end

if graphics
    softfig(h);
    set(h,'name','quality function');
    if ~isempty(x)
        subplot(1,3,1);
    else
        subplot(1,2,1);
    end
    plot((1:it)*iter_div,out.SNRin(1:it));
    xlim([1,max(2,it)]*iter_div);
    title('SNRin')
    if ~isempty(x)
        subplot(1,3,2);
        plot((1:it)*iter_div,out.PSNR(1:it));
        xlim([1,max(2,it)]*iter_div);
        title('PSNR')
        
        subplot(1,3,3);
        plot([1,it]*iter_div,[out.regutarget, out.regutarget],':k');hold on;
    else
        subplot(1,2,2);
    end
    plot((1:it)*iter_div,out.regularity(1:it));hold off;
    xlim([1,max(2,it)]*iter_div);
    title('Regularity');
    

    drawnow;
end

end