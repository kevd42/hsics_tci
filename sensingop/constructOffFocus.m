function [Phi,info] = constructOffFocus(m,n,slm,idLy,SR,upfactor)
% Phi = CONSTRUCTOFFFOCUS(slm,psf,idLy,SR,m,n,s,upfactor)
% Pre-Compute useful arrays and return an opFunction
%  INPUTS:
%   slm : Array of size [si,sj,mS] or of size [si,sj,mS,nf]
%         to be convolved with every band of the 
%         input for every snapshot.
%         No upscaling of slm is performed!
%         A psf can be pre-convolved with this array.
%         The psf can be different for every filter.
%         (A more accurate model would be to convolve before applying SR.)
%   idLy : cell of logical arrays of size ms that contain, for each band,
%          the location of the corresponding filters on the FPA.
%   SR : Spectral response matrix to be applied on every pixel spectrum.
%        Size : nl x nf (nl = n(3) = number of wl, nf = number of filters)
%        If empty, defaults as 1.
%   m : output size [mi, mj, mS] 
%       with ms = [mi mj] spatial dimensions of the FPA.
%       and  mS = number of snapshots.
%   n : input size [ni, nj, nl]
%       with ns = [ni, nj] spatial dimensions of the target image.
%       and  nl = spectral dimension (number of wl).
%   upfactor : multiplicative factor between the target size and the actual
%              slm resolution. The input is upscaled before computing the
%              convolution.
%  OUTPUT:
%   Phi : spot opFunction that models the OffFocus device behavior.
%
%
%  NOTE : This function should work in the Panchromatic case as well, 
%         setting n(3)=1
% 
% Author: Kevin Degraux
%  (c) UCLouvain 2018


ns    = n(1:2);
ns_up = upfactor.*ns;
ms    = m(1:2);
mS    = m(3);

% Crop the slm pattern and pre-compute the convolution kernel fft
s      = [size(slm,1),size(slm,2)];
s_eff  = ns_up + ms - 1;
offset = floor(s/2 - s_eff/2);

% Extract the number of filters.
if isscalar(SR)
    nf = n(3);
else
    nf = size(SR,2);
end

if ndims(slm) < 4
    SLM = zeros([s_eff,mS]);
    for snap = 1:mS
        slm_eff = slm( (1:s_eff(1)) + offset(1),...
                       (1:s_eff(2)) + offset(2),snap);
        SLM(:,:,snap) = fft2(slm_eff,s_eff(1),s_eff(2));
    end
else
    SLM = zeros([s_eff,mS,nf]);
    for snap = 1:mS
        for band = 1:nf
            slm_eff = slm( (1:s_eff(1)) + offset(1),...
                           (1:s_eff(2)) + offset(2),snap,band);
            SLM(:,:,snap,band) = fft2(slm_eff,s_eff(1),s_eff(2));
        end
    end
end

% Valid part of the convolution
idx_valid = {(ns_up(1) : s_eff(1));
             (ns_up(2) : s_eff(2))}; 
         
% Central part of the valid part
offset    = floor((s_eff-ns_up+1)/2 - ms/2 );
fpa_crop  = {(1:m(1)) + offset(1);
             (1:m(2)) + offset(2)}; 
         
% Combine the cropping regions
fused_idx = {idx_valid{1}(fpa_crop{1});
             idx_valid{2}(fpa_crop{2})};



% Initialize the upscaling operator. Upscaling is performed after the
% application of SR (most probably reducing the dimension).
if upfactor == 1
    data.init = [];
else
    data.init = opImresize_init([ns,nf],upfactor,@(x) lanczos(x,3),6,false);
end
% Record other data used at running time.
data.SLM        = SLM;
data.fused_idx  = fused_idx;
data.s_eff      = s_eff;
data.ns_up      = upfactor.*ns;
data.n          = n;
data.m          = m;
data.idLy       = idLy;
data.SR         = SR;
data.SRSRT      = SR*SR';

M = prod(m);
N = prod(n);

    
if SR == 1 && upfactor == 1 
    nn = prod(s_eff)*nf;
    mm = prod(s_eff)*nf*mS;
    info.nn = nn;
    info.mm = mm;
    info.s = s_eff;
    fpa_idx = reshape(1:mm,[s_eff,mS,nf]);
    fpa_idx = fpa_idx(data.fused_idx{1},data.fused_idx{2},:,:);
    info.Rm_idx = zeros(m);
    for snap = 1:mS
        Rm_snap = zeros(ms);
        for band = 1:nf
            fpa_idx_band = fpa_idx(:,:,snap,band);
            Rm_snap(idLy{band}) = fpa_idx_band(idLy{band});
        end
        info.Rm_idx(:,:,snap) = Rm_snap;
    end
    info.Rm_idx = info.Rm_idx(:);
    img_idx = reshape(1:nn,[s_eff,nf]);
    info.Rn_idx = img_idx(1:n(1),1:n(2),:);
    info.Rn_idx = info.Rn_idx(:);
    
    clearvars -except data M N SLM mm nn s_eff nf info;
    info.normPhi2 = max(max(max(abs(sum(SLM.*conj(SLM),3)))));
    if ndims(SLM)<4
        SLM = repmat(SLM,1,1,1,nf);
    end
    info.invPhi_T_Phi_rhoId_h = @(rho) opFunction(nn,nn,@(x,mode) ...
           apply_invPhi_T_Phi_rhoId(x,mode,data,...
           1./(reshape(sum(SLM.*conj(SLM),3),s_eff(1),s_eff(2),nf)+ rho)),...
           false,true);
    clearvars -except data M N mm nn info;
    info.FTDF = opFunction(mm,nn,@(x,mode)apply_FTDF(x,mode,data),false,true);

end


clearvars -except data M N info; % so that useless data is not stored in the spot operator

Phi = opFunction(M,N,@(x,mode)apply_operator(x,mode,data),false,true);


clearvars -except Phi info;

end
%%
function out = apply_operator(in,mode,data)
    
    % Extract the data
    SLM        = data.SLM;
    fused_idx  = data.fused_idx;
    s_eff      = data.s_eff;
    ns_up      = data.ns_up;
    n          = data.n;
    m          = data.m;
    idLy       = data.idLy;
    SR         = data.SR; % Spectral response matrix. Default is 1
    init       = data.init;
    
    % Number of filters
    if isscalar(SR)
        nf  = n(3);
    else
        nf  = size(SR,2);
    end

    
    if mode == 1
        % Apply spectral response matrix on wavelengths
        in   = reshape(in, n(1)*n(2), n(3)) * SR; 
        
        % Upscale x
        if isempty(init)
            x_up = reshape(in,[ns_up,nf]);
        else
            x_up = reshape(opImresize( in , 1 ,init),[ns_up,nf]);
        end
        
        % Compute the fft of the upscaled x
        X_UP = zeros( [s_eff , nf] );
        for band=1:nf
            X_UP(:,:,band) = fft2( x_up(:,:,band) , s_eff(1) , s_eff(2) );
        end
        
        % For every snapshot:
        % - Perform convolution in Fourier domain
        % - go back to spatial domain with ifft
        % - Crop the region of interest
        % - Map the filtered pixel on the final FPA
        y      = zeros(m);
        y_snap = zeros(m(1:2));
        for snap = 1:m(3)
            for band=1:nf
                if ndims(SLM) < 4
                    y_full = ifft2( X_UP(:,:,band) .* SLM(:,:,snap) ,'symmetric');
                else
                    y_full = ifft2( X_UP(:,:,band) .* SLM(:,:,snap,band) ,'symmetric'); % symmetric because real
                end
                y_fpa  = y_full( fused_idx{:} );
                y_snap(idLy{band}) = y_fpa(idLy{band});
            end
            y(:,:,snap) = y_snap;
        end

        out = y(:);
    
    else % mode == 2
        
        y    = reshape(in,m);
        
        % For every snapshot:
        % - Extract the pixel locations corresponding to every filter
        % - Zero-pad the FPA to match convolution dimensions
        % - Perform the conjugate convolution in Fourier domain
        % - Go back to spatial domain, crop the result to valid indices and
        %   add to accumulator
        x_up = zeros([ns_up,nf]);
        y_full = zeros(s_eff);
        for band=1:nf
            X_FULL = zeros(s_eff);
            for snap = 1:m(3)
                y_snap = y(:,:,snap);
                y_fpa  = zeros(m(1:2));
                y_fpa(idLy{band}) = y_snap(idLy{band});
                y_full( fused_idx{:} ) = y_fpa;
                if ndims(SLM) < 4
                    X_FULL = X_FULL + fft2(y_full,s_eff(1) , s_eff(2)) .* conj(SLM(:,:,snap));
                else
                    X_FULL = X_FULL + fft2(y_full,s_eff(1) , s_eff(2)) .* conj(SLM(:,:,snap,band));
                end
            end
            x_full = ifft2( X_FULL,'symmetric'); % symmetric because real
            x_up(:,:,band) =  x_full(1:ns_up(1),1:ns_up(2));
        end

        
        % Perform the adjoint of upscaling
        if isempty(init)
            out = x_up;
        else
            out = opImresize( x_up , 2 ,init);
        end

        % Apply the transpose of SR to spectral response coefficient to go
        % back in the spectral domain.
        out = reshape(out, n(1)*n(2), nf) * SR';
        out = out(:);
        
    end
    
end

function out = apply_FTDF(in,mode,data)

    % This matrix is NOT symmetric if mS>1

    % Extract the data
    SLM    = data.SLM;
    s      = data.s_eff;
    
    nf  = data.n(3);
    mS  = data.m(3);

    if mode == 1

        x = reshape(in, [s,nf]);
        
        % Compute the fft of x
        X = zeros( [s, nf] );
        for band=1:nf
            X(:,:,band) = fft2( x(:,:,band) , s(1) , s(2) );
        end
        
        % For every snapshot:
        % - Perform convolution in Fourier domain
        % - go back to spatial domain with ifft
        y = zeros([s,mS,nf]);
        for snap = 1:mS
            for band=1:nf
                if ndims(SLM) < 4
                    y(:,:,snap,band) = ifft2( X(:,:,band) .* SLM(:,:,snap) ,'symmetric');
                else
                    y(:,:,snap,band) = ifft2( X(:,:,band) .* SLM(:,:,snap,band) ,'symmetric'); % symmetric because real
                end
            end
        end
        out = y(:);
    
    else % mode == 2
        
        y = reshape(in,[s,mS,nf]);
        % For every snapshot:
        % - Perform the conjugate convolution in Fourier domain
        % - Go back to spatial domain, and
        %   add to accumulator
        x = zeros([s, nf]);
        for band=1:nf
            X_BAND = zeros(s);
            for snap = 1:mS
                y_snap = y(:,:,snap,band);
                if ndims(SLM) < 4
                    X_BAND = X_BAND + fft2(y_snap,s(1) , s(2)) .* conj(SLM(:,:,snap));
                else
                    X_BAND = X_BAND + fft2(y_snap,s(1) , s(2)) .* conj(SLM(:,:,snap,band));
                end
            end
            x(:,:,band) = ifft2( X_BAND,'symmetric'); % symmetric because real
        end
        out = x(:);
    end


end


function out = apply_invPhi_T_Phi_rhoId(in,~,data,invSLM_rhoId)

    s   = data.s_eff;
    nf  = data.n(3);

    x = reshape(in, [s,nf]);

    % Compute the fft of x
    y = zeros([s,nf]);
    for band=1:nf
        X_BAND = fft2( x(:,:,band) , s(1) , s(2) );
        y(:,:,band) = ifft2( X_BAND .* invSLM_rhoId(:,:,band) ,'symmetric');
    end

    out = y(:);

end

