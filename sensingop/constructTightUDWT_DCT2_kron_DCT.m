function [A,info] = constructTightUDWT_DCT2_kron_DCT(n)
%% CONSTRUCTTIGHTUDWT_DCT2_KRON_DCT
% [A,info] = CONSTRUCTTIGHTUDWT_DCT2_KRON_DCT(n)
% Constructs the analysis operator corresponding to the kronecker product
% of a DB4 UDWT (+ DCT2 for scaling functions) and a DCT for the spectral
% domain.
% Author : K. Degraux
% Date : 7 Apr. 2017
%  (c) UCLouvain 2018


family    = 'daubechies';
filter    = 4;
lvl       = 3;
redundant = true;
tightframe = true;
padmethod = 'zero'; 
UDWT      = opWavelet2(n(1),n(2),family,filter,lvl,redundant,[],tightframe,padmethod);


Bs        = size(UDWT,1)/(3*lvl+1); % Number of wavelet coefficients per wavelet level-direction ("band")
[bsi,bsj] = findIntegerRoot(Bs); % Each level-direction is reshaped as the closest square

upFactor  = lvl*3+1;
% Ps        = Bs*upFactor;
% A_xy      = opBlockDiag(opDCT2(bsi,bsj),opEye(Ps-Bs)) * UDWT;



% Number of level-directions ("bands") in the wavelet decomposition
% reshaped as the closest integer square.
[pi,pj] = findIntegerRoot(upFactor);

info.bs = [bsi,bsj];
info.ps = [pi,pj];
info.om = [repmat(UDWT.rowscaling,n(3),1)];

data.UDWT  = UDWT;
data.n     = n;
data.Bs    = Bs;
data.Fdct2 = opDCT2(bsi,bsj);
% It is much faster to explicitly compute the 16x16 matrix than using spot
% interface for the wavelet transform
% data.Fdct  = double(opWavelet(n(3),'daubechies',4,3,false));%opDCT(n(3));
% ==> surprisingly seems to give lesser performances

data.Fdct = opDCT(n(3));


M = (size(UDWT,1))*n(3);
N = (size(UDWT,2))*n(3);
% Aref = opBlockDiag(opBlockDiag(fcut,A_xy) , opEye(nf));

clearvars -except data info M N ; % so that useless data is not stored in the spot operator

A = opFunction(M,N,@(x,mode)apply_operator(x,mode,data),false,true);
clearvars -except A info ;


end

function out = apply_operator(in,mode,data)
UDWT = data.UDWT;
Bs   = data.Bs;
n    = data.n;
Fdct2= data.Fdct2;
Fdct = data.Fdct;

if mode == 1
    u = in;
    
    U = reshape(u,[],n(3));
    
    Al = UDWT * U;
    Al(1:Bs,:) = Fdct2 * Al(1:Bs,:);
    
    Al = Al * Fdct';

    out = Al(:);
        
else % mode ==2
    al  = in;
    
    Al = reshape(al,[],n(3));
    
    Al = Al * Fdct;
    
    Al(1:Bs,:) = Fdct2' * Al(1:Bs,:);
    U  = UDWT' * Al;
    
    out = U(:);
end

end















