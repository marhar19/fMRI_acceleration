function [X_hat,U_hat] = MS_RETSINA(X_sampled,Y,S1,S2,S3,ratio_k,ratio_s,F,iter1,iter2)
%MS-RETSINA algorithm
% (c) Charilaos I. Kanatsoulis, University of Minnesota, Sep 24 , 2019
% kanat003@umn.edu
% 
% Reference 1: C.I. Kanatsoulis, X. Fu, N.D. Sidiropoulos and M. Akçakaya, 
%``Tensor Completion from Regular Sub-Nyquist Samples,''
% arXiv preprint

% Reference 2: C.I. Kanatsoulis, N.D. Sidiropoulos, M. Akçakaya and X. Fu, 
%``Regular sampling of tensor signals: Theory and application to fMRI,''
% IEEE International Conference on Acoustics, Speech
% and Signal Processing (ICASSP), 2019
[I,J,K]=size(X_sampled);
%% step 1: Initialization
fprintf('step 1: Initialization \n')
[U_init] = MS_RETSINA_initialization(X_sampled,Y,S1,S2,S3,ratio_k,ratio_s,F,iter1);
%% step 2
fprintf('step 2: Complete the tensor \n')
mask=nan(I,J,K);
ind1=find(X_sampled~=0);
mask(ind1)=1;
[U_hat]=cpd(fmt(mask.*X_sampled),U_init,'MaxIter',iter2);
clear mask
%% reconstruct
X_hat=cpdgen(U_hat);
X_hat(ind1)=X_sampled(ind1);
end

