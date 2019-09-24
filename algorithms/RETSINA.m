function [U_hat] = RETSINA(X_sampled,Y,S1,S2,P1,P2,n,F,iter1,iter2,iter3)
%RETSINA algorithm
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
%% step 1: Initialization
fprintf('step 1: Initialization \n')
U_init = RETSINA_initialization(X_sampled,Y,S1,S2,n,F,iter1);
clear X_sampled
%% step 2
fprintf('step 2: Refinement \n')
[U_ref]=Retsina_refine(Y,U_init,S1,S2,F,iter2,iter2);
clear U_init
%% step 3
fprintf('step 3: Complete the tensor \n')
[ U_hat] = Retsina_3( Y,U_ref,P1,P2,[1 1 1 1],iter3,n);
end

