% Demo for RETSINA algorithm
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

% this demo is developed for 3-fold accelaration. For other types of
% acceleration the code has to be modified.

clear all; close all; clc;
T;% tensor k_y x k_x x coils x frames
[I_init,J_init,K_init,L]=size(T);
for k=1:K_init
    for l=1:L
        X(:,l,k)=reshape(T(:,:,k,l),[I_init*J_init,1]);        
    end
end
[I,J,K]=size(X); %k_space x frame x coil
clear T
%% sample tensor
n=3; % n-fold acceleration sample ratio=1/n
[Y,P1,P2,S1,S2 ] = fiber_sample_fmri( X,n,I_init,J_init);
%% create mask

mask=zeros(I,J,K); %sampling mask (missing entries are equal to 0)
for i=1:n
    mask(S1{i},S2{i},:)=1;
end

%% set tensor rank
F=100; 
%% step 1: Initialization
tic
iter1=50;
iter2=2;
iter3=5;
[X_hat,U_hat] = RETSINA(mask.*X,Y,S1,S2,P1,P2,n,F,iter1,iter2,iter3);
tim=toc;
%% results
tim_min=tim/60;
fprintf('execution time: %3.1f minutes \n',tim_min)
nre=frob(X-X_hat)/frob(X);
fprintf('RETSINA NRE=%3.3f \n',nre)
