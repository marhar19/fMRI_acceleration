% Demo for MS-RETSINA algorithm
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
clear all; close all; clc;
T; % tensor k_y x k_x x coil x slice x frame
[I_init,J_init,K_init,L_init,M_init]=size(T); 
%% covert 5-way array to 3rd order tensor
X=reshape(T,[I_init*J_init,K_init*L_init,M_init]); %k_space x coil-slice x frame 
%(In the TSP paper the tensor modes are k_space x frame x coil-slice)
clear T
[I,J,K]=size(X);
%% sample tensor
s=2;ky=4; % sample 1/k ky lines and 1/s slices (rs)-fold acceleration
[Y,P1,P2,P3,S1,S2,S3 ] = entry_sample_fmri( X,ky,s,[I_init,J_init,K_init,L_init,M_init]);
%% create mask
mask=zeros(I,J,K);
k=0;
for i=1:ky
    for j=1:s
        k=k+1;
        mask(S1{i},S2{j},S3{k})=1;   %sampling mask (missing entries are equal to 0)
    end
end
%% set tensor rank
F=100; %tensor rank
%% MS-RETSINA
iter1=50;
iter2=20;
[X_hat,U_hat] = MS_RETSINA(mask.*X,Y,S1,S2,S3,ky,s,F,iter1,iter2);

% frob(X_hat-X)/frob(X)
sum2=0;sum1=0;
for i=1:K
    sum2=sum2+norm(X(:,:,i)-X_hat(:,:,i),'fro');
    sum1=sum1+norm(X(:,:,i),'fro');
end
fprintf('MS-RETSINA NRE with %d-fold acceleration: %1.6f \n',s*ky,sum2/sum1)

