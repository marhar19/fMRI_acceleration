function [U_hat] = MS_RETSINA_initialization(X_sampled,Y,S1,S2,S3,ky,s,F,iter)
% MS-RETSINA step 1: initialization
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
[~,~,K]=size(X_sampled);
n=ky*s;
%% sum every n vertical slabs
for k=1:floor((K-1)/n)
    X_inter(:,:,k)=squeeze(sum(X_sampled(:,:,(k-1)*n+2:k*n+1),3));
end
%% CPD of the approximate tensor
options.Algorithm=@cpd_als;
U=cpd((X_inter),F,options,'MaxIter',iter);
clear X_inter

A=U{1};B=U{2};
clear U
%% approximate C using Y
k=0;
for i=1:ky
    A_tilde=A(S1{i},:);
    for j=1:s
        B_tilde=B(S2{j},:);
        k=k+1;
        Y3=reshape(Y{k},[size(Y{k},1)*size(Y{k},2),size(Y{k},3)]);
        
        temp=khatri_rao(B_tilde,A_tilde);
        C_hat{k}=(temp\Y3).';
        clear temp Y3 B_tilde
    end
    clear A_tilde
end

C=zeros(K,F);
for i=1:(s*ky)
    C(S3{i},:)=C_hat{i};
end
clear C_hat
U_hat{1}=A;U_hat{2}=B;U_hat{3}=C;
clear A B C
end

