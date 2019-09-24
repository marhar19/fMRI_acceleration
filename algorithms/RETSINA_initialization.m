function [U_hat] = RETSINA_initialization(X_sampled,Y,S1,S2,ratio,F,iter)
% RETSINA step 1: initialization
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
[~,J,~]=size(X_sampled);
%% sum every n vertical slabs
for k=1:floor((J-1)/ratio)-1
    X_inter(:,k,:)=squeeze(sum(X_sampled(:,(k-1)*ratio+4:k*ratio+3,:),2));
end
%% CPD of the approximate tensor
U=cpd((X_inter),F,'MaxIter',iter);
clear X_inter

A=U{1};C=U{3};

%% approximate B using Y
for i=1:ratio
    Y1_2=reshape(permute(Y{i},[2 1 3]),[size(Y{i},2),size(Y{i},3)*size(Y{i},1)]).';
    A_tilde=A(S1{i},:);
    temp=khatri_rao(C,A_tilde);
    B_hat{i}=(temp\Y1_2).';
end
B=zeros(J,F);
for i=1:ratio
B(S2{i},:)=B_hat{i};
end
U_hat{1}=A;U_hat{2}=B;U_hat{3}=C;
clear A B C
end
