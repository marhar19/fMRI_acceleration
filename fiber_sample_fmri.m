function [ Y,P1,P2,S1,S2 ] = fiber_sample_fmri( X,ratio,I_init,J_init)
% fiber_sample_fmri samples 1/n equispaced frequencies in the k-space
% to perform n-fold accelaration.

[I,J,K]=size(X);
PP1=eye(I_init,I_init);
for i=1:ratio
    S{i}=i:ratio:I_init; %sampling index for k_y frequencies
    P{i}=PP1(i:ratio:end,:); %row selection matrix of S
end


PP2=eye(J,J);
S2{1}=1:ratio:J;
P2{1}=PP2(S2{1},:);
for i=2:ratio
    S2{i}=[1,i:ratio:J]; %sampling index for frames (time-slots)
    P2{i}=PP2(S2{i},:); %row selection matrix of S2
end

for i=1:ratio
    P1{i}=kron(eye(J_init),P{i}); %sampling index for concatenated k-space
    S1{i}=find(sum(P1{i})==1); %row selection matrix of S1
end

for i=1:ratio
    Y{i}=X(S1{i},S2{i},:); %the formed sub-tensors
end
end

