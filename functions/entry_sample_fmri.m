function [Y,P1,P2,P5,S1,S2,S5 ] = entry_sample_fmri( X,ky,s,size5)
% entry_sample_fmri samples 1/ky equispaced frequencies in the k-space and
% 1/s slices to perform (sky)-fold accelaration.
I_init=size5(1);J_init=size5(2);K_init=size5(3);L_init=size5(4);M_init=size5(5);

PP1=eye(I_init,I_init);
for i=1:ky
    S1{i}=i:ky:I_init; %sampling index for k_y frequencies
    P{i}=PP1(i:ky:end,:);
    P1{i}=kron(eye(J_init),P{i});
    S1{i}=find(sum(P1{i})==1);
end

PP4=eye(L_init,L_init);
for i=1:s
    S4{i}=i:s:L_init; %sampling index for slices
    P4{i}=PP4(i:s:end,:);
    P2{i}=kron(P4{i},eye(K_init));
    S2{i}=find(sum(P2{i})==1);
end

PP5=eye(M_init,M_init);
for i=1:(ky*s)
    S5{i}=[1,i+1:(ky*s):M_init]; %sampling index for frames (time-slots)
    P5{i}=PP5(S5{i},:);
end


k=0;
for i=1:ky
    for j=1:s
        k=k+1;
        Y{k}=X(S1{i},S2{j},S5{k}); %the formed sub-tensors
    end
end
end

