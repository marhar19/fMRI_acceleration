function [V] = Retsina_refine(Y,U0,S1,S2,F,it1,it2)
%Retsina step 2: refinement
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

%% run CPD of Y{2} to compute factor C
% run only a few iterations to maintain common perm and scaling
T=Y{2};U1{1}=U0{1}(S1{2},:);U1{2}=U0{2}(S2{2},:);U1{3}=U0{3};
U=cpd(T,U1,'MaxIter',it1);
C=U{3};

%% run CPD of Y with C fixed
% run only a few iterations to maintain common scaling
model = struct;
model.factors.A = 'a';
model.factors.B = 'b';
model.factors.C = C; % The third factor is constant.
model.factorizations.myfac.cpd = {'A','B','C'};
I=0;J=0;
for i=1:length(Y)
    I=I+size(Y{i},1);
    J=J+size(Y{i},2);
end
J=J-length(Y)+1;

A=zeros(I,F);B=zeros(J,F);

for i=1:length(Y)
    if i~=2
        model.variables.a = U0{1}(S1{i},:);
        model.variables.b = U0{2}(S2{i},:);
        model.factorizations.myfac.data = Y{i};
        
        % Solve the SDF model.
        options.MaxIter=it2;
        sol = sdf_nls(model,options);
        A(S1{i},:)=sol.variables.a;B(S2{i},:)=sol.variables.b;
        clear sol
    end
end

A(S1{2},:)=U{1};B(S2{2},:)=U{2};

V{1}=A;V{2}=B;V{3}=C;

end
