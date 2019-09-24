function [ U ] = Retsina_3( Y,V_init,P1,P2,lamda,iter,fold )
%RETSINA step 3 solve the Nonlinear least squares problem
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

model.variables{1,1}=V_init{1};
model.variables{1,2}=V_init{2};
model.variables{1,3}=V_init{3};
model.factors.A=1;
model.factors.B=2;
model.factors.C=3;
if fold==4 || fold==3
    model.factors.A_tilde1={1, @(z,task)struct_matvec(z,task,P1{1},[])};
    model.factors.A_tilde2={1, @(z,task)struct_matvec(z,task,P1{2},[])};
    model.factors.A_tilde3={1, @(z,task)struct_matvec(z,task,P1{3},[])};
    model.factors.B_tilde1={2, @(z,task)struct_matvec(z,task,P2{1},[])};
    model.factors.B_tilde2={2, @(z,task)struct_matvec(z,task,P2{2},[])};
    model.factors.B_tilde3={2, @(z,task)struct_matvec(z,task,P2{3},[])};
    
    
    model.factorizations.Y1.data=Y{1};
    model.factorizations.Y1.cpd={'A_tilde1','B_tilde1','C'};
    model.factorizations.Y1.weight=lamda(1);
    model.factorizations.Y2.data=Y{2};
    model.factorizations.Y2.cpd={'A_tilde2','B_tilde2','C'};
    model.factorizations.Y2.weight=lamda(2);
    model.factorizations.Y3.data=Y{3};
    model.factorizations.Y3.cpd={'A_tilde3','B_tilde3','C'};
    model.factorizations.Y3.weight=lamda(3);
    
    if fold==4
        model.factors.A_tilde4={1, @(z,task)struct_matvec(z,task,P1{4},[])};
        model.factors.B_tilde4={2, @(z,task)struct_matvec(z,task,P2{4},[])};
        model.factorizations.Y4.data=Y{4};
        model.factorizations.Y4.cpd={'A_tilde4','B_tilde4','C'};
        model.factorizations.Y4.weight=lamda(4);
    end
else
    fprintf('Error. Use 3-fold or 4-fold accelaration')
end
% sdf_check(model,'print');
[sol,output]=sdf_nls(model,'Display',0,'MaxIter',iter);
U{1}=sol.factors.A;
U{2}=sol.factors.B;
U{3}=sol.factors.C;

end