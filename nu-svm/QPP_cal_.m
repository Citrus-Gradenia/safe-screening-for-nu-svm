 function [alpha]=QPP_cal_(R1,R2,Q,traY,nu)
[l,~]=size(Q);
C=1/l;
eps = 1e-4; 
iter = 100; 
if nnz(R1)+nnz(R2)==0
    lb = zeros(l,1);
    ub = C*ones(l,1);
    [alpha]=DCDM(Q,lb,ub,nu,l,eps,iter);
else
    id_3=(R1~=1 & R2~=1);
    l=nnz(id_3);  
    Qnew=Q(id_3,id_3);
    nu_new = nu - nnz(R2==1) * C;
    lb = zeros(l,1);
    ub = C*ones(l,1);
    [alpha]=DCDM(Qnew,lb,ub,nu_new,l,eps,iter);
end
