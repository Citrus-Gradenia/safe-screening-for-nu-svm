 function [alpha]=QPP_cal(R1,R2,Q,traY,nu)
[l,~]=size(Q);
C=1/l;
if nnz(R1)+nnz(R2)==0
    f = zeros(l,1);
    A = -ones(1,l);
    b = -nu;
    Aeq = [];
    beq = [];
    lb = zeros(l,1);
    ub = C*ones(l,1);
    [alpha]=quadprog(Q,f,A,b,Aeq,beq,lb,ub);
else
    id_3=(R1~=1 & R2~=1);
    l=nnz(id_3);  
    Qnew=Q(id_3,id_3);
    f=C*sum(Q(id_3==1,R2==1),2);
    A = -ones(1,l);
    b = -nu+nnz(R2==1) * C;
    Aeq = [];
    beq = [];
    lb = zeros(l,1);
    ub = C*ones(l,1);
    [alpha]=quadprog(Qnew,f,A,b,Aeq,beq,lb,ub);
end
