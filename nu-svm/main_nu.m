function [R1,R2,alpha_nu,accuracy,auc,specificity,gmeans,Number_SVs,T,PY,rho]=main_nu(tol,kernel_type,kernel_param,nu,traX,traY,tstX,tstY)
[l,~]=size(traX);
C=1/l;
traX=[traX,ones(l,1)];
K = calckernel(kernel_type,kernel_param,traX,traX);
Q=diag(traY)*K*diag(traY);
Lengthnu=length(nu);
alpha_nu=zeros(l,Lengthnu);
Number_SVs=zeros(Lengthnu,3);
R1=zeros(l,Lengthnu);
R2=zeros(l,Lengthnu);
S = zeros(l,1);
rho = zeros(1,Lengthnu);
s=1;

for s=1:Lengthnu 
    tic
    [alpha_nu(:,s)]=QPP_cal([],[],Q,traY,nu(s));
    T(s,1)=toc;
    R1(:,s)=alpha_nu(:,s)>0-tol & alpha_nu(:,s)<0+tol;
    R2(:,s)=alpha_nu(:,s)>C-tol & alpha_nu(:,s)<C+tol;
    S = alpha_nu(:,s)>=0+tol & alpha_nu(:,s)<=C-tol;
    Number_SVs(s,1)=nnz(R1(:,s));
    Number_SVs(s,2)=nnz(R2(:,s));
    Number_SVs(s,3)=sum(Number_SVs(s,1:2));
    rho(s) =  mean(Q(S,:)*alpha_nu(:,s));
end
[accuracy,auc,specificity,gmeans,PY]=prediction(kernel_type,kernel_param,traX,traY,tstX,tstY,alpha_nu,nu);
