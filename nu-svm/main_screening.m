function [RR1,RR2,alpha,Acc,Auc,specificity,gmeans,Num_screening,Time,PY,rho]=main_screening(tol,kernel_type,kernel_param,nu,traX,traY,tstX,tstY)
[l,~]=size(traX);
C=1/l;
traX=[traX,ones(l,1)];
K = calckernel(kernel_type,kernel_param,traX,traX);
Q=diag(traY)*K*diag(traY);
Lengthnu=length(nu);
alpha=-2*ones(l,Lengthnu);
Num_screening=zeros(Lengthnu,3);
RR1=zeros(l,Lengthnu);
RR2=zeros(l,Lengthnu);
S = zeros(l,1);
rho = zeros(1,Lengthnu);
InitalTime=zeros(Lengthnu,1);
DeltaTime=zeros(Lengthnu,1);
RuleTime=zeros(Lengthnu,1);
ReducedQPPTime=zeros(Lengthnu,1);
TotalTime = zeros(Lengthnu,1);


s = 1;
tic
[alpha(:,s)]=QPP_cal([],[],Q,traY,nu(s));
RR1(:,s)=alpha(:,s)>0-tol & alpha(:,s)<0+tol;
RR2(:,s)=alpha(:,s)>C-tol & alpha(:,s)<C+tol;
S = alpha(:,s)>=0+tol & alpha(:,s)<=C-tol;
Num_screening(s,1)=nnz(RR1(:,s));
Num_screening(s,2)=nnz(RR2(:,s));
Num_screening(s,3)=sum(Num_screening(s,1:2));
rho(s) =  mean(Q(S,:)*alpha(:,s));
    
InitalTime(1)=toc;
TotalTime(1) = InitalTime(1);
delta0 = eps*ones(l,1);
for s=2:Lengthnu
    s
    [R1,R2,delta1,DeltaTime(s),RuleTime(s),E(s),Ebar_i0(s),Ebar_j0(s),D(s),length_delta1_Nbar(s)]=DVI(kernel_type,s,delta0,Q,alpha(:,s-1),nu(s),traY);
    alpha(R1==1,s)=0;
    alpha(R2==1,s)=C; 
    nnz(RR1)
    nnz(RR1)
    tic
    if nnz(alpha(:,s)==-2)~=0
        [alpha(alpha(:,s)==-2,s)]=QPP_cal(R1,R2,Q,traY,nu(s)); 
    end
    ReducedQPPTime(s)=toc;
    
    Num_screening(s,1)=nnz(R1);
    Num_screening(s,2)=nnz(R2);
    Num_screening(s,3)=nnz(R1)+nnz(R2);
    RR1(:,s)=R1;
    RR2(:,s)=R2;
    S = alpha(:,s)>=0+tol & alpha(:,s)<=C-tol;
    rho(s) =  mean(Q(S,:)*alpha(:,s));
    delta0 = delta1;
    MoreInf = [Ebar_i0;Ebar_j0;D; E;(Ebar_i0+2*sqrt(D));(Ebar_j0-2*sqrt(D));length_delta1_Nbar];
    TotalTime(s) = InitalTime(s) + DeltaTime(s) + RuleTime(s) + ReducedQPPTime(s);
end
Time = table(InitalTime,DeltaTime,RuleTime,ReducedQPPTime,TotalTime);
[Acc,Auc,specificity,gmeans,PY]=prediction(kernel_type,kernel_param,traX,traY,tstX,tstY,alpha,nu);




