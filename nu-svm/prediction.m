function [accuracy,auc,specificity,gmeans,PY]=prediction(kernel_type,kernel_param,traX,traY,tstX,tstY,alpha,nu)
clc
[l,~]=size(traX);
[lt,~]=size(tstX);
pos=(traY==1);
neg=(traY==-1);
tstX=[tstX,ones(lt,1)];
Lengthnu=length(nu);
PY = zeros(lt,Lengthnu);
accuracy=zeros(Lengthnu,1);
auc=zeros(Lengthnu,1);
Ktst=calckernel(kernel_type,kernel_param,tstX,traX);
for s=1:Lengthnu  
    PY0=sign(Ktst*diag(traY)*alpha(:,s)); 
    PY0(PY0==0)=1;
    PY(:,s)=PY0;
    tp = nnz(tstY == PY(:,s) & PY(:,s) == 1);
    fp = nnz(tstY ~= PY(:,s) & PY(:,s) == 1);
    tn = nnz(tstY == PY(:,s) & PY(:,s) == -1);
    fn = nnz(tstY ~= PY(:,s) & PY(:,s) == -1);
    total = tp+fp+tn+fn;
    accuracy(s)  = (tp + tn)/total;
    precision(s) = tp / (tp + fp);
    recall(s)    = tp / (tp + fn);
    sensitivity(s) = recall(s);
    specificity(s) = tn / (tn + fp);
    gmeans(s)=sqrt(sensitivity(s)*specificity(s));
    [~,I]=sort(PY(:,s));
    M=0;N=0;
    for i=1:length(PY(:,s))
        if(tstY(i)==1)
            M=M+1;
        else
            N=N+1;
        end
    end
    m=0;
    for i=M+N:-1:1
        if(tstY(I(i))==1)
            m=m+i;
        end
    end
    auc(s)=(m-(M+1)*M/2)/(M*N);
end




