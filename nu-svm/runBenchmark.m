clear
clc
rng(2)
[cfilename]=uigetfile({'*.mat','csv file(*.mat)';'*.*','All files(*.*)'},...
    'Select an mat file to caculate');
STRA='./Benchmark data/';
STRB=regexp(cfilename, '.mat', 'split'); 
name =[STRA STRB{1}];
load(name);
outname=STRB{1};
kernel_type='linear';
traX = X1;
traY = Y1;
tstX = tstX1;
tstY = tstY1;

[l,n] = size(traX);
nu = 0.1:0.1:0.9;
kernel_param=2.^(-3:8); 
tol=1e-6;

i = 1;j = 1;
[R1,R2,alpha_nu,Acc_nu,Auc_nu,specificity,gmeans,Number_SVs,T2,PY_nu,rho_nu]=main_nu(tol,kernel_type,kernel_param(i),nu(j),traX,traY,tstX,tstY);
T_nu=sum(T2);


[R1_scr,R2_scr,alpha_scr,Acc_scr,auc_scr,specificity_scr,gmeans_scr,Num_screening,SummaryTime,PY_scr,rho_scr]=main_screening(tol,kernel_type,kernel_param(i),nu(j),traX,traY,tstX,tstY);
T_scr=sum(SummaryTime.TotalTime);
Screening_ratio=Num_screening./Number_SVs;

SummaryTime.OriginTime = T2;
SummaryTime.SpeedUp = SummaryTime.OriginTime./SummaryTime.TotalTime;
SummaryTime.ScreeningRatio = Screening_ratio(:,3);

error_R1=R1-R1_scr<0;
error_R2=R2-R2_scr<0;
Screening_ratio(:,4:5)=[sum(error_R1,1)',sum(error_R2,1)'];
[l,Lc]=size(alpha_oc);
error=sum(sum(abs(alpha_scr-alpha_nu)))/(l*Lc); 
error_screening=sum(Screening_ratio(:,4)+Screening_ratio(:,5)) ;
SummaryTime.ScreeningErrorN = Screening_ratio(:,4)+Screening_ratio(:,5);

Bacc_nu=max(Acc_nu)*100;
Bauc_nu=max(Auc_nu)*100;
MTime_oc=mean(T2);
Bacc_Scr=max(Acc_scr)*100;
Bauc_Scr=max(auc_scr)*100;
MTime_Scr=mean(SummaryTime.TotalTime);

Speedup=T_nu/T_scr;
MSpeedUp = mean(SummaryTime.SpeedUp);
