clear
clc
rng(2)
[cfilename]=uigetfile({'*.mat','csv file(*.mat)';'*.*','All files(*.*)'},...
    'Select an mat file to caculate');
STRA='./Artificial data/';
STRB=regexp(cfilename, '.mat', 'split'); 
name =[STRA STRB{1}];
load(name);
outname=STRB{1};
kernel_type='rbf';
scaleX = X;
traX=scaleX;
traY=Y;
tstX =scaleX;
tstY=Y;
[l,n] = size(traX);
nu = 0.1:0.1:0.9;
kernel_param=2^(-3);
tol=1e-8;
[MoreInf,R1_scr,R2_scr,alpha_scr,Acc_scr,Auc_scr,specificity_scr,gmeans_scr,Num_screening,SummaryTime,PY_scr,rho_scr]=main_screening(tol,kernel_type,kernel_param,nu,traX,traY,tstX,tstY);
T_scr=sum(SummaryTime.TotalTime);  
Screening_ratio=Num_screening./Number_SVs;

SummaryTime.OriginTime = T2;
SummaryTime.SpeedUp = SummaryTime.OriginTime./SummaryTime.TotalTime;
SummaryTime.ScreeningRatio = Screening_ratio(:,3);


error_R1=R1-R1_scr<0;
error_R2=R2-R2_scr<0;
Screening_ratio(:,4:5)=[sum(error_R1,1)',sum(error_R2,1)'];
[l,Lc]=size(alpha_nu);
error = sum(sum(abs(alpha_scr-alpha_nu)))/(l*Lc); %screening errors
error_screening = sum(Screening_ratio(:,4)+Screening_ratio(:,5)); 
SummaryTime.ScreeningErrorN = Screening_ratio(:,4)+Screening_ratio(:,5);

Bacc_oc = max(Acc_nu)*100;
Bauc_oc = max(Auc_nu)*100;
MTime_oc = mean(T2);
Bacc_Scr = max(Acc_scr)*100;
Bauc_Scr = max(Auc_scr)*100;
MTime_Scr = mean(SummaryTime.TotalTime);

Speedup = T_oc/T_scr;
MSpeedUp = mean(SummaryTime.SpeedUp);