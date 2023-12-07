function [X,Y]=genExclusiveData(Center,num_data,dev)
P1=[-Center;Center];
P2=[Center;-Center];
N1=[Center;Center];
N2=[-Center;-Center];

gnd_X1p=repmat(P1,1,num_data*.25);
gnd_X2p=repmat(P2,1,num_data*.25);
gnd_X1n=repmat(N1,1,num_data*.25);
gnd_X2n=repmat(N2,1,num_data*.25);
gnd_X=[gnd_X1p,gnd_X2p,gnd_X1n,gnd_X2n];
X=gnd_X+randn(2, num_data)*dev;
gnd_Y=ones(1,num_data*0.5);
Y=[gnd_Y,-gnd_Y];
X=X';
Y=Y';

