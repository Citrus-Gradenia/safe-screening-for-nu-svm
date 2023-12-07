function [X1,X2, gnd_X]=gencircledata(Center,R,num_data,dev)
phi= 2*pi*rand(1, num_data);
gnd_X=repmat(Center,1,num_data)+R*[cos(phi);sin(phi)];
gnd_Y=ones(1,num_data);

X1=gnd_X+randn(2, num_data)*dev;
X2=repmat(Center,1,num_data)+randn(2, num_data)*dev;