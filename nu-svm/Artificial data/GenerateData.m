mu1 = [2,2];
mu2 = [-2,-2];
sigma = [1,0;0,1]; 
n = 50;
rng (1);
X1 = mvnrnd(mu1,sigma,n);
X2 = mvnrnd(mu2,sigma,n);
plot(X1(:,1),X1(:,2),'bo')
hold on
plot(X2(:,1),X2(:,2),'ro')
X = [X1;X2];
Y = [ones(n,1);-ones(n,1)];
save SimpleData_example X Y