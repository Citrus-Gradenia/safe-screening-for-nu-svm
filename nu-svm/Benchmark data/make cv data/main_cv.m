Kfold = 5;

CV_ID = cv_gen(Y,Kfold);
[X1,Y1,tstX1,tstY1,...
    X2,Y2,tstX2,tstY2,...
    X3,Y3,tstX3,tstY3,...
    X4,Y4,tstX4,tstY4,...
    X5,Y5,tstX5,tstY5] = cross_five(CV_ID,X,Y);