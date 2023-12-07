function [X1,Y1,tstX1,tstY1,...
    X2,Y2,tstX2,tstY2,...
    X3,Y3,tstX3,tstY3,...
    X4,Y4,tstX4,tstY4,...
    X5,Y5,tstX5,tstY5]=cross_five(CV_ID,X,Y)
for iteration = 1 : 5
    %%%%先把X,Y按照CV_ID分折
    train_id = CV_ID ~= iteration;
    test_id = ~train_id;
    Y_train= Y(train_id);
    Y_te = Y(test_id);
    X_train = X(train_id,:);
    X_te = X(test_id,:);
    eval(['X',num2str(iteration),'=',' X_train',';']);
    eval(['tstX',num2str(iteration),'=','X_te',';']);
    eval(['Y',num2str(iteration),'=',' Y_train',';']);
    eval(['tstY',num2str(iteration),'=','Y_te',';']);
end