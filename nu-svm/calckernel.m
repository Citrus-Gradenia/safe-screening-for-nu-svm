function K = calckernel(kernel_type,kernel_param,X2,X1)
n1=size(X1,1);
if nargin>2
    n2=size(X2,1);
end

switch kernel_type

    case 'linear'
        if nargin>2
            K=X2*X1';
        else
            K=X1*X1';
        end

    case 'poly'
        if nargin>2
            K=(X2*X1').^kernel_param;
        else
            K=(X1*X1').^kernel_param;
        end

    case 'rbf'
        if nargin>2
            K = exp(-(repmat(sum(X1.*X1,2)',n2,1) + ...
                repmat(sum(X2.*X2,2),1,n1) - 2*X2*X1') ...
                /(2*kernel_param^2));
        else
            P=sum(X1.*X1,2);
            K = exp(-(repmat(P',n1,1) + repmat(P,1,n1) ...
                - 2*X2*X1')/(2*kernel_param^2));
        end

    otherwise
        error('Unknown kernel function.');
end

