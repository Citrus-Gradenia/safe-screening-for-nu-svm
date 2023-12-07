%-- New data with new response --%
function CV_ID = cv_gen(Y,Kfold)

CV_ID =zeros(length(Y),1);

num_pos = nnz(Y==1);
num_neg = nnz(Y==-1);

cv_ind_pos = crossvalind('Kfold',num_pos,Kfold);
cv_ind_neg = crossvalind('Kfold',num_neg,Kfold);

idx_pos = find(Y==1);
idx_neg = find(Y==-1);

for k = 1:Kfold
   CV_ID(idx_pos(cv_ind_pos==k)) = k; 
   CV_ID(idx_neg(cv_ind_neg==k)) = k; 
end

 