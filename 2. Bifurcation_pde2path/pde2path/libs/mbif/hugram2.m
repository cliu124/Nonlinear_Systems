function V=hugram2(V,W) 
% hugram2: modify V such that V(,i)'*W(,j)=del_ij
m=size(V,2);
for i=1:m
    t=V(:,i);
    for j=1:m; if j~=i; t=t-(t'*W(:,j))*W(:,j); end; end 
    V(:,i)=t/(t'*W(:,i)); 
end