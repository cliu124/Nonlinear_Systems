function V=hugram(V)
% hugram: Gram-Schmidt for q(c)swibra
m=size(V,2);
V(:,1)=V(:,1)/norm(V(:,1));
for i=2:m
    t=V(:,i);
    for j=1:i-1; t=t-(V(:,i)'*V(:,j))*V(:,j);end
    V(:,i)=t/norm(t); 
end