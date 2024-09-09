function res=cbe(al,d,g,f,beta,m)
% cbe: cubic bifurcation equations, set up and solved in cswibra 
res=zeros(m,1); 
for i=1:m
 r=0; 
 for j=1:m
   r=r+6*beta*(g(i,j)+f(i,j))*al(j); 
    for k=1:m
      for l=1:m
        r=r+d(i,j,k,l)*al(j)*al(k)*al(l); 
      end
    end
 end
 res(i)=r; 
end
