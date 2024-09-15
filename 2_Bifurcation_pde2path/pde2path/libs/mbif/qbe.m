function res=qbe(al,a,b,c,al0,m)
% qbe: quadratic bif eqns after Keller77, set up and called in qswibra 
res=zeros(m,1); 
for i=1:m
   r=c(i)*al0^2; 
   for j=1:m
       r=r+2*b(i,j)*al(j)*al0; 
       for k=1:m
         r=r+a(i,j,k)*al(j)*al(k); 
       end
   end
   res(i)=r; 
end