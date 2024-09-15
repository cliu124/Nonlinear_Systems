function Dx=makeDxF(p)
% makeDxF: FD approximation of \pa_x, based on fdcoeffF, here O(h^3)  
% template to set up high order FD differentiation matrices 
% Dx=makeDxF(p)
x=getpte(p); n=size(x,2); Dx=sparse(n,n); 
Dx(1,1:5)=fdcoeffF(1,x(1),x(1:5)); 
Dx(2,1:5)=fdcoeffF(1,x(2),x(1:5));
for i=3:n-2
     Dx(i,i-2:i+2) = fdcoeffF(1,x(i),x(i-2:i+2));    
end 
Dx(n-1,n-4:n)=fdcoeffF(1,x(n-1),x(n-4:n));
Dx(n,n-4:n)=fdcoeffF(1,x(n),x(n-4:n));