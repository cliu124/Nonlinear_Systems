function D=makeDxxF(p)
% makeDxxF: FD approximation of \pa_x^2, based on fdcoeffF, here O(h^3)  
% template to set up high order FD differentiation matrices 
% Dx=makeDxxF(p)
x=getpte(p); n=size(x,2); D=sparse(n,n); 
D(1,1:5)=fdcoeffF(2,x(1),x(1:5)); 
D(2,1:5)=fdcoeffF(2,x(2),x(1:5));
for i=3:n-2
     D(i,i-2:i+2) = fdcoeffF(2,x(i),x(i-2:i+2));    
end 
D(n-1,n-4:n)=fdcoeffF(2,x(n-1),x(n-4:n));
D(n,n-4:n)=fdcoeffF(2,x(n),x(n-4:n));