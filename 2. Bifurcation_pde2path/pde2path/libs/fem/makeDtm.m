function Dt=makeDtm(t)
% makeDtm: dt-matrix 
n=length(t); Dt=sparse(n,n); 
Dt(1,1:5)=fdcoeffF(1,t(1),t(1:5)); 
Dt(2,1:5)=fdcoeffF(1,t(2),t(1:5));
for i=3:n-2
     Dt(i,i-2:i+2) = fdcoeffF(1,t(i),t(i-2:i+2));    
end 
Dt(n-1,n-4:n)=fdcoeffF(1,t(n-1),t(n-4:n));
Dt(n,n-4:n)=fdcoeffF(1,t(n),t(n-4:n));