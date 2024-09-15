function f=nodalf(p,u) % nodal version of nonlinearity
par=u(p.nu+1:end); n=p.np; u1=u(1:n); u2=u(n+1:2*n); 
a=par(1); m=par(2); fs=u2.^m; 
f1=-u1.*fs; f2=-f1; f=[f1;f2];