function f=nodalf(p,u) % nodal version of nonlinearity
par=u(p.nu+1:end); n=p.nu/2; u1=u(1:n); u2=u(n+1:2*n); 
al=par(1); bet=par(2); ga=par(3); del=par(4); 
f1=u1.*(u1-al).*(bet-u1)-u2; 
f2=del*(u1-ga*u2); 
f=[f1;f2];