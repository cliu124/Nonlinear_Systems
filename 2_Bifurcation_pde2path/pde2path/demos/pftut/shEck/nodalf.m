function f=nodalf(p,u) % SH "nonlinearity" for the 2nd-order system formulation
par=u(p.nu+1:end); lam=par(1); nup=par(2); n=p.nu/2; u1=u(1:n); u2=u(n+1:2*n); 
f1=(lam-1)*u1+nup*u1.^2-u1.^3-2*u2; f2=0*u2; % 2nd eqn 0=-lap(u1)+u2 in K 
f=[f1;f2]; 