function f=nodalf(p,u) % SH "nonlinearity" for the 2nd-order system formulation
par=u(p.nu+1:end); lam=par(1); nup=par(2); cp=par(3); 
n=(p.nu-2)/2; u1=u(1:n); u2=u(n+1:2*n); v1=u(p.nu-1); 
if p.qc==1; f1=(lam+cp*v1-1)*u1+nup*u1.^2-u1.^3-2*u2;  % quad-cub
else f1=(lam+cp*v1-1)*u1+nup*u1.^3-u1.^5-2*u2; end
f2=0*u2; % 2nd eqn 0=-lap(u1)+u2 in K 
f=[f1;f2]; 