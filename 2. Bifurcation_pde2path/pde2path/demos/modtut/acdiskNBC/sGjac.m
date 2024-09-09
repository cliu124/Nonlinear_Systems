function Gu=sGjac(p,u) % AC 
n=p.nu; par=u(n+1:end); u=u(1:n); c=par(1); lam=par(2); c2=par(3); c3=par(4); 
fu=lam+2*c2*u+3*c3*u.^2; fu(1:p.na)=0; % zero out at bdry 
Fu=spdiags(fu,0,n,n); % loc. Jac and convert to matrix 
Gu=-c*p.mat.L-Fu+par(5)*p.mat.Dphi; % build Jac 