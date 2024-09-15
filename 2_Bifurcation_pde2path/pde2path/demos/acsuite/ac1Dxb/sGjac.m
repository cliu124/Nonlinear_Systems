function Gu=sGjac(p,u)  % jac for AC with x-dependent terms, 2nd version 
par=u(p.nu+1:end); n=p.nu; u=u(1:n); % split u into parameters and PDE vars
x=getpte(p); x=x'; fu=par(2)+0.5*x+3*u.^2-5*par(3)*u.^4; 
Fu=p.mat.M*spdiags(fu,0,n,n); 
Gu=spdiags(1+0.1*x.^2,0,n,n)*p.mat.K-spdiags(0.2*x,0,n,n)*p.mat.Kx-Fu; 