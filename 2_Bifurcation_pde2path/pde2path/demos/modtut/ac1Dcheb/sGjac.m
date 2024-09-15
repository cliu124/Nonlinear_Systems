function Gu=sGjac(p,u) % Jacobian 
n=p.nu; par=u(n+1:end); u=u(1:n); c=par(1); lam=par(2); c2=par(3); c3=par(4); 
fu=lam+2*c2*u+3*c3*u.^2; fu(1)=0; fu(n)=0; % zero f_u on bdry 
Fu=spdiags(fu,0,n,n); % local Jac, converted to matrix 
Gu=-c*p.mat.L-Fu; % build Jac from bulk 