function Gu=sGjac(p,u) % Jacobian 
n=p.nu; par=u(n+1:end); u=u(1:n); c=par(1); lam=par(2); c2=par(3); c3=par(4); 
fu=lam+2*c2*u+3*c3*u.^2; Fu=spdiags(fu,0,n,n); % local Jac, converted to matrix 
Gu=-c*p.mat.D2(2:n+1,2:n+1)-Fu; % build Jac from bulk 
Gu(:,1)=Gu(:,1)-c*p.mat.D2(2:n+1,1); % add contrib. of boundary points
Gu(:,n)=Gu(:,n)-c*p.mat.D2(2:n+1,n+2); 