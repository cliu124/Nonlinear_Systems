function Gu=sGjac(p,u) % matrix free Jac, except when called with p.needGu=1
global p2pglob; F=p2pglob.F; 
n=p.nu; par=u(n+1:end); uf=u(1:n); lam=par(1); c2=par(2); c3=par(3); 
u=F'*uf; fu=lam+2*c2*u+3*c3*u.^2; p2pglob.fu=fu; 
try needGu=p.needGu; catch needGu=0; end  
if needGu; t1=tic; Fu=F*(spdiags(fu,0,n,n)*F'); t1=toc(t1); fprintf('time for Gu: %g\n',t1);  
   Gu=diag(p2pglob.mu)-Fu; % provide Gu for q(c)swibra 
else Gu=speye(n);  end     % provide dummy 