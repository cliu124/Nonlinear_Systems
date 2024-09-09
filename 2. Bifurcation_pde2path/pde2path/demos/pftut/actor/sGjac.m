function Gu=sGjac(p,u) % AC on torus
par=u(p.nu+1:end); % if R or rho are active, then update LB
if any(ismember(p.nc.ilam,[1 2])); R=par(1); rho=par(2); p.mat.K=LBtor(p,R,rho); end   
lam=par(3); ga=par(4); s=par(5); u=u(1:p.nu); 
fu=lam+2*u-3*ga*u.^2; Fu=spdiags(fu,0,p.nu,p.nu);  
Gu=p.mat.K-p.mat.M*Fu+s*p.mat.Dphi; 