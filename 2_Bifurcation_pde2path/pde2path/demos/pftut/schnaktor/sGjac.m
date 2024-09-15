function Gu=sGjac(p,u)
Ks=p.mat.K; par=u(p.nu+1:end); % lam,sig,d,R,rho,s; now check if R or rho are active: 
if any(ismember(p.nc.ilam,[4 5])); R=par(4); rho=par(5); Ks=LBtor(p,R,rho); end  
n=p.nu/2; s=par(6);  K=[Ks 0*Ks; 0*Ks par(3)*Ks]; 
[f1u,f1v,f2u,f2v]=njac(p,u); % (nodal) Jacobian of 'nonlin' 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Gu=K-p.mat.M*Fu+s*p.mat.Dphi; 