function Gu=sGjac(p,u)
par=u(p.nu+1:end); R=par(4); Ks=p.mat.K/(R^2); 
n=p.nu/2; s=par(5); K=[Ks 0*Ks; 0*Ks par(3)*Ks]; 
[f1u,f1v,f2u,f2v]=njac(p,u); % (nodal) Jacobian of 'nonlin' 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Gu=K-p.mat.M*Fu+s*p.mat.Dphi; 