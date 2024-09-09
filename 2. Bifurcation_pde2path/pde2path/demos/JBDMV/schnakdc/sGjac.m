function Gu=sGjac(p,u)
par=u(p.nu+1:end); n=p.np;
[f1u,f1v,f2u,f2v]=nodaljac(p,u); % (nodal) Jacobian of 'nonlin' 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Ks=p.mat.K; c=par(2); d=par(3); 
Gu=[c*Ks 0*Ks; 0*Ks d*Ks]-p.mat.M*Fu; 