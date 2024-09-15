function Gu=sGjac(p,u)
par=u(p.nu+1:end); n=p.np;
[f1u,f1v,f2u,f2v]=njac(p,u); % (nodal) Jacobian of 'nonlin' 
f1u(p.bdi)=0; f1v(p.bdi)=0; f2u(p.bdi)=0; f2v(p.bdi)=0; % zero fjac on bdry 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
L=-p.mat.L; K=[L,0*L;0*L,par(3)*L]; 
Gu=K-Fu; 