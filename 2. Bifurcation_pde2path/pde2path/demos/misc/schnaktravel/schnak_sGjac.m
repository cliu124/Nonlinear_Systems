function Gu=schnak_sGjac(p,u) % nodal form of Jacobian
par=u(p.nu+1:end); [f1u,f1v,f2u,f2v]=njac(p,u); n=p.nu/p.nc.neq;
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
K=par(3)^2*p.mat.K-par(3)*par(2)*p.mat.Kx; 
Gu=K-p.mat.M*Fu; 
end 

function [f1u,f1v,f2u,f2v]=njac(p,u) % local Jac for Schnakenberg 
u1=u(1:p.nu/p.nc.neq); u2=u(p.nu/p.nc.neq+1:p.nu); 
f1u=-1+2*u1.*u2; 
f1v=u1.^2; 
f2u=-2*u1.*u2; 
f2v=-u1.^2; 
end