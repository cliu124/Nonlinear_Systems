function Gu=sGjac(p,u,r)
% nodal form of Jacobian
par=u(p.nu+1:end); [f1u,f1v,f2u,f2v]=njac(p,u); n=p.nu/p.nc.neq;
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Gu=par(1)^2*p.mat.K-par(2)*p.mat.Kx-p.mat.M*Fu; 
end 

function [f1u,f1v,f2u,f2v]=njac(p,u) 
% Jacobian, nodal version
par=u(p.nu+1:end); u1=u(1:p.nu/p.nc.neq); u2=u(p.nu/p.nc.neq+1:p.nu); 
f1u=1-3*u1.^2;
f1v=-par(1)*(par(4) + par(5)*2*u2 + par(6)*3*u2.^2);
f2u=par(1)^2*ones(p.nu/2,1); 
f2v=-par(1)^2*ones(p.nu/2,1); 
end