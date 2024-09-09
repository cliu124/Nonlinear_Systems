function Gu=sGjac(p,u)
[f1u,f1v,f2u,f2v]=njac(p,u); n=p.np;
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Gu=p.mat.K-p.mat.M*Fu; 
end 
function [f1u,f1v,f2u,f2v]=njac(p,u) 
% Jacobian for Schnakenberg, nodal version
u1=u(1:p.np); u2=u(p.np+1:2*p.np); par=u(p.nu+1:end);  
f1u=-1+2*u1.*u2;
f1v=u1.^2; 
f2u=-2*u1.*u2; 
f2v=-u1.^2; 
end