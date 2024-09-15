function Gu=sGjac(p,u)
[f1u,f1v,f2u,f2v]=njac(p,u); n=p.np; del=u(p.nu+3); 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Gu=del*p.mat.K-p.mat.M*Fu; 
s=u(p.nu+4);  Gu=Gu+s*p.mat.Krot; 
end 

function [f1u,f1v,f2u,f2v]=njac(p,u) 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); al=u(p.nu+2); 
r=u(p.nu+1);  usq=u1.^2+u2.^2; 
f1u=0.5+r-3*u1.^2-u2.^2+2*al*u1.*u2;  
f1v=1-2*u2.*(u1-al*u2)+al*usq; 
f2u=-1-2*u1.*(u2+al*u1)-al*usq; 
f2v=r-3*u2.^2-u1.^2-2*al*u1.*u2; 
end