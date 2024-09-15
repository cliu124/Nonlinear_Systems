function Gu=sGjac(p,u)
[f1u,f1v,f1w,f2u,f2v,f2w,f3u,f3v,f3w]=njac(p,u); n=p.np;
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n),spdiags(f1w,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n),spdiags(f2w,0,n,n)];
    [spdiags(f3u,0,n,n),spdiags(f3v,0,n,n),spdiags(f3w,0,n,n)]];
Gu=p.mat.K-p.mat.M*Fu; 
end 

function [f1u,f1v,f1w,f2u,f2v,f2w, f3u,f3v,f3w]=njac(p,u) 
u1=u(1:p.np); u2=u(p.np+1:2*p.np); u3=u(2*p.np+1:3*p.np); 
par=u(p.nu+1:end); a=par(1); b=par(2); c=par(3); d=par(4); 
f1u=-(1+b)+2*u1.*u2-c; 
f1v=u1.^2;
f1w=d*ones(p.np,1); 
f2u=b-2*u1.*u2;
f2v=-u1.^2; 
f2w=0*ones(p.np,1); 
f3u=c*ones(p.np,1); 
f3v=0*ones(p.np,1); 
f3w=-d*ones(p.np,1); 
end