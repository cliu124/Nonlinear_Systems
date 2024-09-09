function Gu=sGjac(p,u) 
par=u(p.nu+1:end); a=par(1); m=par(2); s=par(3); sf=p.nc.sf; n=p.nu/p.nc.neq;
u1=u(1:n); u2=u(n+1:2*n); 
f1u=-u2.^m; f1v=-m*u1.*u2.^(m-1); 
f2u=-f1u; f2v=-f1v; 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
Qu=p.mat.Qu; Qv=p.mat.Qv; K=p.mat.K; 
bK=[[a*K 0*K];[0*K K]]; 
Gu=bK-s*p.mat.Kx-p.mat.M(1:p.nu,1:p.nu)*Fu+sf*[[Qu 0*Qu];[0*Qu Qv]]; 
end 

