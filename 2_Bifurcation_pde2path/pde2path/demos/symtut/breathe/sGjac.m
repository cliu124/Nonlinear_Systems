function Gu=sGjac(p,u) 
par=u(p.nu+1:end); D=par(5); s=par(6); n=p.nu/p.nc.neq;
al=par(1); bet=par(2); ga=par(3); del=par(4); 
u1=u(1:n); u2=u(n+1:2*n); ov=ones(n,1);
f1u=-3*u1.^2+2*(al+bet)*u1-al*bet;   f1v=-ov; 
f2u=del*ov; f2v=-del*ga*ov; 
Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
    [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
K=p.mat.K; 
bK=[[K 0*K];[0*K D*K]]; 
Gu=bK-s*p.mat.Kx-p.mat.M(1:p.nu,1:p.nu)*Fu; 
end