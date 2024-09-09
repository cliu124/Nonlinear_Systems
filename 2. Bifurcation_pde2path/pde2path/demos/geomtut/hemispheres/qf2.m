function q=qf2(p,u)
par=u(p.nu+1:end);  V0=par(2); V=getV(p,u); q=V-V0;
u=u(1:p.nu); N=getN(p,p.X); 
if p.mpos; M=getM(p); N=M*N; end 
qx=sum(u.*N(:,1)); qy=sum(u.*N(:,2)); 
q=[q; qx; qy];