function  [p,u]=updX(p,u)  % update p.X, return u in p.up, set u to zero
par=u(p.nu+1:end); del=par(6); N=getN(p,p.X); np=p.nu/p.nc.neq; p.up=u; 
p.X=p.X+(p.mat.fill*u(1:np)).*N; u(1:np)=0; 
m=min(p.X(p.idx,3)); M=max(p.X(p.idx,3)); old=M-m;
scaz=del/old;  a=[ones(p.np,2) scaz*ones(p.np,1)]; p.X=a.*p.X; 