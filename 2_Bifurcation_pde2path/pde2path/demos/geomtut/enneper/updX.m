function [p,u]=updX(p,u)  % update p.X, return u in p.up, set u to zero
p.u(p.nu+1:end)=u(p.nu+1:end); Xbc=bcX(p,u); p.X(p.idx,:)=Xbc; % update BCs 
N=getN(p,p.X); np=p.nu/p.nc.neq; p.up=u; p.X=p.X+u(1:np).*N; u(1:np)=0; 