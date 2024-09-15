function [p,u]=updX(p,u)  % local mod of updX.m, with update of BCs. 
if p.bcsw==2; Xbc=bcX(p,u); p.X(p.idx,:)=Xbc; end 
N=getN(p,p.X); np=p.nu/p.nc.neq; p.up=u; p.X=p.X+u(1:np).*N; u(1:np)=0; 