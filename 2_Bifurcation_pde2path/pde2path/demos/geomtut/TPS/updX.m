function [p,u]=updX(p,u)
np=p.nu/p.nc.neq; par=u(p.nu+1:end); del=par(4); p.up=u; u=u(1:np); X=p.X; 
m=min(X(:,3)); M=max(X(:,3)); old=M-m; scaz=del/old;  
%a=[1/sqrt(scaz)*ones(p.np,1) 1/sqrt(scaz)*ones(p.np,1) scaz*ones(p.np,1)]; X=a.*X;
X(:,3)=scaz*X(:,3); %a=scaz*ones(p.np,1); % only scale z 
N0=getN(p,X);  X=p.mat.drop*X; N0=p.mat.drop*N0; X=X+u.*N0; 
p.X=[p.Xfillx*X(:,1) p.Xfilly*X(:,2) p.Xfillz*X(:,3)]; u(1:np)=0; 
par(3)=getA(p,u);  u=[u;par];
