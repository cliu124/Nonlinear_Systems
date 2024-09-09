function A=getA(p,u)
% getA: get areas of triangles in X=p.X+uN
%
% A=getA(p,u)
np=p.nu/p.nc.neq; 
up=u(1:np); u=p.mat.fill*up; N0=getN(p,p.X); X=p.X+u.*N0; 
Af=doublearea(X,p.tri); A=sum(0.5*Af); 
