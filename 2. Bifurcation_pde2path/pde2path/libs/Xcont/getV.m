function V=getV(p,u)
% getV: (algebraic) volume of X=p.X+uN
np=p.nu/p.nc.neq; 
up=u(1:np); u=p.mat.fill*up; N0=getN(p,p.X); 
X=p.X+u.*N0; M=getM(p,X); N=getN(p,X); 
V=abs(sum(M(1:p.np,1:p.np)*(dot(X,N,2)))/3); 
