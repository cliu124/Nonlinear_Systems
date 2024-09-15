function Gu=sGnodpBCjac(p,u) % der of constraints missing! (no matter? check! Alex) 
par=u(p.nu+1:end); u=u(1:p.nu); uf=p.mat.fill*u; % fill 
N0=getN(p,p.X); X=p.X+uf.*N0;
m=min(X(p.idx,3)); M=max(X(p.idx,3)); old=M-m;
del=par(10)/old;  a=[ones(p.np,2) del*ones(p.np,1)]; X=a.*X;
LB=cotmatrix(X,p.tri); LB=p.mat.fill'*LB; M0=getM(p,X); M=p.mat.fill'*M0;
Kr=discrete_gaussian_curvature(X,p.tri); Kr=M0\Kr; Hr=par(1);
F=spdiags((4*Hr.^2-2*Kr),0,p.np,p.np); 
Gu=(-LB-M*F)*p.mat.fill; 