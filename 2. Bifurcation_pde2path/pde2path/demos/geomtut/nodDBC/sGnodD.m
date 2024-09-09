function r=sGnodD(p,u)  % PDE rhs 
par=u(p.nu+1:end); H=par(1); u=u(1:p.np); N0=getN(p,p.X); 
X=p.X+u.*N0; M=getM(p,X); LB=cotmatrix(X,p.tri); N=getN(p,X); 
grX=grad(p.X,p.tri); grXx=grX(1:p.nt,:); grXy=grX(p.nt+1:2*p.nt,:);
E=c2P(p.X,p.tri); grXx=E*grXx; grXy=E*grXy; % interpol.surf grad to nodes 
dphix=-p.X(:,2).*grXx*p.X+p.X(:,1).*grXy*p.X; dphix=dot(dphix,N0,2);
r=dot(-LB*X-M*(2*H*N),N,2)+par(6)*dphix;  % rhs, incl.rotational PC
r(p.idx)=u(p.idx); % DBCs