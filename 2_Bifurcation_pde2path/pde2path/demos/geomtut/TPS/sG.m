function r=sG(p,u)  % PDE rhs 
par=u(p.nu+1:end); H0=par(1); sx=par(5); sy=par(6); sz=par(7); 
Y=p.X; p=updX(p,u); X=p.X; p.X=Y; % X=new scaling, p.X=old 
N0=getN(p,p.X); 
M0=getM(p,X); LB=cotmatrix(X,p.tri); LB=p.mat.fill'*LB; M=p.mat.fill'*M0; 
% Phase conditions 
grX=grad(p.X,p.tri); grXx=grX(1:p.nt,:); 
grXy=grX(p.nt+1:2*p.nt,:);grXz=grX(2*p.nt+1:3*p.nt,:);
E=c2P(p.X,p.tri); grXx=E*grXx; grXy=E*grXy; grXz=E*grXz;
dphix=grXx*p.X; dphix=dot(dphix,N0,2); dphiy=grXy*p.X; dphiy=dot(dphiy,N0,2);
dphiz=grXz*p.X; dphiz=dot(dphiz,N0,2);
N=getN(p,X); 
r=dot(-LB*X+2*H0*(M*N),p.mat.fill'*N,2)+M*(sx*dphiy+sy*dphix+sz*dphiz);