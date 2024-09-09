function r=sG(p,u)  % PDE rhs 
nu=p.nu; np=p.np; par=u(nu+1:end); l1=par(2); l2=par(3); c0=par(4);  
sx=par(7); sy=par(8); sz=par(9); srx=par(10); sry=par(11); srz=par(12);
u1=u(1:np); u2=u(np+1:2*np); X0=p.X; 
N0=getN(p,X0); X=X0+u1.*N0; M=getM(p,X); M=M(1:np,1:np); 
LB=cotmatrix(X,p.tri); 
N1=getN(p,X); H1=0.5*dot(LB*X,N1,2); % mean curv, cot-Lapl. better than discrete_curv
K=discrete_gaussian_curvature(X,p.tri); K=M\K; 
f1=2*u2.*(u2.^2-K)+l1*u2+2*c0*K-2*c0^2*u2-l2; 
r1=LB*u2+M*f1; 
r2=M*u2-H1; 

N=N1; dsx=N(:,1); dsy=N(:,2); dsz=N(:,3); % for PCs 
if p.nc.nq>5 % only compute rot.PC vectors if needed 
grX1dc=cross(p.om.*ones(p.np,3),X,2); uph1=dot(grX1dc,N,2);
grX1dc=cross(p.rh.*ones(p.np,3),X,2); uph2=dot(grX1dc,N,2);
grX1dc=cross(p.w.*ones(p.np,3),X,2); uph3=dot(grX1dc,N,2);
else; uph1=0*dsx; uph2=uph1; uph3=uph1; 
end 
r1=r1+(sx*dsx+sy*dsy+sz*dsz+srx*uph1+sry*uph2+srz*uph3); 
r=[r1;r2]; 