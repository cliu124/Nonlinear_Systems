function r=sG(p,u)  % PDE rhs 
nu=p.nu; np=p.np; id=p.idx; par=u(nu+1:end); al=par(1); l1=par(2); c0=par(3); b=par(4);
u1=u(1:np); u2=u(np+1:2*np); X0=p.X; 
N0=getN(p,X0); X=X0+u1.*N0; M=getM(p,X); M=M(1:np,1:np); LB=cotmatrix(X,p.tri); 
N1=getN(p,X); H1=0.5*dot(LB*X,N1,2); % 
K=discrete_gaussian_curvature(X,p.tri); K=M\K; 
f1=2*u2.*(u2.^2-K)-2*l1*u2+2*c0*K-2*c0^2*u2; 
r1=LB*u2+M*f1; 
r2=M*u2-H1; 
if p.nc.nq>0; % phase cond, only compute when needed 
 grX=grad(X0,p.tri); nt=p.nt; srot=par(5); 
 grXy=grX(nt+1:2*nt,:); grXx=grX(1:nt,:);
 E=c2P(X0,p.tri); grXy=E*grXy; grXx=E*grXx; 
 dphi=-X0(:,2).*grXx*X0+X0(:,1).*grXy*X0; 
%R=X0(:,1).^2+X0(:,2).^2; 
 R=1; dphi=dot(dphi./R,N0,2); r1=r1+srot*dphi;  
end 
r1(id)=u1(id); % BCs for X 
kapn=-dot(N1(id,:),X(id,:),2)/al^2; % kap_N for circle in x-y-plane 
%size(kapn)
r2(id)=u2(id)-c0+b*kapn; 
r=[r1;r2]; 