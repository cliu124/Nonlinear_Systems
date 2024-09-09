function r=sG(p,u)  % PDE rhs 
nu=p.nu; np=p.np; id=p.idx; nt=p.nt; 
par=u(nu+1:end); al=par(1); l1=par(2); c0=par(3); srot=par(4); 
u1=u(1:np); u2=u(np+1:2*np); X0=p.X; 
N0=getN(p,X0); X=X0+u1.*N0; M=getM(p,X); M=M(1:np,1:np); 
LB=cotmatrix(X,p.tri); 
N1=getN(p,X); H1=0.5*dot(LB*X,N1,2); % mean curv, cot-Lapl. better than discrete_curv
K=discrete_gaussian_curvature(X,p.tri); K=M\K; 
f1=2*u2.*(u2.^2-K)-2*l1*u2+2*c0*K-2*c0^2*u2; 
r1=LB*u2+M*f1; 
r2=M*u2-H1; 
r=sqrt(X(id,2).^2+X(id,3).^2); r1(id)=-(r-al); % BCs for X 
grX=grad(X,p.tri); E=c2P(X,p.tri); 
grXx=grX(1:p.nt,:); grXx=E*grXx; w=grXx*u1; %w=grXx*X(:,2); 
r2(id)=w(id); % BCs for \pa_x u 
%w2=grXx*(u2+c0); r1(id)=w2(id); % alternate BCs 
%r1(id)=u1(id); 
if p.nc.nq>0; % phase cond, only compute when needed 
 grX=grad(X0,p.tri); E=c2P(X0,p.tri);  grXy=grX(nt+1:2*nt,:); grXz=grX(2*nt+1:3*nt,:);
 grXy=E*grXy; grXz=E*grXz; dphi=-X0(:,3).*grXy*X0+X0(:,2).*grXz*X0; 
 R=X0(:,2).^2+X0(:,3).^2; dphi=dot(dphi./R,N0,2);    r1=r1+srot*dphi;  
end 
r=[r1;r2]; 