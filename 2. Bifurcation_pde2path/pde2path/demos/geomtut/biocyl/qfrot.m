function q=qfrot(p,u) % rot.PC 
uold=p.up(1:p.np); u=u(1:p.np); X0=p.X; N=getN(p,X0);  
grX=grad(X0,p.tri); nt=p.nt; 
grXy=grX(nt+1:2*nt,:); grXz=grX(2*nt+1:3*nt,:);
E=c2P(X0,p.tri); grXy=E*grXy; grXz=E*grXz; 
dphi=-X0(:,3).*grXy*X0+X0(:,2).*grXz*X0; 
R=X0(:,2).^2+X0(:,3).^2; dphi=dot(dphi./R,N,2); q=dphi'*(u-uold); 
