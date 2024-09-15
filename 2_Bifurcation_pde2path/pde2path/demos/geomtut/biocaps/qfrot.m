function q=qfrot(p,u) % rot.PC 
uold=p.up(1:p.np); u=u(1:p.np); X0=p.X; N=getN(p,X0);  
grX=grad(X0,p.tri); nt=p.nt; 
grXy=grX(nt+1:2*nt,:); grXx=grX(1:nt,:);
E=c2P(X0,p.tri); grXy=E*grXy; grXx=E*grXx; 
dphi=-X0(:,2).*grXx*X0+X0(:,1).*grXy*X0; 
R=X0(:,1).^2+X0(:,2).^2; R=1; dphi=dot(dphi./R,N,2); q=dphi'*u; 
