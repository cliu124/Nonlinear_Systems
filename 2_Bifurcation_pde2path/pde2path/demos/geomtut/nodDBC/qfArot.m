function q=qfArot(p,u)
par=u(p.nu+1:end); A0=par(3); A=getA(p,u); q=A0-A;
u=u(1:p.np); N=getN(p,p.X);  
grX=grad(p.X,p.tri); %grX=grad(p.Xold,p.tri); 
grXx=grX(1:p.nt,:); grXy=grX(p.nt+1:2*p.nt,:);
E=c2P(p.X,p.tri); grXx=E*grXx; grXy=E*grXy; 
dphix=-p.X(:,2).*grXx*p.X+p.X(:,1).*grXy*p.X; dphix=dot(dphix,N,2);
q3=dphix'*u; q=[q;q3];