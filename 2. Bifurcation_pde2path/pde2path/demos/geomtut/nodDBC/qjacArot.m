function qu=qjacArot(p,u) % \pa_u A
par=u(p.np+1:end); H0=par(1); u=u(1:p.np); N0=getN(p,p.X); 
X=p.X+u.*N0; M=getM(p,X);% N=getN(p,X); 
qu1=-2*H0*sum(M,1); % flipped cause sign of N flipped 
grX=grad(p.X,p.tri); grXx=grX(1:p.nt,:); grXy=grX(p.nt+1:2*p.nt,:); 
E=c2P(p.X,p.tri); grXx=E*grXx; grXy=E*grXy; 
dphiy=-p.X(:,2).*grXx*p.X+p.X(:,1).*grXy*p.X; dphiy=dot(dphiy,N0,2);
qu=[qu1;dphiy'];