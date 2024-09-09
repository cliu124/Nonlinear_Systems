function qu=qrotder(p,u) % u-derivative of rotational PC 
grX=grad(p.X,p.tri); nt=p.nt; X0=p.X; N=getN(p,X0);  
grXy=grX(nt+1:2*nt,:); grXx=grX(1:nt,:);
E=c2P(X0,p.tri); grXy=E*grXy; grXx=E*grXx; 
dphi=-X0(:,2).*grXx*X0+X0(:,1).*grXy*X0; 
R=X0(:,1).^2+X0(:,2).^2; R=1; dphi=dot(dphi./R,N,2);
qu=[(dphi./R)', 0*dphi']; 