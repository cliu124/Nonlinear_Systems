function qu=qrotder(p,u) % u-derivative of rotational PC 
grX=grad(p.X,p.tri); nt=p.nt; X0=p.X; N=getN(p,X0);  
grXy=grX(nt+1:2*nt,:); grXz=grX(2*nt+1:3*nt,:);
E=c2P(p.X,p.tri); grXy=E*grXy; grXz=E*grXz; 
dphi=-X0(:,3).*grXy*X0+X0(:,2).*grXz*X0; dphi=dot(dphi,N,2);
R=X0(:,2).^2+X0(:,3).^2;
qu=[(dphi./R)', 0*dphi']; 