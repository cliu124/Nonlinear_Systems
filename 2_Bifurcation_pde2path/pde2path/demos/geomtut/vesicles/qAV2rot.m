function q=qAV2rot(p,u) % 3 trans + 2 rot.PC  
par=u(p.nu+1:end); uold=0; u=u(1:p.np); 
A=getA(p,u); X0=p.X; N=getN(p,X0); X=p.X; 
%rotational PC's
grX1dc=cross(p.om.*ones(p.np,3),X,2); uph1=dot(grX1dc,N,2);
grX1dc=cross(p.rh.*ones(p.np,3),X,2); uph2=dot(grX1dc,N,2);
%translation PC's
dsx=(u-uold).*N(:,1); dsy=(u-uold).*N(:,2); dsz=(u-uold).*N(:,3);
M=getM(p,X0); M=M(1:p.np,1:p.np); V=getV(p,u);  
q=[A-par(5);
   V-par(6); 
   sum(M*dsx); 
   sum(M*dsy);
   sum(M*dsz);...
   sum(M*uph1.*(u-uold)); 
   sum(M*uph2.*(u-uold))];