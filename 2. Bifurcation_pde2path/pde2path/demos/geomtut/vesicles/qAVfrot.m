function q=qAVfrot(p,u) % A,V, transl and full rot.PC 
par=u(p.nu+1:end); u=u(1:p.np); A=getA(p,u); X0=p.X; N=getN(p,X0); X=p.X; 
grX1dc=cross(p.om.*ones(p.np,3),X,2); uph1=dot(grX1dc,N,2);
grX1dc=cross(p.rh.*ones(p.np,3),X,2); uph2=dot(grX1dc,N,2);
grX1dc=cross(p.w.*ones(p.np,3),X,2); uph3=dot(grX1dc,N,2);
uold=0; V=getV(p,u); %p.X=p.X+u.*N; V=altvol(p);
dsx=(u-uold).*N(:,1); dsy=(u-uold).*N(:,2); dsz=(u-uold).*N(:,3);

M=getM(p,(X0)); M=M(1:p.np,1:p.np); 
q=[A-p.A0; 
   V-par(6); 
   sum(M*dsx); 
   sum(M*dsy);
   sum(M*dsz);...
   sum(M*uph1.*(u-uold)); 
   sum(M*uph2.*(u-uold));
   sum(M*uph3.*(u-uold))];
