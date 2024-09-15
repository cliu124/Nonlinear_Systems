function q=qAV(p,u) % A,V, and 3 trans .PC 
par=u(p.nu+1:end); uold=p.up(1:p.np); u=u(1:p.np); 
A=getA(p,u);  V=getV(p,u); X0=p.X; N=getN(p,X0); 
% translational PC's
dsx=(((u-uold).*N(:,1))')'; dsy=(((u-uold).*N(:,2))')'; dsz=(u-uold).*N(:,3);
M=getM(p,X0); M=M(1:p.np,1:p.np); 
q=[A-par(5); 
   V-par(6); 
   sum(M*dsx); 
   sum(M*dsy);
   sum(M*dsz)]; 