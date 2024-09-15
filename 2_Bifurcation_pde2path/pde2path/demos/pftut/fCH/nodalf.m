function f=nodalf(p,u)
par=u(p.nu+1:end); u1=u(1:p.np); u2=u(p.np+1:2*p.np); 
eta1=par(1); ga=par(2); eps=par(3); eta2=par(5); etad=eta1-eta2; 
[w,wp,wpp,wppp]=p.fuha.wfu(u1,p); 
f1=wpp.*u2-eps*eta1*u2-eps*etad*wp-eps*ga;
f2=-wp+0*u2;  
f=[f1;f2]; 