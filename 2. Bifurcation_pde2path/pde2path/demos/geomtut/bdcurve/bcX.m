function Xbc=bcX(p,u) % set BCs; called in sGbdcurv, with u=0 on bdry there 
par=u(p.nu+1:end); al=par(4); k=par(5); phi=angle(p.X(p.idx,1)+1i*p.X(p.idx,2));
b=sqrt(1-(al*cos(k*phi).^2));  Xbc=[b.*cos(phi), b.*sin(phi), al*cos(k*phi)]; 