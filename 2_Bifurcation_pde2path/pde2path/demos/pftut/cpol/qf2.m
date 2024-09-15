function q=qf2(p,u) % phase condition for both, mass and rot invariance 
par=u(p.nu+1:end); R=par(1); m=par(6); u1=u(1:p.nus); u2=u(p.nus+1:p.nu); 
qm=R^2*p.mat.vM1*u1+p.mat.vM2*u2-m; 
uold=p.u(1:p.nus); u0x=p.mat.Dphi*uold; qr=u0x'*u(1:p.nus);
q=[qm;qr]; 