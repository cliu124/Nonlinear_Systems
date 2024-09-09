function r=sG(p,u) % Schnakenberg 
par=u(p.nu+1:end); lam=par(1); c=par(2); d=par(3); m=par(4); 
np=p.np; u1=u(1:np); u2=u(np+1:2*np); if p.zc2; u1=max(u1,0); end 
uo=u1.^(1/m); f1=-uo+u1.^2.*u2; f2=lam-u1.^2.*u2; f=[f1; f2]; 
K1= p.mat.K; K=[c*K1 0*K1; 0*K1 d*K1];  r=K*u(1:p.nu)-p.mat.M*f; 