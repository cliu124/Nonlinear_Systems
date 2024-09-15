function r=sG(p,u) % Schnakenberg 
par=u(p.nu+1:end); f=nodalf(p,u); L=-p.mat.L; % nonlin and -Lapl 
f(p.bdi)=0; f(p.np+p.bdi)=0; % zero nonlin on bdry for both comp
K=[L,0*L;0*L,par(3)*L]; r=K*u(1:p.nu)-f;  % diffusion matrix and residual
