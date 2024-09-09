function r=sG(p,u) % Schnakenberg 
f=nodalf(p,u); par=u(p.nu+1:end); % lam,sig,d 
L=-p.mat.L; K=[L,0*L;0*L,par(3)*L]; % compose diffusion matrix from scalar L 
r=K*u(1:p.nu)-p.mat.M*f; 