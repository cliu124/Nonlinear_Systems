function r=sG(p,u) % KS in 4th order formulation
par=u(p.nu+1:end); al=par(1); eps=par(3); s=par(4); u=u(1:p.nu); 
K=p.mat.K; M0=p.mat.M0; Kx=p.mat.Kx; uxx=K*u; % sys matrices, and u_xx
r=al*K*uxx-M0*uxx+0.5*M0*(Kx*(u.^2))+s*M0*(Kx*u)+eps; % rhs 