function r=sG(p,u) % rhs for AC with Neumann BC, "simple" setting 
par=u(p.nu+1:end); lam=par(1); u=u(1:p.nu); % split u into parameters and PDE variables 
K=p.mat.K; f=lam*u-u.^3; % eff.stiffness and nonlin. 
r=K*u-p.mat.M*f+p.nc.sf*(p.mat.Q*u-p.mat.G); % putting bulk and BCs together 