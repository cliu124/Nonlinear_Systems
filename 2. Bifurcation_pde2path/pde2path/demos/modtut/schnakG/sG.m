function r=sG(p,u)  % (graph) PDE rhs 
f=nodalf(p,u); % compute "nonlinearity" (everything but diffusion) 
par=u(p.nu+1:end); d1=par(3); d2=par(4); u=u(1:p.nu); % split in par and PDE u 
L=p.mat.L; K=[d1*L 0*L; 0*L d2*L];  % compose 2-compo Lapl. from scalar one
r=K*u-f; % the residual 