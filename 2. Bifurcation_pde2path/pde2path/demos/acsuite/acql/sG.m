function r=sG(p,u)  % rhs for ql-AC 
par=u(p.nu+1:end); lam=par(2); ga=par(3); 
u=u(1:p.nu); gr=p.pdeo.grid; ut=(p.mat.p2c*u)'; % interpolate to elem. centers 
c=p.fuha.cfu(ut,par); f=lam*ut+ut.^3-ga*ut.^5; % coefficients for assembly
[K,~,F]=p.pdeo.fem.assema(gr,c,0,f); % assemble K and F (M not used) 
r=K*u-F+p.nc.sf*(p.mat.Q*u-p.mat.G);  