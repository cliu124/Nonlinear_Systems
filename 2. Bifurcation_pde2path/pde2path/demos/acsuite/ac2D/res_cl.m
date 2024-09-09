function r=res_cl(p,u) % compute res. of a OOPDE soln in classical way 
par=u(p.nu+1:end);u=u(1:p.np); lam=par(2);ga=par(3); % separate u and pars 
ut=p.pdeo.grid.point2Center(u); % u interpolated to element (triangle) centers
f=lam*ut+ut.^3-ga*ut.^5; % evaluate f on centers 
[K,~,F]=p.pdeo.fem.assema(p.pdeo.grid,par(1),0,f); % assemble K and F (K not 
% strictly needed since c=par(1) is const and hence we could use par(1)*p.mat.K) 
r=K*u+p.nc.sf*(p.mat.Q*u-par(4)*p.mat.G)-F; % the rhs 