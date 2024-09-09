function r=sG(p,u) % rhs for AC with Neumann BC 
par=u(p.nu+1:end); u=u(1:p.nu); % split u into parameters and PDE variables 
f=par(2)*u+u.^3-par(3)*u.^5; % "nonlinearity", i.e., everything but diffusion 
r=par(1)*p.mat.K*u-p.mat.M*f; % the residual 