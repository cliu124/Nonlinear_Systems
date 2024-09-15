function [cj,aj,bj]=acfront_Gjac(p,u)  
% jacobian for AC 
% separate pde and auxiliary variables, here "par"
par=u(p.nu+1:length(u));
u=pdeintrp(p.mesh.p,p.mesh.t,p.mat.fill*u(1:p.nu));  % interpolate to x-coordinates

cj=1; bj=[0;par(3)];                                      % diffusion and convection are linear in u
fu=par(1)*par(2)-3*par(1)*u.^2 + 2*par(1)*(1-par(2))*u;   % u-derivative of nonlinearity
aj=-fu;                                          % matrix part for pde
