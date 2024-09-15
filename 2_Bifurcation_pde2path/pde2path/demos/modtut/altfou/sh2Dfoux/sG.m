function r=sG(p,u)  % PDE rhs 
f=nodalf(p,u); % compute "nonlinearity" (everything but diffusion) 
r=p.mat.L*u(1:p.nu)-f;    % bulk part of PDE 
  