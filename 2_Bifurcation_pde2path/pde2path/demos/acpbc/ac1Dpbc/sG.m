function r=sG(p,u)  % AC with periodic BC 
par=u(p.nu+1:end); up=u(1:p.nu); % params, and u on periodic domain 
u=p.mat.fill*up; % extend ('fill') u to full domain 
x=getpte(p); x=x'; % extract the point coordinates from p 
f=par(2)*u+u.^3-par(3)*u.^5+0.5*x.*u; % f, with x-dependent term
F=p.mat.M0*f; % multiply by M, map back to active nodes of periodic domain 
r=p.mat.K*up-F;  % bulk part of PDE 
  