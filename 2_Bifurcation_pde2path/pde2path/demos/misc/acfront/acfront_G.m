function [c,a,f,b]=acfront_G(p,u)  
% separate pde and auxiliary variables, here "par"
par=u(p.nu+1:end);  
u=pdeintrp(p.mesh.p,p.mesh.t,p.mat.fill*u(1:p.nu)); % interpolate to x-coordinates

c=1; a=0; b=[0;par(3)]; f=par(1)*u.*(1-u).*(par(2)+u);
