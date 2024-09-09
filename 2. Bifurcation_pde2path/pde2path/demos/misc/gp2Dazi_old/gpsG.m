function r=gpsG(p,u) % GP, speedfem version ... TO DO (23.6.14) !
par=u(p.nu+1:end); lam=par(1); mu=par(2); pa=par(3); c=[1;0;0;1;1;0;0;1]; 
pot=pa*p.mat.poti; x=p.mesh.p(1,:)'; y=p.mesh.p(2,:)'; 
up=p.mat.fill*u(1:p.nu); 
u=u(1:p.np); v=up(p.np+1:2*p.np); ua=u.^2+v.^2; 
% to do !  
f=(-p.mat.pot+mu).*u-ga*u.^3;  f=p.mat.drop*f; 
r=p.mat.K*u-p.mat.M*f; 
