function [f1u,f1v,f2u,f2v]=nodaljac(p,u) % Jacobian for Schnakenberg
par=u(p.nu+1:end); m=par(4); ga=1/m; 
u1=u(1:p.nu/2); u2=u(p.nu/2+1:p.nu); del=p.nc.del; 
iz=find(u1<2*p.nc.del); 
fsd=ga*u1.^(ga-1); 
fsd(iz)=((u1(iz)+del).^ga-u1(iz).^ga)/del; % singular derivative, FD-like 
f1u=-fsd+2*u1.*u2; 
f1v=u1.^2; f2u=-2*u1.*u2; f2v=-f1v;  % regular part 