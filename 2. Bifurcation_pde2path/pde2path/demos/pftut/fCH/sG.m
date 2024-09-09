function r=sG(p,u)  % pearling for fCH, nodal version 
% -eps^2 Del u+W'(u)+v=0
% -eps^2 Del v+W''(u)v-eps*eta1*v-eps*etad*W'(u)-eps*ga=0 
%p.nu, p.np, p.mat, pause 
f=nodalf(p,u); par=u(p.nu+1:end); eps2=par(3)^2; 
Ks=p.mat.Ks; Ms=p.mat.Ms; 
Mnl=[[Ms 0*Ms];[0*Ms Ms]];
Ksys=[[par(6)*p.mat.Kx -eps2*Ks]; [eps2*Ks Ms]]; 
r=Ksys*u(1:p.nu)-Mnl*f; 