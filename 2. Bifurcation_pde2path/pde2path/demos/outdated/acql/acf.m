function [c,a,f,b]=acf(p,u)  
% separate upde and auxiliary variables, here "par"
par=u(p.nu+1:end); u=u(1:p.nu); al=par(1); del=par(2); ga=par(3); lam=par(4); 
%[ux,uy]=pdegrad(p.mesh.p,p.mesh.t,u); 
% pde for AC 
u=pdeintrp(p.mesh.p,p.mesh.t,u); % interpolate to x-coordinates
c=al+del*u+ga*u.^2; a=0; b=0; % b=-2*par(2)*[u.*ux;u.*uy]; 
f=lam*u+u.^3-u.^5;
