function [E1, E2]=chE(p,u)  % energy for CH  
par=u(p.nu+1:end); eps=par(2); u=u(1:p.np); 
ux=p.mat.Dx*u; uy=p.mat.Dy*u; uz=p.mat.Dz*u; % for |grad u|^2
W=0.25*(u.^2-1).^2; ux2=ux.^2+uy.^2+uz.^2; sig=sqrt(2)/3; 
dens1=0.5*eps^2*ux2+W-0.25; E1=sum(p.mat.M*dens1); % physical E 
dens2=0.5*eps*ux2+W/eps; E2=sum(p.mat.M*dens2)/(2*sig); % interface-length E