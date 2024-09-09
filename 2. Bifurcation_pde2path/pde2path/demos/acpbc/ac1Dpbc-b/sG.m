function r=sG(p,u)  % AC with periodic BC 
par=u(p.nu+1:end); up=u(1:p.nu); % params, and u on periodic domain 
u=p.mat.fill*up; f=par(2)*u+u.^3-par(3)*u.^5; % fill, then comp nonlin. 
F=p.mat.M0*f; % multiply by M, map back to active nodes of periodic domain 
r=p.mat.K*up+par(4)*p.mat.Kx*up-F;  % resi, with par(4)=wave-speed