function r=sG(p,u)  % ac3Dpbc, breaking transl.symmetry by adding inhom. in (x,y,z)
par=u(p.nu+1:end); up=u(1:p.nu); u=p.mat.fill*up; 
po=getpte(p); % compute the inhomogeneity (x+1)^2+(y+1)^2+(z+1)^2 
r2=(po(1,:)+1).^2+(po(2,:)+1).^2+(po(3,:)+1).^2; r2=r2'; 
f=par(2)*u+u.^3-par(3)*u.^5+1./(1+r2).*u; 
F=p.mat.M0*f; r=par(1)*p.mat.K*up-F; % bulk part of PDE 