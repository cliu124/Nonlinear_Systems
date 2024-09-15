function r=sG(p,u)  % ac2D, with optional "classical FEM" evaluation 
if p.cl; r=res_cl(p,u); return; end; % to compare to classical FEM 
par=u(p.nu+1:end); u=u(1:p.nu); f=par(2)*u+u.^3-par(3)*u.^5; 
r=par(1)*p.mat.K*u-p.mat.M*f+p.nc.sf*(p.mat.Q*u-par(4)*p.mat.G); 