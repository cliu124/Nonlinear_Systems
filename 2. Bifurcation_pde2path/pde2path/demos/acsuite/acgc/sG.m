function r=sG(p,u)  % pde for AC with global coupling of type c(u)*<h(u)>
fnonloc=fnl(p,u); par=u(p.nu+1:end); u=u(1:p.nu); 
K=par(1)*p.mat.K; f=par(2)*u+u.^3-par(3)*u.^5+fnonloc; 
r=K*u-p.mat.M*f+p.nc.sf*p.mat.Q*u; 