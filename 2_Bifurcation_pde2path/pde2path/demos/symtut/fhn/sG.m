function r=sG(p,u) % nodal form of residual
par=u(p.nu+1:end); 
r=(par(1)^2*p.mat.K-par(2)*p.mat.Kx)*u(1:p.nu)-p.mat.M*nodalf(p,u); 
end