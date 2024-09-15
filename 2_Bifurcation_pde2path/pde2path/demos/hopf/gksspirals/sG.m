function r=sG(p,u) % compute pde-part of residual
f=nodalf(p,u); del=u(p.nu+3); % s=del (for cont. in domain size) 
r=del*p.mat.K*u(1:p.nu)-p.mat.M*f; 
s=u(p.nu+4); r=r+s*p.mat.Krot*u(1:p.nu);  
    