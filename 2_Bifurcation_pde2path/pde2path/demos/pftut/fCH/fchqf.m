function q=fchqf(p,u)
par=u(p.nu+1:end); u=u(1:p.np); 
q=sum(p.mat.Ms*u)/p.vol-par(4);
%q=sum(p.mat.Ms*u)-par(4);
