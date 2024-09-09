function q=qf(p,u) % mass (and phase) constraint for KS 
par=u(p.nu+1:end); u=u(1:p.nu); % extract pars and u-vars 
q=sum(p.mat.M0*u)/p.vol-par(2); % mass constraint 
if p.nc.nq==2; % if active, then add phase constraint
  if isfield(p,'u0x'); u0x=p.u0x(1:p.nu); else u0x=p.mat.Kx*p.u(1:p.nu); end
  q=[q;u0x'*u];
end 