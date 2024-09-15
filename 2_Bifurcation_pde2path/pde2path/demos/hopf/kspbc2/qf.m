function q=qf(p,u)
par=u(p.nu+1:end); n=p.nu/2; u1=u(1:n); 
M=p.mat.Ms; q1=sum(M*u1)/p.vol-par(2); 
switch p.nc.nq % see if phase cond. is switched on or off
  case 1; q=q1;
  case 2; if isfield(p,'u0x'); u0x=p.u0x(1:n); 
          else u0x=p.mat.Kx*p.u(1:n); end 
          q2=u0x'*u1; q=[q1;q2];
end 