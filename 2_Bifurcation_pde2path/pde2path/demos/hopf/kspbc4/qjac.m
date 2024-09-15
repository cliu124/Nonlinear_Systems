function qu=qjac(p,~)
M=p.mat.M0; q1u=(M*ones(length(M),1)/p.vol)';
switch p.nc.nq % see if phase cond. is switched on or off
  case 1; qu=q1u;
  case 2; 
    if isfield(p,'u0x'); u0x=p.u0x(1:p.nu); else u0x=p.mat.Kx*p.u(1:p.nu); end
    q2u=u0x'; qu=[q1u;q2u];
end