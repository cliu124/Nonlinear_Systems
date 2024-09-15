function qu=qjac(p,~)
%uo=p.u(1:p.nu);  Kx=p.mat.Kx; 
M=p.mat.Ms; n=p.nu/2; 
q1u=([M*ones(length(M),1)/p.vol; 0*ones(length(M),1)])';
switch p.nc.nq % Decide whether phase cond. is switched on or off
  case 1; qu=q1u;
  case 2; if isfield(p,'u0x'); u0x=p.u0x(1:n); 
    else u0x=p.mat.Kx*p.u(1:n); end 
    q2u=[u0x', 0*u0x']; qu=[q1u;q2u];
end