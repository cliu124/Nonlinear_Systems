function out=hpcbra(p,u)
% HPCBRA: output for HPC, in particular frequency om 
%   [pars, max(abs(u1)), min(abs(u1)), om]
%
nu=p.nu/3; upde=u(1:nu); np=3*nu/p.nc.neq; 
om=u(3*nu+p.nc.nq+3); 
out=[u(p.nu+1:p.nu+p.naux); max(abs(upde(1:np))); min(abs(upde(1:np)));om]; 